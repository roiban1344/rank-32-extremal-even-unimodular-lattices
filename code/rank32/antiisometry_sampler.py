"""Tables of 2-elementary discriminant forms, and random anti-isometries between two of them.

Used by the involution gluing (computation C4 of the first paper): the components F = L^sigma (the
fixed lattice) and K = L^{-sigma} (the anti-fixed lattice) of an involution sigma are 2-elementary,
and the gluing is an anti-isometry phi of their discriminant forms.
disc(G) tabulates D = G^#/G (an F_2-space of rank k; an element is an integer mask, bit i is the
coefficient of basis vector x_i): 2q mod 4 of every element, 2 mu (the least norm of the class,
capped at 8), and 4 b(x_i, x_j) mod 4.
random_antiisometry chooses the images y_i in D_K of the basis x_i of D_F one after the other,
uniformly among the candidates with 2q_K(y_i) = -2q_F(x_i) and 4b_K(y_i, y_j) = -4b_F(x_i, x_j)
(mod 4) that are linearly independent of the earlier ones. It returns the list of the numbers of
candidates as well (the weight of the sample in importance sampling); a dead end returns None.

2-初等的な判別式形式の表と、その 2 つの間の乱択の反等長。

対合の貼り合わせ（第 1 の論文の計算 C4）が使う。対合 sigma の成分 F = L^sigma（固定格子）と
K = L^{-sigma}（反固定格子）は 2-初等的で、貼り合わせは判別式形式の反等長 phi である。
disc(G) は D = G^#/G（階数 k の F_2 上の空間。元は整数の mask で、ビット i が基底 x_i の係数）の
表を作る：各元の 2q mod 4、2 mu（類の最小ノルム。8 で頭打ち）、4 b(x_i, x_j) mod 4。
random_antiisometry は、D_F の基底 x_i の像 y_i in D_K を順に選ぶ。候補は
2q_K(y_i) = -2q_F(x_i)、4b_K(y_i, y_j) = -4b_F(x_i, x_j)（mod 4）を満たし、前の像と一次独立な
ものであり、その中から一様に選ぶ。候補の数の列（重点抽出での標本の重み）も返す。行き止まりなら None。
"""
import numpy as np
from sage.all import ZZ, QQ, GF, matrix, vector

from . import forms


def disc(G):
    """Table of the discriminant group G^#/G (2-elementary): basis x_i (rational vectors), 2q mod 4 of all elements, 2 mu (capped at 8), 4b(x_i, x_j) mod 4.

    判別式群 G^#/G（2-初等）の表：基底 x_i（有理ベクトル）、全元の 2q mod 4、2 mu（8 で頭打ち）、4b(x_i, x_j) mod 4。"""
    d = G.nrows(); Xi = G.inverse()
    gens = []; span = matrix(GF(2), 0, d)
    for j in range(d):
        c = Xi.column(j); v2 = vector(GF(2), [int(2 * x) % 2 for x in c])
        if span.stack(v2).rank() > span.rank():
            gens.append(c); span = span.stack(v2)
    k = len(gens); n = 1 << k
    X = matrix(QQ, gens); Q = X * G * X.transpose()
    B4 = np.array([[int(4 * Q[i, j]) % 4 for j in range(k)] for i in range(k)], dtype=np.int64)
    q2g = np.array([int(2 * Q[i, i]) % 4 for i in range(k)], dtype=np.int64)
    masks = np.arange(n, dtype=np.int64)
    bits = ((masks[:, None] >> np.arange(k)) & 1).astype(np.int64)          # n x k
    # 2|sum x_i|^2 = sum 2|x_i|^2 + 4 sum_{i<j} (x_i, x_j)  (mod 4)
    q2 = (bits @ q2g + ((bits @ np.triu(B4, 1)) * bits).sum(1)) % 4
    # least norm of a class: short vectors a of 2 G^{-1} (2|x|^2 = a^T (2G^{-1}) a <= 7), sorted into classes
    # 類の最小ノルム：2 G^{-1} の短いベクトル a（2|x|^2 = a^T (2G^{-1}) a <= 7）を類に直す
    K2 = (2 * Xi).change_ring(ZZ)
    V = forms.short_vectors(K2, 7)
    mu2 = np.full(n, 8, dtype=np.int64); mu2[0] = 0
    if len(V):
        K2n = np.array([[int(x) for x in r] for r in K2.rows()], dtype=np.int64)
        nrm = np.einsum("ij,jk,ik->i", V, K2n, V)
        Y = (V @ K2n.T) % 2                                               # coordinates of 2x mod 2 / 2x の座標 mod 2
        S = np.array([[int(x) for x in r] for r in span.rows()], dtype=np.int64)   # k x d
        J = matrix(GF(2), span).pivots()
        Sinv = np.array([[int(x) for x in r] for r in matrix(GF(2), span).matrix_from_columns(J).inverse().rows()], dtype=np.int64)
        C = (Y[:, list(J)] @ Sinv) % 2                                    # coordinates of the class / 類の座標
        cm = C @ (1 << np.arange(k))
        np.minimum.at(mu2, cm, np.minimum(nrm, 8))
    return dict(k=k, q2=q2, mu2=mu2, B4=B4, bits=bits)


def random_antiisometry(F, K, rng):
    """Choose the images y_i in D_K of the basis x_i of F uniformly one after the other, with 2q_K(y_i) = -2q_F(x_i), 4b_K(y_i, y_j) = -4b_F(x_i, x_j) (mod 4) and linear independence.

    F の基底 x_i の像 y_i in D_K を順に一様に選ぶ。条件：2q_K(y_i) = -2q_F(x_i)、4b_K(y_i, y_j) = -4b_F(x_i, x_j)（mod 4）、一次独立。"""
    k = F["k"]; n = 1 << k
    qtarget = (-F["q2"][1 << np.arange(k)]) % 4
    ys = []; counts = []
    cand_ok = np.ones(n, dtype=bool); cand_ok[0] = False
    span = {0}
    for i in range(k):
        ok = cand_ok & (K["q2"] == qtarget[i])
        for j, yj in enumerate(ys):
            g = (K["B4"] @ K["bits"][yj]) % 4                                # 4b_K(e_l, y_j)
            bv = (K["bits"] @ g) % 4                                         # 4b_K(y, y_j)
            ok &= bv == (-F["B4"][i, j]) % 4
        ok[list(span)] = False
        idx = np.nonzero(ok)[0]
        counts.append(len(idx))
        if not len(idx):
            return None, counts
        y = int(idx[rng.randrange(len(idx))])
        ys.append(y); span |= {s ^ y for s in span}
    return ys, counts


def phi_table(ys, k):
    """The image of every element of D_F under the map with basis images ys.

    基底の像が ys である写像による、D_F の全元の像。"""
    n = 1 << k; m = np.arange(n, dtype=np.int64); t = np.zeros(n, dtype=np.int64)
    for i, y in enumerate(ys):
        t ^= ((m >> i) & 1) * y
    return t
