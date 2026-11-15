"""Theta series of even p-elementary lattices for the linear programs on types (p = 3, 5, 7).

For a type p-(z, d)-s, the fixed part F (rank d) and the cyclic part Z (rank z (p - 1)) are both even,
p-elementary and of determinant p^s. The theta series of such a lattice X of rank n (k = n/2) lies in
M_k(Gamma0(p), chi), where chi is the quadratic character fixed by the parity of s, and the Fricke
involution maps theta_X to c_s theta_{sqrt(p) X^#}, c_s = (a common constant) p^((k - s)/2).
This module builds the pairs (theta_i, g_i) from known lattices X_i (orthogonal sums of the blocks
below, s_i = s mod 2), g_i = p^((s - s_i)/2) theta_{sqrt(p) X_i^#}, until the theta_i span the space:
  blocks(p)               the blocks (name, Gram, Gram of sqrt(p) X^#, rank, s)
                            p = 3: A2 (2, 1), E6 (6, 1), E8 (8, 0), sqrt(3) E8 (8, 8), sqrt(3) E6^# (6, 5)
                            p = 5: A4 (4, 1), sqrt(5) A4^# (4, 3), E8 (8, 0), sqrt(5) E8 (8, 8)
                            p = 7: A6 (6, 1), sqrt(7) A6^# (6, 5), E8 (8, 0), sqrt(7) E8 (8, 8), [[2,1],[1,4]] (2, 1)
  theta_q(G, scale)       the theta series as a power series in q^((x,x)/2), precision P
  spaces_for(p, n, s, th) the pairs, the rank they reach and the dimension of the space
  consistent(pairs)       the linear relations among the theta_i also hold for the g_i
Conditions for correctness of a linear program built on these pairs: (i) the theta_i span
M_k(Gamma0(p), chi) (rank = dimension); (ii) the g_i satisfy the relations of the theta_i (a check of
the linearity of W and of the constants). Truncating the q-expansion only relaxes the constraints.
Used by computation C2 of the first paper (rank32.joint_theta_lp).

型の線形計画のための、偶 p-基本格子のテータ級数（p = 3, 5, 7）。

型 p-(z, d)-s では、固定部分 F（階数 d）と巡回部分 Z（階数 z (p - 1)）は、どちらも偶・p-基本・
行列式 p^s である。そういう階数 n（k = n/2）の格子 X のテータ級数は M_k(Gamma0(p), chi) にある
（chi は s の偶奇で決まる 2 次指標）。Fricke 対合は theta_X を c_s theta_{sqrt(p) X^#} に写す。
c_s = (共通の定数) p^((k - s)/2) である。
このモジュールは、既知の格子 X_i（下のブロックの直和、s_i = s mod 2）から組 (theta_i, g_i)、
g_i = p^((s - s_i)/2) theta_{sqrt(p) X_i^#} を、theta_i が空間を張るまで集める。
  blocks(p)               ブロック (名前, Gram, sqrt(p) X^# の Gram, 階数, s)
                            p = 3: A2 (2, 1)、E6 (6, 1)、E8 (8, 0)、sqrt(3) E8 (8, 8)、sqrt(3) E6^# (6, 5)
                            p = 5: A4 (4, 1)、sqrt(5) A4^# (4, 3)、E8 (8, 0)、sqrt(5) E8 (8, 8)
                            p = 7: A6 (6, 1)、sqrt(7) A6^# (6, 5)、E8 (8, 0)、sqrt(7) E8 (8, 8)、[[2,1],[1,4]] (2, 1)
  theta_q(G, scale)       テータ級数を q^((x,x)/2) の冪級数で（精度 P）
  spaces_for(p, n, s, th) 組、届いた階数、空間の次元
  consistent(pairs)       theta_i の一次関係を g_i も満たす
この組の上の線形計画が正しい条件：(i) theta_i が M_k(Gamma0(p), chi) を張る（階数 = 次元）。
(ii) g_i が theta_i の一次関係を満たす（W の線形性と定数の検算）。q 展開の打ち切りは制約を緩める
だけである。
第 1 の論文の計算 C2 が使う（rank32.joint_theta_lp）。
"""
import itertools

from sage.all import (matrix, ZZ, QQ, pari, PowerSeriesRing, CartanMatrix,
                      ModularForms, DirichletGroup, kronecker_character)

from . import density_sieve as sieve

P = 64  # precision of the q-expansions (coefficients q^0..q^(P-1)) / q 展開の精度（係数 q^0..q^(P-1)）
Rq = PowerSeriesRing(QQ, "q", default_prec=P)


def theta_q(G, scale):
    """The theta series, as a series in q^((x,x)/2), of the lattice with real Gram matrix G/scale (G even, integral; scale 1 or p).

    The real norm is N = (G-norm)/scale. The callers pass the Gram matrix of sqrt(p) X^# made integral, with scale = 1.

    実 Gram = G/scale（G は偶な整数行列、scale は 1 か p）の格子のテータ級数を、q^((x,x)/2) の級数で返す。
    実ノルムは N = (G-ノルム)/scale。呼び出し側は sqrt(p) X^# の Gram を整数にして scale = 1 で渡す。"""
    G = matrix(ZZ, G)
    maxg = 2 * (P - 1) * scale
    v = pari(G).qfrep(maxg // 2, 1)
    c = [0] * P
    c[0] = 1
    for j in range(1, maxg // 2 + 1):
        realhalf = QQ(j) / scale  # real norm / 2 / 実ノルム/2
        if realhalf.denominator() == 1 and realhalf < P:
            c[int(realhalf)] += 2 * int(v[j - 1])
    return Rq(c)


def blocks(p):
    """The blocks for p: a list of (name, real Gram (integral), real Gram of sqrt(p) X^# (integral), rank, s).

    p のブロック：(名前, 実 Gram（整数）, sqrt(p) X^# の実 Gram（整数）, 階数, s) の列。"""
    E8 = matrix(ZZ, CartanMatrix(["E", 8]))
    out = []  # (name, real Gram, real Gram of sqrt(p) X^#, rank, s) / (名前, 実 Gram, sqrt(p) X^# の実 Gram, 次元, s)
    def add(name, G):
        G = matrix(ZZ, G)
        Gd = p * G.inverse()  # the real Gram matrix of sqrt(p) X^# is p G^-1 / sqrt(p) X^# の実 Gram = p G^-1
        assert Gd.denominator() == 1
        s = G.det().valuation(p)
        assert G.det() == p ** s
        out.append((name, G, matrix(ZZ, Gd), G.nrows(), s))
    add("E8", E8)
    add(f"rt{p}E8", p * E8)
    if p == 3:
        A2 = matrix(ZZ, [[2, -1], [-1, 2]])
        E6 = matrix(ZZ, CartanMatrix(["E", 6]))
        add("A2", A2)
        add("E6", E6)
        add("rt3E6d", 3 * E6.inverse())
    if p == 5:
        A4 = matrix(ZZ, CartanMatrix(["A", 4]))
        add("A4", A4)
        add("rt5A4d", 5 * A4.inverse())
    if p == 7:
        A6 = matrix(ZZ, CartanMatrix(["A", 6]))
        add("A6", A6)
        add("rt7A6d", 7 * A6.inverse())
        add("Q7", matrix(ZZ, [[2, 1], [1, 4]]))
    for b in out:
        G, Gd = b[1], b[2]
        assert all(G[i, i] % 2 == 0 for i in range(G.nrows())) and all(Gd[i, i] % 2 == 0 for i in range(Gd.nrows())), b[0]
    return out


def char_for(p, k, s):
    """The character of the space of theta: d -> ((-1)^k p^s / d); off the d divisible by p it is the quadratic character mod p or trivial.

    None when it is not a character mod p (then there is no genus).

    theta が入る空間の指標：d -> ((-1)^k p^s / d)。p で割れる d 以外では、mod p の 2 次指標か自明になる。
    mod p の指標にならないときは None（そのときは種が無い）。"""
    D = (-1) ** k * p ** (s % 2)
    if D == 1:
        return DirichletGroup(p)[0]
    eps = kronecker_character(D)
    G = DirichletGroup(p)
    for e in G:
        if all(e(a) == eps(a) for a in range(1, 4 * p) if a % p):
            return e
    return None


def spaces_for(p, n, s, th):
    """Collect the pairs (theta_i, g_i) of lattices of rank n and of the parity of s until the rank is full; returns (pairs, rank, dim).

    階数 n・s の偶奇の格子の (theta_i, g_i) の組を、階数が満ちるまで集める。(pairs, rank, dim) を返す。"""
    k = n // 2
    chi = char_for(p, k, s)
    if chi is None:
        return None, None, None
    dim = ModularForms(chi, k).dimension()
    names = list(th)
    pairs = []
    A_rows = []
    rank = 0
    # the numbers of blocks, in increasing order (rank n, s_i = s mod 2) / ブロックの個数の組を小さい順に（次元 n、s_i = s mod 2）
    ranges = [range(n // th[nm][2] + 1) for nm in names]
    for counts in itertools.product(*ranges):
        if sum(c * th[nm][2] for c, nm in zip(counts, names)) != n:
            continue
        si = sum(c * th[nm][3] for c, nm in zip(counts, names))
        if (si - s) % 2:
            continue
        f, fd = Rq(1), Rq(1)
        for c, nm in zip(counts, names):
            if c:
                f *= th[nm][0] ** c
                fd *= th[nm][1] ** c
        g = QQ(p) ** QQ((s - si) / 2) * fd if (s - si) % 2 == 0 else None
        if g is None:
            continue
        row = [f[j] for j in range(P)]
        newA = matrix(QQ, A_rows + [row])
        r = newA.rank()
        if r > rank:
            A_rows.append(row)
            pairs.append((f, g, dict(zip(names, counts)), si))
            rank = r
        else:
            pairs.append((f, g, dict(zip(names, counts)), si))  # the dependent ones are kept for the check / 検算用に従属なものも持つ
        if rank == dim and len(pairs) >= dim + 10:
            break
    return pairs, rank, dim


def consistent(pairs):
    """True when the g_i satisfy the linear relations of the theta_i (a check of the linearity of W and of the constants).

    theta_i の一次関係を g_i も満たすなら True（W の線形性と定数の検算）。"""
    A = matrix(QQ, [[f[j] for j in range(P)] for f, _, _, _ in pairs])
    B = matrix(QQ, [[g[j] for j in range(P)] for _, g, _, _ in pairs])
    return A.rank() == A.augment(B).rank()
