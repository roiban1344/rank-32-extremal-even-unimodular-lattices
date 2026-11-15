"""Gluing F and Z along an anti-isometry of their discriminant groups F_p^s, for elements of odd prime order p.

Lemma: if F and Z have minimum at least 4, the glued lattice L = L(phi) has a root if and only if a class a
of F^# containing a vector of norm 2/p and a class b of Z^# containing a vector of norm 2 - 2/p satisfy
phi(a) = b. So A = {classes of norm 2/p} and B = {classes of norm 2 - 2/p} decide the extremality of L.
  Disc(G, p)                  coordinates on X^#/X = F_p^s, the form Q, and the classes of given norm
  nullspace(M, p)             a basis of the right kernel of M over F_p
  random_phi(DF, DZ, rng)     a random anti-isometry phi (Q_Z(phi a) = -Q_F(a)), or None
  build_L(F, Z, DF, DZ, T)    the Gram matrix of the glued lattice (checked: determinant 1, even), LLL-reduced
  F_lattices(p, d, s), Z_lattices(z, s)   the lists of F and Z of minimum >= 4 (outputs of C2), read from
                              src_dir (default OUT)
  TYPES3                      the types 3-(z, d)-s with a fixed lattice that remain after the linear program
Used by computation C2 of the first paper (code/produce/element_types/prime_order_gluing_sample.py) and
by later computations of the types of order 3.

奇素数位数 p の元について、F と Z を判別式の群 F_p^s の反等長写像で貼り合わせる。

補題：F と Z の最小が 4 以上なら、貼り合わせの格子 L = L(phi) がルートを持つのは、ノルム 2/p の
ベクトルを含む F^# の類 a と、ノルム 2 - 2/p のベクトルを含む Z^# の類 b で phi(a) = b となるものが
あるときに限る。よって A = {ノルム 2/p の類} と B = {ノルム 2 - 2/p の類} で L の極値性が決まる。
  Disc(G, p)                  X^#/X = F_p^s の座標、形式 Q、指定したノルムの類
  nullspace(M, p)             F_p 上の M の右カーネルの基底
  random_phi(DF, DZ, rng)     乱択の反等長写像 phi（Q_Z(phi a) = -Q_F(a)）。作れなければ None
  build_L(F, Z, DF, DZ, T)    貼り合わせの格子の Gram 行列（行列式 1・偶を確かめる）。LLL 簡約
  F_lattices(p, d, s), Z_lattices(z, s)   最小 >= 4 の F と Z の一覧（C2 の出力）。src_dir（既定は OUT）から読む
  TYPES3                      線形計画の後に残る、固定格子のある型 3-(z, d)-s
第 1 の論文の計算 C2（code/produce/element_types/prime_order_gluing_sample.py）と、位数 3 の型の
後の計算が使う。
"""
import json
from pathlib import Path

import numpy as np

from sage.all import ZZ, QQ, GF, matrix, vector

from .forms import short_vectors
from .io import asset

OUT = asset("computation-outputs/C2")  # default directory of the lists of F and Z / F と Z の一覧の既定のディレクトリ


class Disc:
    """Coordinates on X^#/X = F_p^s with the form Q, and the classes that contain vectors of a given norm.

    X^#/X = F_p^s の座標と形式 Q、指定したノルムのベクトルを含む類。"""

    def __init__(self, G, p):
        self.G, self.p = G, p
        n = G.nrows()
        self.H = (p * G.inverse()).change_ring(ZZ)
        E = matrix(GF(p), self.H).echelon_form()
        self.rows = [r for r in E.rows() if not r.is_zero()]
        self.s = len(self.rows)
        self.piv = [next(i for i in range(n) if r[i] != 0) for r in self.rows]
        self.Y = matrix(QQ, [[int(x) for x in r] for r in self.rows]) / p
        B = p * self.Y * G * self.Y.transpose()
        assert B.denominator() == 1
        self.B = np.array([[int(x) % p for x in r] for r in B.rows()], dtype=np.int64)
        self.inv2 = pow(2, -1, p)

    def Q(self, a):
        """The value of the quadratic form Q on the coordinate vector a.

        座標のベクトル a での 2 次形式 Q の値。"""
        a = np.asarray(a, dtype=np.int64)
        return int(a @ self.B @ a % self.p) * self.inv2 % self.p

    def norm_classes(self, hnorm):
        """The coordinates (a set of tuples) of the classes of the vectors of sqrt(p) X^# of norm hnorm.

        sqrt(p) X^# のノルム hnorm のベクトルの類の座標（タプルの集合）。"""
        p = self.p
        C = short_vectors(self.H, hnorm)
        if C.shape[0] == 0:
            return set()
        Hn = np.array([[int(x) for x in r] for r in self.H.rows()], dtype=np.int64)
        norms = np.einsum("ij,jk,ik->i", C, Hn, C)
        Y = (C[norms == hnorm] @ Hn) % p
        out = set()
        for row in Y[:, self.piv]:
            t = tuple(int(x) for x in row)
            out.add(t)
            out.add(tuple((-x) % p for x in t))
        out.discard(tuple([0] * self.s))
        return out


def nullspace(M, p):
    """A basis of the right kernel of the matrix M (rows) over F_p.

    行列 M（行）の F_p 上の右カーネルの基底。"""
    M = M.copy() % p
    r, c = M.shape if M.size else (0, M.shape[1])
    piv = []
    row = 0
    for col in range(c):
        pr = None
        for i in range(row, r):
            if M[i, col]:
                pr = i
                break
        if pr is None:
            continue
        M[[row, pr]] = M[[pr, row]]
        M[row] = M[row] * pow(int(M[row, col]), -1, p) % p
        for i in range(r):
            if i != row and M[i, col]:
                M[i] = (M[i] - M[i, col] * M[row]) % p
        piv.append(col)
        row += 1
        if row == r:
            break
    free = [j for j in range(c) if j not in piv]
    basis = []
    for f in free:
        v = np.zeros(c, dtype=np.int64)
        v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-M[i, f]) % p
        basis.append(v)
    return basis


def random_phi(DF, DZ, rng):
    """A random phi: F_p^s -> F_p^s with Q_Z(phi a) = -Q_F(a), as the matrix T (acting on row vectors from the right); None on failure.

    乱択の phi: F_p^s -> F_p^s、Q_Z(phi a) = -Q_F(a)。行列 T（行ベクトルに右から）を返す。無理なら None。"""
    p, s = DF.p, DF.s
    V, W = [], []
    for i in range(s):
        bF = nullspace((np.array(V, dtype=np.int64) @ DF.B) % p, p) if V else [np.eye(s, dtype=np.int64)[k] for k in range(s)]
        bZ = nullspace((np.array(W, dtype=np.int64) @ DZ.B) % p, p) if W else [np.eye(s, dtype=np.int64)[k] for k in range(s)]
        found = False
        for _ in range(200):
            v = sum(rng.randrange(p) * b for b in bF) % p
            q = DF.Q(v)
            if q == 0:
                continue
            for _ in range(200):
                w = sum(rng.randrange(p) * b for b in bZ) % p
                if DZ.Q(w) == (-q) % p:
                    V.append(v)
                    W.append(w)
                    found = True
                    break
            if found:
                break
        if not found:
            return None
    Vm = matrix(GF(p), [[int(x) for x in v] for v in V])
    Wm = matrix(GF(p), [[int(x) for x in w] for w in W])
    T = Vm.inverse() * Wm
    return np.array([[int(x) for x in r] for r in T.rows()], dtype=np.int64)


def build_L(F, Z, DF, DZ, T):
    """The Gram matrix of the lattice glued along T (checked: determinant 1 and even), in an LLL-reduced basis.

    T で貼り合わせた格子の Gram 行列（行列式 1・偶を確かめる）。LLL 簡約の基底で。"""
    p = DF.p
    d, n2 = F.nrows(), Z.nrows()
    gens = [[p if i == j else 0 for j in range(d + n2)] for i in range(d + n2)]
    for i in range(DF.s):
        y = DF.Y[i]
        w = sum((QQ(int(T[i, j])) * DZ.Y[j] for j in range(DF.s)), vector(QQ, [0] * n2))
        gens.append([int(p * t) for t in list(y) + list(w)])
    H = matrix(ZZ, gens).hermite_form(include_zero_rows=False)
    B = H.change_ring(QQ) / p
    block = matrix.block_diagonal([F, Z], subdivide=False)
    L = B * block * B.transpose()
    assert L.denominator() == 1
    L = L.change_ring(ZZ)
    assert L.det() == 1 and all(L[i, i] % 2 == 0 for i in range(L.nrows()))
    U = L.LLL_gram().transpose()
    return U * L * U.transpose()


def F_lattices(p, d, s, src_dir=None):
    """The lattices F of minimum >= 4 of the genus (p, d, s), as (Gram, |Aut|); None if the list is missing.

    種 (p, d, s) の最小 >= 4 の格子 F の一覧 (Gram, |Aut|)。一覧が無ければ None。"""
    path = (OUT if src_dir is None else Path(src_dir)) / f"genus-p{p}-n{d}-s{s}.json"
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    return [(matrix(ZZ, c["gram"]), int(c["aut"])) for c in data["classes"] if c["min"] >= 4]


def Z_lattices(z, s, src_dir=None):
    """The lattices Z (p = 3, rank z over Z[omega], determinant 3^s) of minimum >= 4, as (Gram, |U|); None if the list is missing.

    The rootless unimodular lattices for s = z, otherwise the lattices with this s from the overlattice chain.

    最小 >= 4 の格子 Z（p = 3、Z[omega] 上の階数 z、行列式 3^s）の一覧 (Gram, |U|)。一覧が無ければ None。
    s = z ならルート無しのユニモジュラー格子、そうでなければ上格子の鎖の、この s の格子。"""
    base = OUT if src_dir is None else Path(src_dir)
    if s == z:
        path = base / f"unimod-z{z}.json"
        if not path.exists():
            return None
        d = json.loads(path.read_text())
        return [(matrix(ZZ, c["G"]), int(c["aut"])) for c in d["classes"] if c["inv"][0] == 0]
    path = base / f"overlat-z{z}.json"
    if not path.exists():
        return None
    d = json.loads(path.read_text())
    for lv in d["levels"]:
        if lv["s"] == s:
            return [(matrix(ZZ, c["gram"]), int(c["aut"])) for c in lv["classes"]]
    return None


TYPES3 = [(6, 20, 6), (8, 16, 8), (9, 14, 7), (9, 14, 9), (10, 12, 6), (10, 12, 8), (10, 12, 10), (12, 8, 8), (13, 6, 5)]
