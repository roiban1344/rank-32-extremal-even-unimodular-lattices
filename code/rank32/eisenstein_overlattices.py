"""Overlattices of Z[omega]-lattices along isotropic lines of the discriminant group (p = 3).

One step Z subset Z': Z^#/Z = F_3^s, on which omega acts trivially, so every isotropic line [y]
((y, y) = 0 mod 2) gives a Z[omega]-lattice Z' = Z + Z y with s' = s - 2. The minimum of Z' is at least 4
if and only if the class of y contains no vector of norm 2 (no other norm below 4 can occur).
  Disc(G)                    coordinates on Z^#/Z = F_3^s and the action of isometries on it
  good_lines(G)              (Disc, all isotropic lines, the lines whose class has no vector of norm 2)
  orbits(D, gens, keys)      the orbits of a group on a set of lines (representative, size)
  overlattice(G, W, D, key)  the overlattice Z + Z y for the line key, with omega, in an LLL basis
  qfsafe                     rank32.forms with retries after a SystemError of cypari2 (a race when the alarm
                             fires just before sig_on); a PARI error (stack overflow) becomes QfTimeout
Lines are stored as integers: the coordinates in base 3, normalised so that the first nonzero one is 1.
The PARI stack is set to 1 GB (at most 2 GB) on import: a qfisom that does not finish keeps using it
(it reached 5 GB); when it overflows, the caller turns to another method.
Used by computation C2 of the first paper (code/produce/element_types/eisenstein_overlattice_chain.py)
and by later computations.

判別式の群の等方直線に沿った、Z[omega]-格子の上格子（p = 3）。

Z subset Z' の 1 歩：Z^#/Z = F_3^s で omega は自明に働くので、等方直線 [y]（(y, y) = 0 mod 2）は
どれも Z[omega]-格子 Z' = Z + Z y を与え、s' = s - 2 となる。Z' の最小が 4 以上であることは、y の
剰余類にノルム 2 のベクトルが無いことと同値である（4 未満のノルムは他に起きない）。
  Disc(G)                    Z^#/Z = F_3^s の座標と、その上の等長写像の作用
  good_lines(G)              (Disc, 全ての等方直線, 類にノルム 2 のベクトルを持たない直線)
  orbits(D, gens, keys)      群の直線の集合への作用の軌道（代表, 大きさ）
  overlattice(G, W, D, key)  直線 key の上格子 Z + Z y を、omega と一緒に LLL の基底で
  qfsafe                     rank32.forms に、cypari2 の SystemError（alarm が sig_on の直前に切れた
                             ときの競合）の後の再試行を足したもの。PARI のエラー（スタック溢れ）は
                             QfTimeout にする
直線は整数で持つ：座標を 3 進で並べ、0 でない最初の座標が 1 になるように正規化する。
import 時に PARI のスタックを 1 GB（最大 2 GB）にする。終わらない qfisom はスタックを使い続ける
（5 GB に達した）。溢れたら、呼び出し側が別の手段に移る。
第 1 の論文の計算 C2（code/produce/element_types/eisenstein_overlattice_chain.py）と、後の計算が使う。
"""
import numpy as np

from sage.all import ZZ, QQ, GF, matrix, vector, pari

from cypari2.handle_error import PariError

from . import forms as _qfsafe

pari.allocatemem(1 << 30, 2 << 30)  # a qfisom that does not finish keeps using the stack / 終わらない qfisom がスタックを使い続ける


class qfsafe:
    """rank32.forms with retries after a SystemError of cypari2 (the race when the alarm fires just before sig_on).

    rank32.forms に、cypari2 の SystemError（alarm が sig_on の直前に切れたときの競合）の再試行を足したもの。"""
    QfTimeout = _qfsafe.QfTimeout
    short_vectors = staticmethod(_qfsafe.short_vectors)

    @staticmethod
    def _again(f, *a, **k):
        for seed in range(1, 9):
            try:
                return f(*a, seed=seed, **k)
            except SystemError as ex:
                print("SystemError of cypari2, retrying:", str(ex)[:60], flush=True)
            except PariError as ex:
                print("PARI error (stack overflow or similar):", str(ex)[:80], flush=True)
                raise _qfsafe.QfTimeout("PARI stack overflow")
        raise _qfsafe.QfTimeout("SystemError repeated")

    @staticmethod
    def qfauto(*a, **k):
        return qfsafe._again(_qfsafe.qfauto, *a, **k)

    @staticmethod
    def qfisom(*a, **k):
        return qfsafe._again(_qfsafe.qfisom, *a, **k)


F3 = GF(3)


def forms(G, W):
    """The pair of forms [G, G W^T] whose common isometries are the Z[omega]-linear isometries.

    共通の等長写像が Z[omega]-線形な等長写像になる形式の組 [G, G W^T]。"""
    return [G, G * W.transpose()]


class Disc:
    """Coordinates on Z^#/Z = F_3^s.

    Z^#/Z = F_3^s の座標。"""

    def __init__(self, G):
        self.G = G
        n = G.nrows()
        self.Ginv = G.inverse()
        D3 = (3 * self.Ginv).change_ring(ZZ)
        assert all(x % 1 == 0 for x in D3.list())
        M = matrix(F3, D3)
        E = M.echelon_form()
        self.rows = [r for r in E.rows() if not r.is_zero()]
        self.s = len(self.rows)
        self.piv = [next(i for i in range(n) if r[i] != 0) for r in self.rows]
        self.E = matrix(F3, self.rows)
        # lifts of the basis y_i = (lift of r_i)/3 / 基底の持ち上げ y_i = (r_i のリフト)/3
        Y = matrix(QQ, [[int(x) for x in r] for r in self.rows]) / 3
        self.Y = Y
        B = 3 * Y * G * Y.transpose()
        assert B.denominator() == 1
        self.B = np.array(matrix(ZZ, B).apply_map(lambda x: int(x) % 3).list(), dtype=np.int64).reshape(self.s, self.s)

    def coords(self, y3):
        """The coordinates (in F_3^s) of y = y3/3, for an integer row vector y3.

        y = y3/3（整数の行ベクトル y3）の座標（F_3^s）。"""
        a = [int(y3[p]) % 3 for p in self.piv]
        return a

    def action(self, T):
        """The matrix on F_3^s of x -> x T^t (acting on row vectors from the right).

        x -> x T^t の F_3^s 上の行列（行ベクトルに右から）。"""
        A = []
        for i in range(self.s):
            y3 = (3 * self.Y[i] * T.transpose())
            assert y3.denominator() == 1
            A.append(self.coords(y3.change_ring(ZZ)))
        return np.array(A, dtype=np.int64)


def line_key(a):
    """The integer key of the line spanned by the vector a of F_3^s.

    F_3^s のベクトル a が張る直線の整数のキー。"""
    for x in a:
        if x:
            if x == 2:
                a = [(2 * t) % 3 for t in a]
            break
    k = 0
    for x in a:
        k = 3 * k + int(x)
    return k


def key_to_vec(k, s):
    """The normalised vector of F_3^s of a line key.

    直線のキーの、正規化した F_3^s のベクトル。"""
    v = [0] * s
    for i in range(s - 1, -1, -1):
        v[i] = k % 3
        k //= 3
    return v


def good_lines(G):
    """All isotropic lines and those that meet no vector of norm 2, as sets of keys: returns (Disc, all, good).

    等方直線の全体と、ノルム 2 のベクトルに当たらないものの集合（キーの集合）。(Disc, 全体, 良い直線) を返す。"""
    D = Disc(G)
    s = D.s
    # all isotropic lines / 全等方直線
    if s == 0:
        return D, set(), set()
    idx = np.arange(3 ** s, dtype=np.int64)
    vecs = np.zeros((3 ** s, s), dtype=np.int8)
    t = idx.copy()
    for i in range(s - 1, -1, -1):
        vecs[:, i] = t % 3
        t //= 3
    quad = np.einsum("ij,jk,ik->i", vecs.astype(np.int64), D.B, vecs.astype(np.int64)) % 3
    first = np.argmax(vecs != 0, axis=1)
    lead = vecs[np.arange(3 ** s), first]
    mask = (quad == 0) & (lead == 1)
    all_iso = set(int(x) for x in idx[mask])
    del vecs, quad, first, lead, mask, idx, t
    # the vectors of norm 2 (norm 6 after scaling by sqrt(3)) / ノルム 2 のベクトル（sqrt(3) 倍してノルム 6）
    H = (3 * D.Ginv).change_ring(ZZ)
    half = qfsafe.short_vectors(H, 6)
    bad = set()
    Hnp = np.array([[int(x) for x in r] for r in H.rows()], dtype=np.int64)
    if half.shape[0]:
        norms = np.einsum("ij,jk,ik->i", half, Hnp, half)
        Y3 = (half[norms == 6] @ Hnp) % 3
        for row in Y3:
            a = [int(row[p]) for p in D.piv]
            if not any(a):
                continue  # a vector of Z (does not occur since min Z >= 4) / Z の中のベクトル（min Z >= 4 なので起きない）
            bad.add(line_key(a))
    return D, all_iso, all_iso - bad


def orbits(D, gens, keys):
    """The orbits of the group generated by gens on the set of line keys, as a list of (representative, size).

    gens の生成する群の、直線のキーの集合への作用の軌道。(代表, 大きさ) の列。"""
    keys = set(keys)
    acts = [D.action(T) for T in gens]
    reps = []
    while keys:
        k0 = keys.pop()
        reps.append(k0)
        stack = [k0]
        size = 1
        while stack:
            k = stack.pop()
            v = np.array(key_to_vec(k, D.s), dtype=np.int64)
            for A in acts:
                w = (v @ A) % 3
                kk = line_key(list(w))
                if kk in keys:
                    keys.discard(kk)
                    stack.append(kk)
                    size += 1
        reps[-1] = (k0, size)
    return reps


def overlattice(G, W, D, key):
    """The overlattice Z + Z y of the line key, with the matrix of omega, in an LLL-reduced basis: returns (Gram, W).

    直線 key の上格子 Z + Z y を、omega の行列と一緒に LLL 簡約の基底で返す：(Gram, W)。"""
    a = key_to_vec(key, D.s)
    y = sum((QQ(int(c)) * D.Y[i] for i, c in enumerate(a)), vector(QQ, [0] * G.nrows()))
    n = G.nrows()
    gens = matrix(ZZ, [[3 if i == j else 0 for j in range(n)] for i in range(n)] + [list(3 * y)])
    H = gens.hermite_form(include_zero_rows=False)
    assert H.nrows() == n
    B = H.change_ring(QQ) / 3
    Gn = B * G * B.transpose()
    assert Gn.denominator() == 1
    Gn = Gn.change_ring(ZZ)
    assert all(Gn[i, i] % 2 == 0 for i in range(n))
    Wn = B * W * B.inverse()
    assert Wn.denominator() == 1, "not invariant under omega"
    Wn = Wn.change_ring(ZZ)
    T = Gn.LLL_gram().transpose()
    Tinv = T.inverse()
    Gl = T * Gn * T.transpose()
    Wl = (T * Wn * Tinv).change_ring(ZZ)
    assert Wl * Gl * Wl.transpose() == Gl
    return Gl, Wl
