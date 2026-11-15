"""Invariants, canonical forms and unitary groups of Gaussian Hermitian lattices of rank 16 (the keys of the enumeration).

Paper: computation C12 of the first paper (the enumeration of the rootless Gaussian lattices);
also used by C1 and by the check of the shell counts. The census of the walk and the check of |U|
against centralizers are code/produce/gaussian/gaussian_census.py and
code/crosscheck/gaussian/unitary_order_vs_centralizer.py.
L is a Hermitian lattice of rank 16 over Z[i] (g is x -> x W), h(x, y) = (x, y) + i (x, i y), pi = 1 + i.
  V = L/piL is F_2^16 and omega(x, y) = (x, y + iy) mod 2 is a nondegenerate alternating form.
  n(c) = number of Z[i]-lines of norm 4 in the class c, r(c) = number of lines of norm 2.
  From n and r we build U-equivariant functions on V by Walsh-transform convolutions
  (n*n, n*n*n, r*n); the colour of a class is the tuple of their values, and the number of
  classes of each colour is an invariant (digests A, C, D).
Canonical form (rootless): take the lines of norm 4 in classes with n = 1 whose colour occurs once,
in increasing colour, until they generate L, and write the matrix of h normalised by units.
  If they generate L, equal matrices prove U-isometry, different ones prove non-isometry, and
  |U| = 4 follows (each line is fixed and the unit is fixed by h). This is the status DISCRETE.
  Otherwise a backtracking search (matching h among the lines of the same colour) gives |U| and isometry.
One step of the pi-neighbour walk: v in L not in piL, y = pi v, M = {x : (x, y) even}, N = M + Z y/2.
The C program _fincke_pohst.c is compiled on first use (gcc -O3 -march=native) into $RANK32_BUILD
(default ~/.cache/rank32).

Z[i] 上の階数 16 のエルミート格子の不変量・標準形・ユニタリ群（列挙のキー）。

論文：第 1 の論文の計算 C12（ルート無しの Gauss の格子の列挙）。C1 と、殻の本数の検算でも使う。
ランダムウォークの調査と、|U| と中心化群の照合は、code/produce/gaussian/gaussian_census.py と
code/crosscheck/gaussian/unitary_order_vs_centralizer.py にある。
L は Z[i] 上の階数 16 のエルミート格子（g は x -> x W）。h(x, y) = (x, y) + i (x, i y)、pi = 1 + i。
  V = L/piL は F_2^16 で、omega(x, y) = (x, y + iy) mod 2 は非退化な交代形式である。
  n(c) = 類 c に入るノルム 4 の Z[i] 直線の本数、r(c) = ノルム 2 の直線の本数。
  n、r から Walsh 変換の畳み込みで作る V 上の関数（n*n、n*n*n、r*n）は U-同変である。類の色を
  それらの値の組とし、色ごとの類の個数を不変量とする（ダイジェスト A、C、D）。
標準形（ルート無し）：色が一度しか現れない n = 1 の類の、ノルム 4 の直線を色の順に並べ、L を生成するまで
取る。単位で正規化した h の行列を書く。
  生成するなら、行列の一致は U-同型の証明、不一致は非同型の証明で、|U| = 4 も従う（各直線が固定され、
  単位は h で決まる）。これが状態 DISCRETE である。
  生成しないときは、後戻りの探索（同じ色の直線の中で h を合わせる）で |U| と同型を決める。
pi 近傍のランダムウォークの 1 歩：v in L で v not in piL、y = pi v、M = {x : (x, y) は偶数}、N = M + Z y/2。
C のプログラム _fincke_pohst.c は、初めて使うときに $RANK32_BUILD（既定は ~/.cache/rank32）へ
コンパイルする（gcc -O3 -march=native）。
"""
import ctypes
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

_SRC = Path(__file__).with_name("_fincke_pohst.c")
LINE_KEY = np.random.default_rng(20261031).integers(-(1 << 40), 1 << 40, size=32, dtype=np.int64)
CAP = 1 << 21
_lib = None
_buf = None

# Values of the field "status" of analyse() / analyse() の欄 status の値
DISCRETE = "discrete"            # the canonical form exists (|U| = 4) / 標準形がある（|U| = 4）
NON_DISCRETE = "non-discrete"    # fewer than 16 lines with a unique colour / 色が唯一の直線が 16 本未満
NOT_GENERATING = "not-generating"  # those lines do not generate L / それらの直線が L を生成しない
HAS_ROOTS = "has-roots"          # the lattice has roots; no canonical form is tried / ルートがある。標準形は作らない


def _build_dir():
    """The directory for compiled helpers.

    コンパイルしたものを置くディレクトリ。"""
    d = os.environ.get("RANK32_BUILD")
    if d:
        return Path(d)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "rank32"


def enum_short(G):
    """The vectors of norm <= 4, one of each pair +-v (Fincke-Pohst in C; the same sets as PARI qfminim on 3 lattices).

    ノルム <= 4 のベクトルを +-v の片方ずつ返す（C の Fincke-Pohst。3 つの格子で PARI の qfminim と集合が一致した）。"""
    global _lib, _buf
    if _lib is None:
        build = _build_dir()
        so = build / ("fincke_pohst-" + hashlib.sha256(_SRC.read_bytes()).hexdigest()[:12] + ".so")
        if not so.exists():
            build.mkdir(parents=True, exist_ok=True)
            tmp = build / f"fincke_pohst-{os.getpid()}.so.tmp"
            subprocess.run(["gcc", "-O3", "-march=native", "-shared", "-fPIC", "-o", str(tmp), str(_SRC), "-lm"], check=True)
            os.rename(tmp, so)
        _lib = ctypes.CDLL(str(so))
        _lib.fp_enum.restype = ctypes.c_longlong
        _buf = np.zeros((CAP, 32), dtype=np.int64)
    g = np.array(G, dtype=np.int64)
    nodes = ctypes.c_longlong(0)
    c = _lib.fp_enum(ctypes.c_int(32), g.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(4), _buf.ctypes.data_as(ctypes.c_void_p), ctypes.c_longlong(CAP), ctypes.byref(nodes))
    if c < 0:
        raise ArithmeticError(f"enumeration failed {c}")
    return _buf[:c].copy()


def f2_right_kernel(A):
    """A basis (rows) of the right kernel of the matrix A over F_2.

    F_2 上の行列 A の右カーネルの基底（行）。"""
    A = A.copy() % 2
    m, n = A.shape
    piv = []
    r = 0
    for j in range(n):
        p = next((i for i in range(r, m) if A[i, j]), None)
        if p is None:
            continue
        A[[r, p]] = A[[p, r]]
        for i in range(m):
            if i != r and A[i, j]:
                A[i] ^= A[r]
        piv.append(j)
        r += 1
        if r == m:
            break
    free = [j for j in range(n) if j not in piv]
    K = np.zeros((len(free), n), dtype=np.int64)
    for t, j in enumerate(free):
        K[t, j] = 1
        for i, q in enumerate(piv):
            K[t, q] = A[i, j]
    return K


def fwht(a):
    a = a.copy()
    h = 1
    n = len(a)
    while h < n:
        a = a.reshape(-1, 2 * h)
        x, y = a[:, :h].copy(), a[:, h:].copy()
        a[:, :h], a[:, h:] = x + y, x - y
        a = a.reshape(-1)
        h *= 2
    return a


class Structure:
    """The lines of one lattice and the structure of V = L/piL.

    1 つの格子の直線と、V = L/piL の構造。"""

    def __init__(self, G, W):
        self.G, self.W = np.array(G, dtype=np.int64), np.array(W, dtype=np.int64)
        V = enum_short(self.G)
        norms = np.einsum("ij,ij->i", V @ self.G, V)
        self.root_reps = self._lines(V[norms == 2])
        self.reps = self._lines(V[norms == 4])
        self.roots = 4 * len(self.root_reps)
        I = np.eye(32, dtype=np.int64)
        self.K = f2_right_kernel((I + self.W) % 2)  # 16 rows; x -> K x are the coordinates on L/piL / 16 行。x -> K x が L/piL の座標
        assert self.K.shape == (16, 32)
        self.pow2 = 1 << np.arange(16, dtype=np.int64)
        self.cls = self.classify(self.reps)
        self.root_cls = self.classify(self.root_reps)

    def _lines(self, V):
        if len(V) == 0:
            return V
        k, kk = V @ LINE_KEY, (V @ self.W) @ LINE_KEY
        assert not (np.abs(k) == np.abs(kk)).any()
        R = V[np.abs(k) < np.abs(kk)]
        assert 2 * len(R) == len(V)
        return R

    def classify(self, X):
        if len(X) == 0:
            return np.zeros(0, dtype=np.int64)
        return (((X @ self.K.T) % 2) * self.pow2).sum(axis=1)

    def omega_table(self):
        """The table u -> Omega^T u (as indices), with omega(u, c) = (Omega^T u) . c.

        u -> Omega^T u の添字の表。omega(u, c) = (Omega^T u) . c。"""
        A = np.concatenate([self.K.copy(), np.eye(16, dtype=np.int64)], axis=1) % 2
        piv = []
        r = 0
        for j in range(32):
            p = next((i for i in range(r, 16) if A[i, j]), None)
            if p is None:
                continue
            A[[r, p]] = A[[p, r]]
            for i in range(16):
                if i != r and A[i, j]:
                    A[i] ^= A[r]
            piv.append(j)
            r += 1
            if r == 16:
                break
        lift = np.zeros((16, 32), dtype=np.int64)
        for e in range(16):
            for i, q in enumerate(piv):
                lift[e, q] = A[i, 32 + e]
        assert ((lift @ self.K.T) % 2 == np.eye(16, dtype=np.int64)).all()
        Om = (lift @ self.G @ ((lift + lift @ self.W).T)) % 2  # Omega_{jk} = omega(e_j, e_k)
        assert (Om == Om.T).all() and not np.diag(Om).any()
        table = np.zeros(1 << 16, dtype=np.int64)
        rows = (Om * self.pow2).sum(axis=1)
        for j in range(16):
            s = 1 << j
            table[s:2 * s] = table[0:s] ^ rows[j]
        return table

    def herm(self, X, Y):
        """The matrices of the real and imaginary parts of h(x, y).

        h(x, y) の実部と虚部の行列。"""
        return X @ self.G @ Y.T, X @ self.G @ (Y @ self.W).T


def functions(S):
    """The functions on V. The transform of n is returned only for a check (without roots it is a function of n).

    V 上の関数。n の変換は検算のためだけに返す（ルート無しでは n の関数になる）。"""
    n = np.bincount(S.cls, minlength=1 << 16).astype(np.int64)
    r = np.bincount(S.root_cls, minlength=1 << 16).astype(np.int64)
    F = fwht(n)
    hat = F[S.omega_table()]
    nn = fwht(F * F) >> 16
    nnn = fwht(F * F * F) >> 16
    if len(S.root_cls):
        rn = fwht(fwht(r) * F) >> 16
    else:
        rn = np.zeros(1 << 16, dtype=np.int64)
    return n, r, hat, nn, nnn, rn


def mix(*arrs):
    h = np.zeros(1 << 16, dtype=np.uint64)
    for a in arrs:
        h = (h ^ a.astype(np.uint64)) * np.uint64(0x9E3779B97F4A7C15)
        h ^= h >> np.uint64(29)
    return h


def digest(*arrs):
    """Hash of the multiset of tuples of colours.

    色の組の多重集合のハッシュ。"""
    u, k = np.unique(mix(*arrs), return_counts=True)
    return hashlib.sha256(u.tobytes() + k.tobytes()).hexdigest()[:16], len(u)


def normalise(S, X, sign):
    """Rotate the lines of X by units and return the matrix of h (real part, imaginary part); None if they are not connected.
    sign = -1 gives the complex conjugate form.

    X の直線を単位で回し、h の行列（実部, 虚部）を返す。つながらなければ None。sign = -1 は複素共役の形。"""
    t = len(X)
    a, b = S.herm(X, X)
    b = sign * b
    # x_j -> eps x_j in the order reached from x_0 through h != 0, smallest index first (canonical)
    # x_0 から h != 0 でたどる順に x_j -> eps x_j。添字の最も小さいものから（標準的）
    seen = [0]
    rest = list(range(1, t))
    while rest:
        j = k = None
        for c in rest:
            k = next((k for k in seen if a[k, c] or b[k, c]), None)
            if k is not None:
                j = c
                break
        if j is None:
            return None
        rest.remove(j)
        seen.append(j)
        x, y = int(a[k, j]), int(b[k, j])
        e = 0
        while not (x > 0 and y >= 0):
            x, y = y, -x  # multiply by i^{-1} / i^{-1} を掛ける
            e += 1
        for _ in range(e):  # column j times i^{-1}, row j times i / 第 j 列に i^{-1}、第 j 行に i を掛ける
            a[:, j], b[:, j] = b[:, j].copy(), -a[:, j].copy()
            a[j, :], b[j, :] = -b[j, :].copy(), a[j, :].copy()
    assert not np.diag(b).any()
    return a, b


def generates(S, X):
    """Whether the Z-span of the x and i x is L.

    x と ix の Z-スパンが L か。"""
    from sage.all import ZZ, matrix
    A = np.concatenate([X, X @ S.W], axis=0)
    M = matrix(ZZ, A.tolist())
    if M.rank() < 32:
        return False
    return M.hermite_form(include_zero_rows=False).det().abs() == 1


def canon(S, m, n):
    """The canonical form from the lines in classes with n = 1 whose colour occurs once (for rootless lattices).

    色が一度しか現れない n = 1 の類の直線から、標準形を作る（ルート無しのとき）。"""
    u, idx, k = np.unique(m, return_index=True, return_counts=True)
    single = idx[(k == 1) & (n[idx] == 1)]  # increasing colour / 色の昇順
    if len(single) < 16:
        return dict(status=NON_DISCRETE, singles=int(len(single)))
    pos = np.full(1 << 16, -1, dtype=np.int64)
    pos[S.cls] = np.arange(len(S.cls))
    X = S.reps[pos[single]]
    t = 16
    limit = min(len(X), 64)  # if they do not generate, |U| > 4 or there is an antilinear element; do not test ranks for long
    # 生成しないときは |U| > 4 か反線形の元がある。行列の階数の検査を長く続けない
    while t <= limit and not generates(S, X[:t]):
        t += 4
    if t > limit:
        return dict(status=NOT_GENERATING, singles=int(len(single)))
    plus, minus = normalise(S, X[:t], 1), normalise(S, X[:t], -1)
    assert plus is not None and minus is not None  # L is not an orthogonal sum / L は直和に分かれない

    def h(ab):
        return hashlib.sha256(ab[0].astype(np.int16).tobytes() + ab[1].astype(np.int16).tobytes()).hexdigest()[:24]

    cp, cm = h(plus), h(minus)
    return dict(status=DISCRETE, singles=int(len(single)), t=t, plus=cp, minus=cm, key=min(cp, cm), self_conj=cp == cm)


def analyse(G, W):
    """Invariants and canonical form of one lattice (Gram G, g: x -> x W).

    1 つの格子（Gram G、g は x -> x W）の不変量と標準形。"""
    t0 = time.process_time()
    S = Structure(G, W)
    t1 = time.process_time()
    n, r, hat, nn, nnn, rn = functions(S)
    selfdual = bool((hat[1:] == 256 * n[1:] - 144).all()) and int(hat[0]) == len(S.reps) if S.roots == 0 else None
    t2 = time.process_time()
    dA = digest(n, r)
    dC = digest(n, r, nn, rn)
    dD = digest(n, r, nn, nnn, rn)
    t3 = time.process_time()
    c = canon(S, mix(n, nn, nnn), n) if S.roots == 0 else dict(status=HAS_ROOTS)
    t4 = time.process_time()
    return dict(roots=S.roots, lines=len(S.reps), selfdual=selfdual, hist=np.bincount(n).tolist(), rhist=np.bincount(r).tolist(), A=dA[0], C=dC[0], D=dD[0], nA=dA[1], nC=dC[1], nD=dD[1], **c,
                sec=dict(short=round(t1 - t0, 3), walsh=round(t2 - t1, 4), digest=round(t3 - t2, 4), canon=round(t4 - t3, 3)))


class Limit(Exception):
    pass


class Search:
    """Search for the images of the lines of base, restricted by colours and by h.

    base の直線の像を、色と h で絞って探す。"""

    def __init__(self, S):
        self.S = S
        n, r, _, nn, nnn, rn = functions(S)
        m = mix(n, r, nn, nnn, rn)
        self.X = np.concatenate([S.root_reps, S.reps], axis=0)
        cls = np.concatenate([S.root_cls, S.cls])
        norm = np.concatenate([np.full(len(S.root_reps), 2, dtype=np.uint64), np.full(len(S.reps), 4, dtype=np.uint64)])
        self.col = m[cls] ^ (norm * np.uint64(0xC2B2AE3D27D4EB4F))
        self.XW = self.X @ S.W

    def cells(self):
        u, inv, k = np.unique(self.col, return_inverse=True, return_counts=True)
        return inv, k

    def base(self):
        """Choose lines, smallest colour classes first, to get a connected sequence that generates L.

        小さい色の組から順に直線を選び、L を生成する、つながった列を作る。"""
        from sage.all import ZZ, QQ, matrix
        S = self.S
        inv, k = self.cells()
        order = np.lexsort((self.col, k[inv]))
        chosen = []
        p = 2147483629
        ech = []  # (pivot, row) mod p / (ピボット, 行) mod p

        def reduce(v):
            v = [int(x) % p for x in v]
            for piv, row in ech:
                if v[piv]:
                    c = v[piv]
                    v = [(a - c * b) % p for a, b in zip(v, row)]
            return v

        def add(v):
            v = reduce(v)
            piv = next((i for i, a in enumerate(v) if a), None)
            if piv is None:
                return False
            inv_ = pow(v[piv], p - 2, p)
            ech.append((piv, [a * inv_ % p for a in v]))
            return True

        pending = list(order)
        # first up to rank 32 / まず階数 32 まで
        while len(ech) < 32:
            progress = False
            for pos, i in enumerate(pending):
                x = self.X[i]
                if chosen and not (self.X[chosen] @ S.G @ x != 0).any() and not (self.X[chosen] @ S.G @ self.XW[i] != 0).any():
                    continue
                if not any(reduce(x)):
                    pending.pop(pos)
                    progress = True
                    break
                add(x)
                add(self.XW[i])
                chosen.append(i)
                pending.pop(pos)
                progress = True
                break
            if not progress:
                raise ArithmeticError("no generating sequence")
        # then add lines until the index is 1 / 指数が 1 になるまで足す
        M = matrix(ZZ, np.concatenate([self.X[chosen], self.XW[chosen]]).tolist()).hermite_form(include_zero_rows=False)
        for i in pending:
            if M.det().abs() == 1:
                break
            x = matrix(QQ, [list(map(int, self.X[i]))]) * M.change_ring(QQ).inverse()
            if x.denominator() == 1:
                continue
            chosen.append(i)
            M = matrix(ZZ, np.concatenate([self.X[chosen], self.XW[chosen]]).tolist()).hermite_form(include_zero_rows=False)
        assert M.det().abs() == 1
        return chosen


class Plan:
    """Realise the h of the base of source in target; conj = True uses the complex conjugate h (antilinear maps).

    source の base の h を target で実現する。conj = True なら複素共役の h（反線形の写像）。"""

    def __init__(self, source, target, conj=False, limit=200000):
        self.src, self.tgt = source, target
        self.basis = source.base()
        B = source.X[self.basis]
        self.A, self.B = source.S.herm(B, B)
        if conj:
            self.B = -self.B
        self.colours = source.col[self.basis]
        self.nodes = 0
        self.limit = limit
        self.cells = {}
        self.parent = []
        for j in range(len(self.basis)):
            self.parent.append(next((k for k in range(j) if self.A[k, j] or self.B[k, j]), None))
            assert j == 0 or self.parent[j] is not None

    def candidates(self, j, images):
        """Candidates for the j-th base vector (vectors, the unit included).

        j 番目の基底ベクトルの候補（単位まで決めたベクトル）。"""
        T = self.tgt
        c = self.colours[j]
        if c not in self.cells:
            self.cells[c] = np.nonzero(T.col == c)[0]
        idx = self.cells[c]
        if j == 0:
            return T.X[idx]
        gy = images[self.parent[j]] @ T.S.G
        a, b = T.X[idx] @ gy, T.XW[idx] @ gy
        ta, tb = int(self.A[self.parent[j], j]), int(self.B[self.parent[j], j])
        keep = np.nonzero(a * a + b * b == ta * ta + tb * tb)[0]
        out = []
        x, y = a[keep], b[keep]
        for e in range(4):  # h(y, i^e x) = i^{-e} h(y, x)
            hit = np.nonzero((x == ta) & (y == tb))[0]
            if len(hit):
                ii = idx[keep[hit]]
                v = {0: T.X[ii], 1: T.XW[ii], 2: -T.X[ii], 3: -T.XW[ii]}[e]
                out.append(v)
            x, y = y, -x
        if not out:
            return np.zeros((0, 32), dtype=np.int64)
        C = np.concatenate(out)
        G, W = T.S.G, T.S.W
        for k in range(j):
            if k == self.parent[j]:
                continue
            gy = images[k] @ G
            a = C @ gy
            b = (C @ W) @ gy
            f = (a == self.A[k, j]) & (b == self.B[k, j])
            C = C[f]
            if not len(C):
                break
        return C

    def extend(self, j, images):
        self.nodes += 1
        if self.nodes > self.limit:
            raise Limit
        if j == len(self.basis):
            return True
        for v in self.candidates(j, images):
            if self.extend(j + 1, images + [v]):
                return True
        return False

    def exists(self):
        return self.extend(0, [])

    def order(self):
        """|U| = 4 * product of the orbit sizes at each position (when source = target); position 0 is an orbit of lines.

        |U| = 4 * 各位置の軌道の大きさの積（source = target のとき）。位置 0 は直線の軌道。"""
        res = 4
        prefix = []
        B = self.src.X[self.basis]
        for j in range(len(self.basis)):
            orbit = 0
            for v in self.candidates(j, prefix):
                if self.extend(j + 1, prefix + [v]):
                    orbit += 1
            assert orbit >= 1
            res *= orbit
            prefix.append(B[j])
        return res


def unitary_order(G, W, limit=200000):
    """|U| by the backtracking search, and whether an antilinear automorphism exists.

    後戻りの探索による |U| と、反線形の自己同型があるか。"""
    s = Search(Structure(G, W))
    p = Plan(s, s, limit=limit)
    u = p.order()
    q = Plan(s, s, conj=True, limit=limit)
    return dict(U=u, antilinear=q.exists(), t=len(p.basis), nodes=p.nodes + q.nodes)


def verify_unitary(G, W, limit=2000000):
    """Check of |U|: the search without colours (split by norm only) gives the same order; the same base is used.

    |U| の検算：色を使わない（ノルムだけで分ける）探索でも同じ位数になるか。base は同じものを使う。"""
    s = Search(Structure(G, W))
    with_colour = Plan(s, s, limit=limit).order()
    basis = s.base()
    plain = Search(s.S)
    plain.base = lambda: basis  # the base from before the colours are dropped (the order does not depend on it)
    # 色を消す前の base を使う（base の取り方は位数を変えない）
    norm = np.concatenate([np.full(len(s.S.root_reps), 2, dtype=np.uint64), np.full(len(s.S.reps), 4, dtype=np.uint64)])
    plain.col = norm
    p = Plan(plain, plain, limit=limit)
    without = p.order()
    return dict(U=with_colour, U_plain=without, plain_nodes=p.nodes, agree=with_colour == without)


def iso(r1, r2, limit=200000):
    """Whether two rows (gram, W) are isometric by a linear or an antilinear map.

    2 つの行 (gram, W) が線形または反線形の写像で同型か。"""
    s1, s2 = Search(Structure(r1["gram"], r1["W"])), Search(Structure(r2["gram"], r2["W"]))
    return dict(linear=Plan(s1, s2, limit=limit).exists(), antilinear=Plan(s1, s2, conj=True, limit=limit).exists())


def read_jsonl(path):
    """Rows of a JSON Lines file ([] if it does not exist).

    JSON Lines のファイルの行（無ければ []）。"""
    if not Path(path).exists():
        return []
    with open(path) as f:
        return [json.loads(l) for l in f]


# ---------------------------------------------------------------- pi-neighbours (one step of the walk) / pi 近傍（ランダムウォークの 1 歩）

def reduce_(G, W):
    from sage.all import ZZ
    T = G.LLL_gram().transpose()
    Gl = T * G * T.transpose()
    Wl = (T * W * T.inverse()).change_ring(ZZ)
    return Gl, Wl


def neighbour(G, W, v):
    """The pi-neighbour defined by v; None if v is in pi L.

    v の定める pi 近傍。v が pi L に入るなら None。"""
    from sage.all import ZZ, QQ, matrix, vector
    n = G.nrows()
    y = v + v * W
    a = G * y
    k = next((i for i in range(n) if a[i] % 2), None)
    if k is None:
        return None
    if (v * G * v) % 4:
        # lift v to v + pi e_k with x = e_k such that (x, y) is odd; (v, v) changes by 2 (a + b)
        # (x, y) が奇数の x = e_k で v を v + pi e_k に持ち上げる。(v, v) は 2 (a + b) だけ変わる
        e = vector(ZZ, [1 if i == k else 0 for i in range(n)])
        v = v + e + e * W
        y = v + v * W
        a = G * y
        k = next(i for i in range(n) if a[i] % 2)
    assert (y * G * y) % 8 == 0
    rows = []
    for i in range(n):
        e = [0] * n
        if i == k:
            e[k] = 4
        else:
            e[i] = 2
            e[k] = 2 * (int(a[i]) % 2)
        rows.append(e)
    rows.append(list(y))
    B = matrix(ZZ, rows).hermite_form(include_zero_rows=False).change_ring(QQ) / 2
    assert B.nrows() == n
    Gn = (B * G * B.transpose()).change_ring(ZZ)
    Wn = (B * W * B.inverse()).change_ring(ZZ)
    assert Gn.det() == 1 and all(Gn[i, i] % 2 == 0 for i in range(n))
    return reduce_(Gn, Wn)


def step(G, W, rng):
    """One step of the uniform walk on pi-neighbours (v uniform in F_2^n, retried while v is in pi L).

    pi 近傍の一様なランダムウォークの 1 歩（v は F_2^n で一様、v が pi L に入る間は引き直す）。"""
    from sage.all import ZZ, vector
    n = G.nrows()
    while True:
        v = vector(ZZ, [rng.randrange(2) for _ in range(n)])
        r = neighbour(G, W, v)
        if r is not None:
            return r
