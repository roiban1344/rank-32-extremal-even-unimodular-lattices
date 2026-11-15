"""Functions of the 2-neighbour method that walks only representatives of the orbits of the automorphism group.

A random neighbour walk meets a class with probability proportional to 1/|Aut|, so it tends to miss the
classes with a large group (the mass stays slightly short). Splitting the classes of L/2L into orbits of
Aut L and building neighbours only from representatives makes the classes with a large group cheap, and
their neighbours are exhausted.
  isotropic_reps(G, gens, rng, cap)  the Aut-orbit representatives of the nonzero classes v of L/2L with
                                     (v, v) = 0 mod 4 (as integers: bit i is coordinate i)
  real_neighbour(G, v)               the 2-neighbour of an even Z-lattice of odd determinant (one per v)
  herm_neighbours(G, W, v)           the two neighbours of a Z[omega]-lattice at the isotropic line [v] for
                                     the inert prime 2
  real_seeds(p, n, s, k12_dir)       starting lattices of the genus of even p-elementary Z-lattices (sums
                                     of blocks)
  short_basis, counts, quick_isom, safe_canonical, forms, to_vec   helpers
  qfsafe                             rank32.forms with retries after a SystemError of cypari2
PARI is called through rank32.forms (time limits and retries with a new basis). The PARI stack is set to
1 GB (at most 2 GB) on import: a qfisom that does not finish keeps using it.
Used by computation C2 of the first paper (code/produce/element_types/orbit_representative_walk.py) and by
later computations (herm_neighbours).

自己同型群の軌道の代表だけを歩く 2-近傍法の関数。

乱択の近傍法は、類に当たる確率が 1/|Aut| に比例するので、群の大きい類を見落としやすい（質量が僅かに
足りないまま止まる）。L/2L の類を Aut L の軌道に分け、代表だけから近傍を作ると、群の大きい類ほど調べる
数が減り、その類の近傍を漏れなく尽くせる。
  isotropic_reps(G, gens, rng, cap)  L/2L の 0 でない類 v で (v, v) = 0 mod 4 のものの Aut 軌道の代表
                                     （整数で持つ。ビット i が座標 i）
  real_neighbour(G, v)               判別式が奇数の偶 Z-格子の 2-近傍（v ごとに 1 個）
  herm_neighbours(G, W, v)           Z[omega]-格子の、不活性素数 2 の等方直線 [v] での近傍 2 個
  real_seeds(p, n, s, k12_dir)       偶・p-基本の Z-格子の種の出発点（ブロックの直和）
  short_basis, counts, quick_isom, safe_canonical, forms, to_vec   補助
  qfsafe                             rank32.forms に cypari2 の SystemError の後の再試行を足したもの
PARI は rank32.forms 経由で呼ぶ（時間の上限と、基底を替えての再試行）。import 時に PARI のスタックを
1 GB（最大 2 GB）にする。終わらない qfisom はスタックを使い続ける。
第 1 の論文の計算 C2（code/produce/element_types/orbit_representative_walk.py）と、後の計算
（herm_neighbours）が使う。
"""
import itertools
import json
import resource

import numpy as np

from cysignals.alarm import AlarmInterrupt, alarm, cancel_alarm
from sage.all import ZZ, QQ, GF, matrix, vector, pari, CartanMatrix, block_diagonal_matrix

from cypari2.handle_error import PariError

from . import canonical as autgraph
from . import forms as _qfsafe
from . import rootless_unimodular_mass as ptroot

ptmass = ptroot.ptmass

pari.allocatemem(1 << 30, 2 << 30)  # a qfisom that does not finish keeps using the stack / 終わらない qfisom はスタックを使い続ける


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


LIMIT = 300
CANON_MAX = 30000  # the largest shell used for a canonical form (counted as by qfminim, +-v both) / 標準形を作るノルムの殻の本数の上限（qfminim の数え方、+- 両方）
A2 = matrix(ZZ, [[2, -1], [-1, 2]])
W2 = matrix(ZZ, [[1, 1], [-1, 0]])
C = lambda t, n: matrix(ZZ, CartanMatrix([t, n]))


def forms(G, W):
    """[G], or [G, G W^T] for a Z[omega]-lattice.

    [G]。Z[omega]-格子なら [G, G W^T]。"""
    return [G] if W is None else [G, G * W.transpose()]


def quick_isom(c, G, W, limit=8):
    """One qfisom with a short time limit: True or False, or None when it runs out of time.

    短い上限で 1 回だけ qfisom。True か False、時間切れなら None。"""
    fa, fb = forms(c["G"], c["W"]), forms(G, W)
    pa = pari(fa[0]) if W is None else pari(fa)
    pb = pari(fb[0]) if W is None else pari(fb)
    try:
        alarm(limit)
        r = pari.qfisom(pa, pb, [1, 0]) != 0  # fl = [1, 0] as in rank32.forms / rank32.forms と同じ fl = [1, 0]
        cancel_alarm()
        return r
    except (AlarmInterrupt, SystemError, PariError):
        return None
    finally:
        cancel_alarm()


def safe_canonical(G, W):
    """The canonical form (rank32.canonical), or False when the shell is too large or extra automorphisms appear.

    標準形（rank32.canonical）。殻が大きすぎるか、余分な自己同型が出るなら False。"""
    try:
        return autgraph.canonical(G, W, maxvec=CANON_MAX) or False
    except ValueError as ex:
        print("canonical form not usable:", str(ex)[:80], flush=True)
        return False


def counts(G):
    """The numbers of vectors of norm <= 2, <= 4, <= 6 (as counted by qfminim).

    ノルム <= 2, <= 4, <= 6 のベクトルの数（qfminim の数え方）。"""
    q = pari(G)
    return tuple(int(q.qfminim(m, 0)[0]) for m in (2, 4, 6))


def rss():
    """The current and the largest resident memory (GB), for the log.

    いまの RSS と最大値（GB）。ログ用。"""
    cur = 0.0
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS"):
                cur = int(line.split()[1]) / 2 ** 20
    return f"RSS {cur:.2f} GB, max {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20:.2f} GB"


def short_basis(G, W, rng, bound=6):
    """If a diagonal entry exceeds bound, change to a basis of vectors of norm <= bound; unchanged if they do not span.

    PARI enumerates everything up to the largest diagonal entry, which made the memory grow beyond 8 GB.

    対角が bound を超えるとき、基底をノルム <= bound のベクトルに取り替える。張らなければそのまま返す。
    PARI は対角の最大値まで全部を列挙するので、メモリが 8 GB を超えた。"""
    n = G.nrows()
    if max(G[i, i] for i in range(n)) <= bound:
        return G, W
    X = qfsafe.short_vectors(G, bound)
    if X.shape[0] < n:
        return G, W
    X = X[rng.sample(range(X.shape[0]), X.shape[0])]
    rows = []
    for r in X:
        c = matrix(ZZ, rows + [list(map(int, r))])
        if c.rank() == len(rows) + 1:
            rows.append(list(map(int, r)))
        if len(rows) == n:
            break
    if len(rows) < n:
        return G, W
    B = np.array(rows, dtype=np.int64)
    for _ in range(400):
        Bm = matrix(ZZ, B.tolist())
        if abs(Bm.det()) == 1:
            G2 = Bm * G * Bm.transpose()
            return G2, None if W is None else (Bm * W * Bm.inverse()).change_ring(ZZ)
        C = X @ np.array(Bm.inverse().numpy(dtype="float64"))
        hits = np.argwhere((np.abs(C) > 1e-6) & (np.abs(C) < 1 - 1e-6))
        if not len(hits):
            break
        i, j = hits[0]
        B[j] = X[i]
    return G, W


def to_vec(x, n):
    """The vector of Z^n with coordinates the bits of the integer x.

    整数 x のビットを座標とする Z^n のベクトル。"""
    return vector(ZZ, [(x >> i) & 1 for i in range(n)])


def isotropic_reps(G, gens, rng, cap):
    """The Aut-orbit representatives of the nonzero classes v of L/2L with (v, v) = 0 mod 4: returns (reps, total).

    The classes are integers (bit i = coordinate i). If there are more than cap representatives, a random
    sample of cap of them is returned.

    L/2L の 0 でない類 v で (v, v) = 0 mod 4 のものの Aut 軌道の代表：(代表, 総数) を返す。
    類は整数（ビット i = 座標 i）。代表が cap を超えれば、cap 個を乱択で返す。"""
    n = G.nrows()
    N = 1 << n
    q = np.zeros(1, dtype=np.uint8)
    for i in range(n):
        gm = sum(((int(G[i, j]) & 1) << j) for j in range(i))
        # parity(x & gm)
        par = np.zeros(1, dtype=np.uint8)
        for j in range(i):
            par = np.concatenate([par, par ^ ((gm >> j) & 1)])
        q = np.concatenate([q, q ^ par ^ ((int(G[i, i]) // 2) & 1)])
    # the images of the generators are not stored; they are computed from two small tables (storing them all took more than 8 GB for n = 26)
    # 生成元の像は持たず、2 つの小さい表から都度作る（n = 26 で像を全部持つと 8 GB を超えた）
    dtype = np.int32 if n < 31 else np.int64
    lo = n // 2
    tables = []
    for T in gens:
        masks = [sum(((int(T[j, i]) & 1) << j) for j in range(n)) for i in range(n)]  # e_i -> column i (the image of e_i under v -> v T^t) / e_i -> 第 i 列（v -> v T^t の e_i の像）
        t1 = np.zeros(1, dtype=dtype)
        for i in range(lo):
            t1 = np.concatenate([t1, t1 ^ masks[i]])
        t2 = np.zeros(1, dtype=dtype)
        for i in range(lo, n):
            t2 = np.concatenate([t2, t2 ^ masks[i]])
        tables.append((t1, t2))
    label = np.arange(N, dtype=dtype)
    CH = 1 << 22
    lomask = (1 << lo) - 1
    first = True
    while True:
        changed = False
        for t1, t2 in tables:
            for s in range(0, N, CH):
                x = np.arange(s, min(s + CH, N), dtype=dtype)
                y = t1[x & lomask] ^ t2[x >> lo]
                if first and not np.array_equal(q[y], q[x]):
                    raise ValueError("a generator does not preserve q")
                ly, lx = label[y], label[x]
                m = np.minimum(lx, ly)
                if not np.array_equal(m, lx):
                    label[x] = m
                    changed = True
                if not np.array_equal(m, ly):
                    label[y] = m
                    changed = True
        first = False
        if not changed:
            break
    idx = np.nonzero((label == np.arange(N)) & (q == 0))[0]
    idx = idx[idx != 0]
    del label, tables, q
    if cap and idx.shape[0] > cap:
        return rng.sample(list(map(int, idx)), cap), int(idx.shape[0])
    return [int(x) for x in idx], int(idx.shape[0])


def real_neighbour(G, v):
    """The 2-neighbour of the even Z-lattice G (odd determinant) at the class v: a list of at most one (Gram, None).

    偶 Z-格子 G（判別式は奇数）の類 v での 2-近傍：高々 1 個の (Gram, None) の列。"""
    n = G.nrows()
    nv = v * G * v
    if nv % 4:
        return []
    a = G * v
    k = next((i for i in range(n) if a[i] % 2), None)
    if k is None:
        return []
    if nv % 8:
        v = v + 2 * vector(ZZ, [1 if i == k else 0 for i in range(n)])
        a = G * v
    assert (v * G * v) % 8 == 0
    rows = []
    for i in range(n):
        e = [0] * n
        if i == k:
            e[k] = 2
        else:
            e[i] = 1
            e[k] = int(a[i]) % 2
        rows.append([2 * x for x in e])
    rows.append(list(v))
    B = matrix(ZZ, rows).hermite_form(include_zero_rows=False).change_ring(QQ) / 2
    Gn = (B * G * B.transpose()).change_ring(ZZ)
    assert Gn.det() == G.det() and all(Gn[i, i] % 2 == 0 for i in range(n))
    U = Gn.LLL_gram().transpose()
    return [(U * Gn * U.transpose(), None)]


def herm_neighbours(G, W, v):
    """The two neighbours at the isotropic line [v] for the inert prime 2: a list of (Gram, W).

    不活性素数 2 の等方直線 [v] の近傍 2 個：(Gram, W) の列。"""
    n = G.nrows()
    F2 = GF(2)
    if (v * G * v) % 4:
        return []
    wv = v * W
    a1, a2 = G * v, G * wv
    rows2 = matrix(F2, [list(a1), list(a2)])
    if rows2.rank() < 2:
        return []
    if (v * G * v) % 8:
        k = next(i for i in range(n) if a1[i] % 2)
        v = v + 2 * vector(ZZ, [1 if i == k else 0 for i in range(n)])
    assert (v * G * v) % 8 == 0
    K = rows2.right_kernel().basis_matrix()
    # the other lift: u with (v, u) even and (omega v, u) odd / もう 1 つの持ち上げ：(v, u) が偶数で (omega v, u) が奇数の u
    u = next(r for r in matrix(F2, 1, n, list(a1)).right_kernel().basis_matrix().rows() if sum(int(x) * int(y) for x, y in zip(r, a2)) % 2)
    v2 = v + 2 * vector(ZZ, [int(x) for x in u])
    assert (v2 * G * v2) % 8 == 0
    out = []
    for vv in (v, v2):
        ww = vv * W
        big = [[2 * int(x) for x in r] for r in K.rows()] + [[4 if i == j else 0 for j in range(n)] for i in range(n)] + [list(vv), list(ww)]
        B = matrix(ZZ, big).hermite_form(include_zero_rows=False).change_ring(QQ) / 2
        assert B.nrows() == n
        Gn = B * G * B.transpose()
        assert Gn.denominator() == 1
        Gn = Gn.change_ring(ZZ)
        assert Gn.det() == G.det()
        Wn = (B * W * B.inverse()).change_ring(ZZ)
        T = Gn.LLL_gram().transpose()
        Gl = T * Gn * T.transpose()
        Wl = (T * Wn * T.inverse()).change_ring(ZZ)
        assert Wl * Gl * Wl.transpose() == Gl
        out.append((Gl, Wl))
    return out


def real_seeds(p, n, s, k12_dir=None):
    """Starting lattices of the genus of even p-elementary Z-lattices of rank n and determinant p^s (sums of blocks).

    For p = 3, K_12 is one of the blocks when k12_dir contains unimod-z6.json.

    偶・p-基本・階数 n・行列式 p^s の Z-格子の種の出発点（ブロックの直和）。
    p = 3 では、k12_dir に unimod-z6.json があれば K_12 もブロックに入れる。"""
    E8 = C("E", 8)
    P = [("E8", E8), ("rtE8", p * E8)]
    if p == 3:
        E6 = C("E", 6)
        P += [("A2", C("A", 2)), ("E6", E6), ("rt3E6d", (3 * E6.inverse()).change_ring(ZZ))]
        k12 = None if k12_dir is None else k12_dir / "unimod-z6.json"
        if k12 is not None and k12.exists():
            P.append(("K12", matrix(ZZ, [c for c in json.loads(k12.read_text())["classes"] if c["inv"][0] == 0][0]["G"])))
    if p == 5:
        A4 = C("A", 4)
        P += [("A4", A4), ("rt5A4d", (5 * A4.inverse()).change_ring(ZZ))]
    if p == 7:
        A6 = C("A", 6)
        P += [("A6", A6), ("rt7A6d", (7 * A6.inverse()).change_ring(ZZ)), ("Q7", matrix(ZZ, [[2, 1], [1, 4]]))]
    dims = [g.nrows() for _, g in P]
    dets = [ZZ(g.det()).valuation(p) for _, g in P]
    res = []
    for combo in itertools.product(*[range(n // d + 1) for d in dims]):
        if sum(c * d for c, d in zip(combo, dims)) == n and sum(c * e for c, e in zip(combo, dets)) == s:
            res.append(block_diagonal_matrix([g for c, (_, g) in zip(combo, P) for _ in range(c)], subdivide=False))
    return res
