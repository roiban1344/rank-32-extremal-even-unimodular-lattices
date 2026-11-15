"""The type (8, 24) of involutions: the lattices M of rank 24 (det 2^8, minimum 4) from subspaces W of a Niemeier lattice, their isometry classes and automorphism groups.

Used by computation C6 of the first paper (Table 11).
Notation of Appendix B.3 of the paper: M, M^# and the Niemeier lattice N of M_from_W are K, K^# and P
there; in analyse, subgroup_orbit and census, N = sqrt(2) M^# and N' are J and J'.
Paths: UNIMOD and ORBITS (the orbits <P>/W.jsonl) are read by default from computation-outputs/C6 of
the release assets; OUTDIR (default results/C6 under the current directory) receives the classes.
M_from_W(G, W) is the lattice {x in N : (x, w) = 0 mod 2 for w in W} + 2N, LLL-reduced, for the
Niemeier lattice N with Gram matrix G.
analyse(lib, M) decides isometry through the canonical form of the graph of the 720 vectors of
norm 2 of M^# (norm 4 in sqrt(2) M^#), with the colours of the vertices from the diagonal. When the
720 vectors generate sqrt(2) M^# over Z (the HNF of lattice_index has determinant 1),
Aut(M) = Aut(graph) and the canonical form is a complete invariant. Otherwise the orbit of
N/N' in N'^#/N' (N' spanned by the 720 vectors) under Aut(graph) gives Aut(M) and a complete
invariant (subgroup_orbit). For lattices where that orbit is too large (the most symmetric ones),
|Aut M| comes from qfauto and isometry is decided later by qfisom on the lattices with equal
invariants. census records, for the 270 maximal singular subspaces T' of M^#/M, the distribution
of the numbers of roots of the Niemeier lattices over them.
extremal_lattice(M) glues M with sqrt(2) E_8 to an extremal lattice of rank 32 (any
anti-isometry between hyperbolic bases: W(E_8) -> O^+(8, 2) is onto).
The C part _colored_gram_canon.c is compiled on first use (gcc -O3 -shared -fPIC, linked to bliss)
into $RANK32_BUILD (default ~/.cache/rank32); bliss is looked up as in rank32.canonical.

対合の (8, 24) 型：Niemeier 格子の部分空間 W から得る階数 24 の格子 M（判別式 2^8、最小 4）と、その同型類と自己同型群。

第 1 の論文の計算 C6（表 11）が使う。
論文の付録 B.3 の記号との対応：M、M^#、M_from_W の Niemeier 格子 N は、そこでの K、K^#、P である。
analyse、subgroup_orbit、census の中の N = sqrt(2) M^# と N' は、J と J' である。
パス：UNIMOD と ORBITS（軌道 <P>/W.jsonl）は、既定では Release のアセットの computation-outputs/C6 から読む。
OUTDIR（既定は今のディレクトリの下の results/C6）に類を書く。
M_from_W(G, W) は、Gram 行列 G の Niemeier 格子 N について、{x in N : w in W に対し (x, w) = 0 mod 2} + 2N
を LLL したものである。
analyse(lib, M) は、M^# のノルム 2 のベクトル 720 本（sqrt(2) M^# ではノルム 4）のグラフ（頂点の色は対角）の
標準形で同型を判定する。720 本が sqrt(2) M^# を Z 上で生成する（lattice_index の HNF の行列式が 1）ときは
Aut(M) = Aut(グラフ) で、標準形が完全不変量になる。生成しないときは、N/N' in N'^#/N'（N' は 720 本の張る格子）の
Aut(グラフ) の軌道から Aut(M) と完全不変量を求める（subgroup_orbit）。その軌道が大きすぎる格子（対称性の
最も高いもの）は qfauto で位数を求め、同型は後で不変量の一致する組を qfisom で照合する。census は、
M^#/M の極大特異部分空間 T'（270 個）の上の Niemeier 格子のルートの本数の分布を記録する。
extremal_lattice(M) は M と sqrt(2) E_8 を貼り合わせて階数 32 の極値格子にする（双曲基底どうしの反等長で
よい。W(E_8) -> O^+(8, 2) が全射だから）。
C の部分 _colored_gram_canon.c は、初めて使うときに $RANK32_BUILD（既定は ~/.cache/rank32）へコンパイル
する（gcc -O3 -shared -fPIC、bliss をリンクする）。bliss の探し方は rank32.canonical と同じ。
"""
import ctypes
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from sage.all import GF, QQ, ZZ, matrix, pari
from sage.libs.gap.libgap import libgap

from . import forms as qfsafe
from .io import asset

UNIMOD = asset("computation-outputs/C6/unimodular-min2.json")
ORBITS = asset("computation-outputs/C6")  # <P>/W.jsonl, the orbits W (input) / 軌道 W の <P>/W.jsonl（入力）
OUTDIR = Path("results/C6")             # <P>/classes-*.jsonl (output) / <P>/classes-*.jsonl（出力）
NAMES = {"Leech": 179, "A1": 177, "A2": 178}
D = 24
_SRC = Path(__file__).with_name("_colored_gram_canon.c")


def _build_dir():
    d = os.environ.get("RANK32_BUILD")
    if d:
        return Path(d)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "rank32"


def _bliss_prefix():
    return Path(os.environ.get("BLISS_PREFIX") or os.environ.get("CONDA_PREFIX") or sys.prefix)


def load_canon():
    """Compile _colored_gram_canon.c if needed and load it with ctypes.

    _colored_gram_canon.c を必要ならコンパイルし、ctypes で読み込む。"""
    build = _build_dir()
    exe = build / ("colored-gram-canon-" + hashlib.sha256(_SRC.read_bytes()).hexdigest()[:12] + ".so")
    if not exe.exists():
        env = _bliss_prefix()
        build.mkdir(parents=True, exist_ok=True)
        tmp = str(exe) + f".{os.getpid()}.tmp"
        subprocess.run(["gcc", "-O3", "-shared", "-fPIC", str(_SRC), "-o", tmp, "-I" + str(env / "include"), "-L" + str(env / "lib"),
                        "-Wl,-rpath," + str(env / "lib"), "-lbliss", "-lstdc++"], check=True)
        os.replace(tmp, exe)
    lib = ctypes.CDLL(str(exe))
    lib.canon.restype = ctypes.c_int
    return lib


def certificate(lib, Gv):
    """The hash of the canonical form of a Gram matrix (the diagonal colours the vertices), the canonical labelling, and generators of the automorphisms (arrays of permutations).

    Gram（対角が頂点の色）の標準形のハッシュ、標準番号、自己同型の生成元（置換の配列）。"""
    n = Gv.shape[0]
    q = Gv.astype(np.int8)
    col = np.diag(Gv).astype(np.uint32)
    lab = np.zeros(n, dtype=np.uint32)
    maxg = 1000
    gens = np.zeros((maxg, n), dtype=np.uint32)
    ptr = lambda a: a.ctypes.data_as(ctypes.c_void_p)
    k = lib.canon(ctypes.c_uint32(n), ptr(q), ptr(col), ptr(lab), ptr(gens), ctypes.c_uint32(maxg))
    assert k >= 0
    assert sorted(lab.tolist()) == list(range(n))
    inv = np.argsort(lab)
    canon = q[inv][:, inv]
    return hashlib.sha256(canon.tobytes()).hexdigest()[:24], lab, gens[:k]


def group_order(gens, n):
    """The order of the permutation group generated by gens (GAP).

    gens の生成する置換群の位数（GAP）。"""
    if len(gens) == 0:
        return 1
    return int(libgap.Group([libgap.PermList([int(x) + 1 for x in g]) for g in gens]).Size())


def M_from_W(G, W):
    """The lattice M (LLL-reduced Gram matrix) of the subspace W (masks) of N/2N, N with Gram matrix G.

    Gram 行列 G の N について、N/2N の部分空間 W（mask）の格子 M（LLL 済みの Gram 行列）。"""
    A = matrix(GF(2), [[(w >> i) & 1 for i in range(D)] for w in W]) * G.change_ring(GF(2))
    K = A.right_kernel().basis_matrix().change_ring(ZZ)
    H = K.stack(2 * matrix.identity(ZZ, D)).hermite_form(include_zero_rows=False)
    assert H.nrows() == D and abs(H.det()) == 16
    Mg = H * G * H.transpose()
    U = Mg.LLL_gram()
    return U.transpose() * Mg * U


def standard_subspaces():
    """The 270 maximal singular subspaces of F_2^8 with q(x) = x_0 x_4 + x_1 x_5 + x_2 x_6 + x_3 x_7 (arrays of their 15 points).

    F_2^8、q(x) = x_0 x_4 + x_1 x_5 + x_2 x_6 + x_3 x_7 の極大特異部分空間 270 個（各 15 点の整数の表）。"""
    def q(x):
        return bin((x & 15) & (x >> 4)).count("1") & 1
    def b(x, y):
        return (q(x ^ y) + q(x) + q(y)) & 1
    pts = [x for x in range(1, 256) if q(x) == 0]
    assert len(pts) == 135
    level = {frozenset([p]) for p in pts}
    for _ in range(3):
        nxt = set()
        for sp in level:
            for d in pts:
                if d not in sp and all(b(d, e) == 0 for e in sp):
                    nxt.add(frozenset(sp | {d} | {d ^ e for e in sp}))
        level = nxt
    assert len(level) == 270
    return np.array(sorted(sorted(s) for s in level), dtype=np.int32)


SUBS = standard_subspaces()
PACK = 1 << np.arange(D, dtype=np.int64)


def hyperbolic_basis(Dn):
    """A hyperbolic basis e_1..e_4, f_1..f_4 of the discriminant group (rank 8, type +, delta = 0), as integer coordinates a (v = (Dn/2) a), with q(v) = (a^T Dn a)/2 mod 2.

    判別式群（階数 8、+ 型、delta = 0）の双曲基底 e_1..e_4, f_1..f_4 を、整数座標 a（v = (Dn/2) a）で返す。q(v) = (a^T Dn a)/2 mod 2。"""
    d = Dn.shape[0]
    pack = 1 << np.arange(d, dtype=np.int64)
    Dmod = Dn % 2
    # coordinates of 2 M^# / 2 M: a -> Dn a mod 2 / 2 M^# / 2 M の座標：a -> Dn a mod 2
    cols = (Dmod.T @ pack).astype(np.int64)
    pick = []
    ech = {}
    for j in range(d):
        x = int(cols[j])
        for p, v in ech.items():
            if x >> p & 1:
                x ^= v
        if x:
            p = x.bit_length() - 1
            for t in list(ech):
                if ech[t] >> p & 1:
                    ech[t] ^= x
            ech[p] = x
            pick.append(j)
    assert len(pick) == 8
    E = np.eye(d, dtype=np.int64)[pick]
    # to the standard form: take e_i, f_i in turn (q is (a^T Dn a)/2 mod 2, b is a^T Dn a' mod 2)
    # 標準形へ：e_i, f_i を順に取る（q は (a^T Dn a)/2 mod 2、b は a^T Dn a' mod 2）
    def qv(a):
        return int(a @ Dn @ a) // 2 & 1
    def bv(a, c):
        return int(a @ Dn @ c) & 1
    pool = [E[i] for i in range(8)]
    es, fs = [], []
    for _ in range(4):
        # pool is orthogonal to the pairs already made and of type +. Take a singular a and c with b(a, c) = 1,
        # make c singular, and project the rest
        # pool は既に作った対に直交で + 型。特異な元 a と b(a, c) = 1 の c を取り、c を特異に直し、残りを射影する
        a = None
        for mask in range(1, 1 << len(pool)):
            x = sum((pool[i] for i in range(len(pool)) if mask >> i & 1), np.zeros(d, dtype=np.int64))
            if qv(x) == 0:
                a = x
                break
        c = next(x for x in pool if bv(a, x))
        if qv(c):
            c = c + a
        assert qv(a) == 0 and qv(c) == 0 and bv(a, c) == 1
        rest = [x + bv(x, c) * a + bv(x, a) * c for x in pool]
        # drop two dimensions / 2 次元ぶん落とす
        keep = []
        ech = {}
        for x in rest:
            k = int(((Dmod @ (x % 2)) % 2) @ pack)
            for pv, v in ech.items():
                if k >> pv & 1:
                    k ^= v
            if k:
                pv = k.bit_length() - 1
                for t in list(ech):
                    if ech[t] >> pv & 1:
                        ech[t] ^= k
                ech[pv] = k
                keep.append(x)
        pool = keep
        assert all(bv(x, a) == 0 and bv(x, c) == 0 for x in pool)
        assert len(pool) == 6 - 2 * len(es)
        es.append(a)
        fs.append(c)
    return es + fs


def census(Dn, V):
    """The numbers n(p) of vectors of V (vectors of norm 4 of N in coordinates) per class, and the distribution of the numbers of roots over the T'.

    V（N の座標のノルム 4 のベクトル）の類ごとの本数 n(p) と、T' ごとのルートの本数の分布。"""
    Dmod = Dn % 2
    basis = hyperbolic_basis(Dn)
    table = {}
    for i in range(256):
        a = sum((basis[j] for j in range(8) if i >> j & 1), np.zeros(D, dtype=np.int64))
        table[int(((Dmod @ (a % 2)) % 2) @ PACK)] = i
    assert len(table) == 256
    keys = ((V @ Dmod) % 2) @ PACK
    n = np.zeros(256, dtype=np.int64)
    for k, c in Counter(keys.tolist()).items():
        n[table[k]] = c
    assert n[0] == 0
    per = n[SUBS].sum(1)
    assert per.sum() == 80 * 270
    return dict(zip(*[x.tolist() for x in np.unique(per, return_counts=True)])), dict(zip(*[x.tolist() for x in np.unique(n[[x for x in range(1, 256) if bin((x & 15) & (x >> 4)).count("1") % 2 == 0]], return_counts=True)]))


def extremal_lattice(Mg):
    """The Gram matrix of the lattice of rank 32 glued from M and sqrt(2) E_8 (the anti-isometry of discriminants maps hyperbolic basis to hyperbolic basis; W(E_8) -> O^+(8, 2) is onto, so any choice gives the same lattice).

    M と sqrt(2) E_8 を貼り合わせた階数 32 の格子の Gram（判別式の反等長は双曲基底どうしで取る。W(E_8) -> O^+(8, 2) の全射でどれでも同型）。"""
    from sage.all import CartanMatrix, block_diagonal_matrix
    C = matrix(ZZ, CartanMatrix(["E", 8]))
    Dm = (2 * Mg.inverse()).change_ring(ZZ)
    De = C.inverse().change_ring(ZZ)
    bm = hyperbolic_basis(np.array(Dm.list(), dtype=np.int64).reshape(D, D))
    be = hyperbolic_basis(np.array(De.list(), dtype=np.int64).reshape(8, 8))
    # coordinates: the basis of sqrt(2) E_8 (8 vectors) and of M (24 vectors); v = (Dn/2) a
    # 座標：sqrt(2) E_8 の基底（8 本）と M の基底（24 本）。v = (Dn/2) a
    X = [list((De * matrix(ZZ, [x.tolist()]).transpose()).list()) for x in be]
    Y = [list((Dm * matrix(ZZ, [x.tolist()]).transpose()).list()) for x in bm]
    rows = [[QQ(t) / 2 for t in x] + [QQ(t) / 2 for t in y] for x, y in zip(X, Y)]
    Bq = matrix(QQ, rows)
    Id = matrix.identity(QQ, 32)
    A = Id.stack(Bq)
    H = (2 * A).change_ring(ZZ).hermite_form(include_zero_rows=False).change_ring(QQ) / 2
    assert H.nrows() == 32
    T = block_diagonal_matrix([2 * C, Mg], subdivide=False)
    L = H * T * H.transpose()
    assert all(t in ZZ for t in L.list())
    L = L.change_ring(ZZ)
    U = L.LLL_gram()
    return U.transpose() * L * U


def lattice_index(V):
    """The index in Z^24 of the lattice spanned by the rows V (0 if the rank is too small).

    行ベクトル V の張る格子の、Z^24 での指数（階数が足りなければ 0）。"""
    H = matrix(ZZ, V.tolist()).hermite_form(include_zero_rows=False)
    return int(abs(H.det())) if H.nrows() == D else 0


def subgroup_orbit(Dn, V, lab, gens, index, cap=20000):
    """The orbit of N/N' in N'^#/N' under A = Aut(graph), where N' (spanned by the 720 vectors) has index `index` in N.

    Aut(M) is the stabilizer of N/N', and [A : Aut(M)] is the size of the orbit. The least element of
    the orbit is written in canonical coordinates, so together with the canonical form of the graph
    it is a complete invariant of M. Coordinates: x -> ((x, v_j))_{j in J}, J the first 24
    independent vectors in canonical order. Returns (size of the orbit, hash of the least element),
    or None beyond cap.

    720 本が張る N' は N の指数 index の部分格子。A = Aut(グラフ) の、N/N' in N'^#/N' の軌道を求める。

    Aut(M) は N/N' の固定部分群で、[A : Aut(M)] = 軌道の大きさ。軌道の最小元は標準の座標で書くので、グラフの
    標準形と合わせて M の完全不変量になる。座標：x -> ((x, v_j))_{j in J}（J は標準の順で最初の 24 本の独立な
    ベクトル）。戻り値は (軌道の大きさ, 最小元のハッシュ)。cap を超えたら None。"""
    inv = np.argsort(lab)
    Vc = V[inv]
    Gc = Vc @ Dn @ Vc.T
    J = []
    ech = matrix(QQ, 0, D)
    for j in range(len(Vc)):
        t = ech.stack(matrix(QQ, [Vc[j].tolist()]))
        if t.rank() > ech.nrows():
            ech = t
            J.append(j)
            if len(J) == D:
                break
    J = np.array(J)
    R = np.array(matrix(ZZ, Gc[:, J].tolist()).hermite_form(include_zero_rows=False).list(), dtype=np.int64).reshape(D, D)
    piv = [int(np.flatnonzero(R[i])[0]) for i in range(D)]

    def reduce(y):
        y = y.copy()
        for i in range(D):
            y -= (y[piv[i]] // R[i, piv[i]]) * R[i]
        return tuple(y.tolist())

    # all elements of N/N' / N/N' の元を全部作る
    unit = [reduce((np.eye(D, dtype=np.int64)[i] @ Dn @ Vc[J].T)) for i in range(D)]
    elems = {tuple([0] * D)}
    frontier = list(elems)
    while frontier:
        nxt = []
        for e in frontier:
            for u in unit:
                z = reduce(np.array(e, dtype=np.int64) + np.array(u, dtype=np.int64))
                if z not in elems:
                    elems.add(z)
                    nxt.append(z)
        frontier = nxt
    assert len(elems) == index
    # action of the generators on the coordinates: (x g, v_j) = (x, v_{pi^{-1}(j)}); write v_{pi^{-1}(j)} in terms of v_J
    # 生成元の座標への作用：(x g, v_j) = (x, v_{pi^{-1}(j)})。v_{pi^{-1}(j)} を v_J で書く
    Ginv = matrix(QQ, Gc[np.ix_(J, J)].tolist()).inverse()
    acts = []
    for g in gens:
        pc = np.zeros(len(Vc), dtype=np.int64)
        pc[lab] = lab[g]                      # the permutation pi in canonical numbering / 標準の番号での置換 pi
        pinv = np.argsort(pc)
        C = matrix(QQ, Gc[np.ix_(pinv[J], J)].tolist()) * Ginv
        den = int(C.denominator())
        acts.append((np.array((den * C).change_ring(ZZ).list(), dtype=np.int64).reshape(D, D), den))

    def apply(act, S):
        C, den = act
        out = []
        for e in S:
            y = C @ np.array(e, dtype=np.int64)
            assert not (y % den).any()
            out.append(reduce(y // den))
        return frozenset(out)

    start = frozenset(elems)
    orbit = {start}
    frontier = [start]
    while frontier:
        nxt = []
        for S in frontier:
            for act in acts:
                T = apply(act, S)
                if T not in orbit:
                    orbit.add(T)
                    nxt.append(T)
                    if len(orbit) > cap:
                        return None
        frontier = nxt
    best = min(tuple(sorted(S)) for S in orbit)
    return len(orbit), hashlib.sha256(json.dumps(best).encode()).hexdigest()[:24]


def analyse(lib, Mg):
    """The invariant (cert), |Aut M| and the census of a lattice M of type (8, 24); returns (record, Mg).

    (8, 24) 型の格子 M の不変量（cert）、|Aut M|、census。(記録, Mg) を返す。"""
    p = pari(Mg)
    assert Mg.det() == 256 and all(Mg[i, i] % 2 == 0 for i in range(D))
    assert int(p.qfminim(2, 0)[0]) == 0
    assert int(p.qfminim(4, 0)[0]) == 12240
    Dn_s = (2 * Mg.inverse()).change_ring(ZZ)
    assert all(Dn_s[i, i] % 2 == 0 for i in range(D))
    V0 = qfsafe.short_vectors(Dn_s, 4)
    Dn = np.array(Dn_s.list(), dtype=np.int64).reshape(D, D)
    assert len(V0) == 360 and ((V0 @ Dn * V0).sum(1) == 4).all()
    V = np.vstack([V0, -V0])
    index = lattice_index(V0)
    assert index >= 1, "the 720 vectors have rank < 24"
    t_census, n_dist = census(Dn, V)
    cert, lab, gens = certificate(lib, V @ Dn @ V.T)
    aut = group_order(gens, len(V))
    method, orbit = "graph", 1
    if index > 1:
        res = subgroup_orbit(Dn, V, lab, gens, index)
        if res is None:
            # a large orbit means the most symmetric lattices, where qfauto is quick; isometry is checked by qfisom at the merge
            # 軌道が大きいのは対称性の最も高い格子で、qfauto がすぐ終わる。同型は merge のステージで qfisom で確かめる
            aut = int(qfsafe.qfauto(Mg)[0])
            method, orbit = "qfauto", 0
            cert = "hi:" + hashlib.sha256(json.dumps([index, aut, sorted(t_census.items()), sorted(n_dist.items())]).encode()).hexdigest()[:20]
        else:
            orbit, h = res
            assert aut % orbit == 0
            aut //= orbit
            method = "orbit"
            cert = cert + ":" + h
    return dict(cert=cert, index=index, method=method, orbit=orbit, aut=aut, census=t_census, np=n_dist), Mg
