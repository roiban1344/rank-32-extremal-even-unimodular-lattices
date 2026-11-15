"""Direct and complete enumeration of the 2-elementary even lattices of rank d and minimum >= 4 (the candidates for the fixed and anti-fixed lattices of involutions of extremal lattices of rank 32).

Used by computations C5 and C7 of the first paper (Table 11). For an involution sigma of L, the fixed
lattice is F = L^sigma and the anti-fixed lattice is K = L^{-sigma}.
Completeness: such a lattice M satisfies sqrt(2) U <= M for a unimodular lattice U of rank d without
vectors of norm 1 [Neb13b, Lemma 4.9]; there are finitely many U [SPLAG, Table 16.7]. There is a
chain sqrt(2) U = Y_0 < Y_1 < ... < M with index 2 at each step (Lemma B.1 of the first paper);
every Y_j is a sublattice of M, hence even of minimum >= 4, and 2-elementary since
Y_j^# <= U / sqrt(2). So a breadth-first search from sqrt(2) U that adds one isotropic class of the
discriminant group (|x|^2 in 2Z) and keeps the minimum 4 reaches every M. At each step isometric
lattices are merged with qfisom (rank32.forms); the candidate classes are tried once per orbit of
Aut(M).
enumerate_classes(d, smin) returns the records written to min4-twoel-<d>[-s<smin>].json (the
Gram matrix, s, delta, |Aut|, the number of vectors of norm 4, the origin). For d >= 17 the lists of
unimodular lattices UNIMODULAR_FILES are read.

最小 4 以上の階数 d の 2-初等的な偶格子を、直接・完全に列挙する（階数 32 の極値格子の対合の固定格子・反固定格子の候補）。

第 1 の論文の計算 C5・C7（表 11）が使う。L の対合 sigma について、固定格子は F = L^sigma、反固定格子は
K = L^{-sigma} である。
完全性：そのような格子 M は、ノルム 1 を持たない階数 d のユニモジュラー格子 U について sqrt(2) U <= M を
満たす [Neb13b, Lemma 4.9]。U は有限個である [SPLAG, Table 16.7]。さらに sqrt(2) U = Y_0 < Y_1 < ... < M
（各ステージの指数 2）と鎖が取れ（第 1 の論文の補題 B.1）、途中の Y_j は M の部分格子なので偶・最小 4 以上、
かつ Y_j^# <= U / sqrt(2) から 2-初等的である。よって sqrt(2) U から始めて「判別式群の等方的な類
（|x|^2 in 2Z）を 1 本足して、最小 4 が保たれる拡大」を幅優先でたどれば全部に届く。各ステージで
qfisom（rank32.forms）により同型をまとめる。候補の類は Aut(M) の軌道ごとに 1 本だけ試す。
enumerate_classes(d, smin) は min4-twoel-<d>[-s<smin>].json に書く記録（Gram 行列、s、delta、|Aut|、
ノルム 4 の本数、出どころ）を返す。d >= 17 ではユニモジュラー格子の一覧 UNIMODULAR_FILES を読む。
"""
import json
import time
from pathlib import Path

from sage.all import ZZ, QQ, GF, matrix, vector, pari, CartanMatrix

from . import forms
from .io import asset

# lists of unimodular lattices without norm 1 for d >= 17, read in this order
# d >= 17 のノルム 1 を持たないユニモジュラー格子の一覧。この順に読む
UNIMODULAR_FILES = [asset("computation-outputs/C6/unimodular-min2.json"), asset("computation-outputs/C7/unimodular-min2.json")]
C = lambda t, n: matrix(ZZ, CartanMatrix([t, n]))


def dn_plus(n):
    """D_n^+ (also for n not divisible by 4): D_n with the spinor vector (1/2, ..., 1/2) added; built in standard coordinates, then LLL.

    D_n^+（n が 4 で割れないときも含む）：D_n にスピノル類のベクトル (1/2, ..., 1/2) を足したもの。座標は標準基底で作り LLL。"""
    rows = []
    for i in range(n - 1):
        r = [0] * n; r[i] = 1; r[i + 1] = -1; rows.append(r)
    r = [0] * n; r[n - 2] = 1; r[n - 1] = 1; rows.append(r)
    rows.append([QQ(1) / 2] * n)
    B = matrix(QQ, rows)
    H = (2 * B).change_ring(ZZ).hermite_form(include_zero_rows=False).change_ring(QQ) / 2
    G = (H * H.transpose()).change_ring(ZZ)
    U = G.LLL_gram(); return U.transpose() * G * U


def glue(G, xs):
    """The Gram matrix (after LLL) of the lattice of G with the rational vectors xs (basis coordinates) added.

    G の格子に有理座標のベクトル xs（基底座標）を足した格子の Gram（LLL 済み）。"""
    d = G.nrows()
    B = matrix(QQ, list(matrix.identity(QQ, d).rows()) + list(xs))
    den = B.denominator()
    H = (den * B).change_ring(ZZ).hermite_form(include_zero_rows=False).change_ring(QQ) / den
    Gn = H * G * H.transpose()
    assert all(t in ZZ for t in Gn.list())
    Gn = Gn.change_ring(ZZ)
    U = Gn.LLL_gram(); return U.transpose() * Gn * U


def unimodular_min2(d):
    """The unimodular lattices of rank d without vectors of norm 1 [SPLAG, Table 16.7]; d = 8, 12, 14, 15, 16 are built here, other ranks are read from UNIMODULAR_FILES.

    階数 d のノルム 1 を持たないユニモジュラー格子 [SPLAG, Table 16.7]。d = 8, 12, 14, 15, 16 はここで作り、他の階数は UNIMODULAR_FILES から読む。"""
    E8 = C("E", 8)
    if d == 8:
        return {"E8": E8}
    if d == 12:
        return {"D12+": dn_plus(12)}
    if d == 14:
        E7 = C("E", 7); R = matrix.block_diagonal([E7, E7], subdivide=False)
        Dinv = R.inverse(); x = Dinv.column(0) * 0
        # one vector of the non-trivial class of the dual of E_7 (norm 3/2) from each component: (E_7^2)^+
        # E_7 の双対の非自明な類（ノルム 3/2）を両方の成分から 1 本ずつ：(E_7^2)^+
        w = [j for j in range(7) if Dinv[j, j] not in ZZ][0]
        x = Dinv.column(w) + Dinv.column(7 + w)
        return {"E7^2+": glue(R, [x])}
    if d == 15:
        # A_15^+: A_15 with the glue [4] (norm 4*12/16 = 3) / A_15^+：A_15 に [4]（ノルム 4*12/16 = 3）の貼り合わせ
        A = C("A", 15); Dinv = A.inverse()
        return {"A15+": glue(A, [Dinv.column(3)])}
    if d == 16:
        D8 = C("D", 8); R = matrix.block_diagonal([D8, D8], subdivide=False); Dinv = R.inverse()
        # the glue (s, v), (v, s) of D_8^2 (s = spinor class, v = vector class) / D_8^2 の貼り合わせ (s, v), (v, s)（s = スピノル類、v = ベクトル類）
        sp, vc = 7, 0                                          # Bourbaki numbering: omega_8 is a spinor class (norm 2), omega_1 the vector class (norm 1) / Bourbaki の番号：omega_8 はスピノル類（ノルム 2）、omega_1 はベクトル類（ノルム 1）
        assert Dinv[sp, sp] == 2 and Dinv[vc, vc] == 1
        x1 = Dinv.column(sp) + Dinv.column(8 + vc); x2 = Dinv.column(vc) + Dinv.column(8 + sp)
        return {"E8^2": matrix.block_diagonal([E8, E8], subdivide=False), "D16+": dn_plus(16), "D8^2+": glue(R, [x1, x2])}
    for f in UNIMODULAR_FILES:
        if Path(f).exists():
            dat = json.load(open(f))
            if str(d) in dat:
                return {f"U{d}_{i}": matrix(ZZ, r["gram"]) for i, r in enumerate(dat[str(d)])}
    raise SystemExit(f"rank {d} is not available (no list of unimodular lattices of this rank in UNIMODULAR_FILES)")


def disc_data(G):
    """Generators x_i of the discriminant group G^#/G (2-elementary), as rational vectors in basis coordinates.

    判別式群 G^#/G（2-初等）の生成元 x_i（基底座標の有理ベクトル）。"""
    d = G.nrows(); D = G.inverse()
    gens = []; span = matrix(GF(2), 0, d)
    for j in range(d):
        c = D.column(j)
        v2 = vector(GF(2), [int(2 * x) % 2 for x in c])
        M2 = span.stack(v2)
        if M2.rank() > span.rank():
            gens.append(c); span = M2
    return gens


def extensions(G):
    """The even extensions of index 2 of G with minimum >= 4 (one candidate vector per orbit of Aut(G); isometric extensions may repeat).

    G の偶・指数 2 の拡大で最小 4 以上のもの（Aut(G) の軌道ごとに候補のベクトル 1 本。同型な拡大は重複しうる）。"""
    d = G.nrows(); D = G.inverse()
    gens = disc_data(G); k = len(gens)
    if k == 0:
        return []
    X = matrix(QQ, gens)                       # the rows are the generators / 行が生成元
    Q = X * G * X.transpose()                  # inner products of the generators (rational) / 生成元どうしの内積（有理）
    # isotropy: |sum c_i x_i|^2 in 2Z. q(c) = sum c_i Q_ii + 2 sum_{i<j} c_i c_j Q_ij mod 2 (Q_ii in Z/2, 2Q_ij in Z)
    # 等方性：|sum c_i x_i|^2 in 2Z。q(c) = sum c_i Q_ii + 2 sum_{i<j} c_i c_j Q_ij を mod 2 で（Q_ii in Z/2、2Q_ij in Z）
    q2 = [int(4 * Q[i, i]) % 8 for i in range(k)]           # 4|x_i|^2 mod 8
    b2 = [[int(4 * Q[i, j]) % 8 for j in range(k)] for i in range(k)]
    # classes of least norm 2 (bad classes): vectors of G^# of norm <= 2, sorted into classes
    # 最小ノルム 2 の類（悪い類）：G^# のノルム <= 2 のベクトルを類に直す
    K = 2 * D                                   # 2 (Gram of G^#, dual basis coordinates) / 2 (G^# の Gram、双対基底の座標)
    vs = pari(K.change_ring(ZZ)).qfminim(4)[2]                # 2|x|^2 <= 4 iff |x|^2 <= 2
    V = matrix(ZZ, vs).transpose() if len(vs) else matrix(ZZ, 0, d)
    bad = set()
    # from dual-basis coordinates a to the bits of the class: 2 D a mod 2 in the F_2 coordinates of the generators
    # 双対基底の座標 a から類のビット：2 D a mod 2 を生成元の F_2 座標に直す
    Spanmat = matrix(GF(2), [[int(2 * x) % 2 for x in g] for g in gens])
    for a in V.rows():
        y = vector(GF(2), [int(t) % 2 for t in (2 * D * vector(ZZ, a))])
        c = Spanmat.solve_left(y)
        bad.add(sum(int(c[i]) << i for i in range(k)))
    # action of generators of Aut(M) on the discriminant group F_2^k (as matrices). x and g x give isometric extensions,
    # so only one class per orbit is tried
    # Aut(M) の生成元の、判別式群 F_2^k への作用（行列）。x と g x は同型な拡大を与えるので軌道ごとに 1 本だけ試す
    acts = []
    try:
        A = forms.qfauto(G)
        for g in A[1]:
            gm = matrix(ZZ, g)
            cols = []
            for i in range(k):
                y = vector(GF(2), [int(t) % 2 for t in (2 * (gm * gens[i]))])
                c = Spanmat.solve_left(y)
                cols.append(sum(int(c[t]) << t for t in range(k)))
            acts.append(cols)
    except Exception as e:
        print(f"  (qfauto not available: {e}; no reduction by orbits)", flush=True)
    def act(cols, mask):
        r = 0; i = 0
        while mask:
            if mask & 1:
                r ^= cols[i]
            mask >>= 1; i += 1
        return r
    seen = bytearray(1 << k)
    out = []
    for mask in range(1, 1 << k):
        if seen[mask] or mask in bad:
            continue
        idx = [i for i in range(k) if mask >> i & 1]
        v = sum(q2[i] for i in idx) + 2 * sum(b2[i][j] for a_, i in enumerate(idx) for j in idx[a_ + 1:])
        if v % 8:                              # 4|x|^2 = 0 mod 8 iff |x|^2 in 2Z
            continue
        # colour the orbit / 軌道を塗る
        stack = [mask]; seen[mask] = 1
        while stack:
            u = stack.pop()
            for cols in acts:
                w = act(cols, u)
                if not seen[w]:
                    seen[w] = 1; stack.append(w)
        x = sum((gens[i] for i in idx), vector(QQ, [0] * d))
        out.append(x)
    return out


def invariant(G):
    """(det, number of vectors of norm <= 4, of norm <= 6) as counted by qfminim; used to sort lattices before qfisom.

    qfminim の数え方の (det、ノルム 4 以下の本数、ノルム 6 以下の本数)。qfisom の前に格子を分けるのに使う。"""
    q = pari(G)
    return (int(G.det()), int(q.qfminim(4, 0)[0]), int(q.qfminim(6, 0)[0]))


def output_name(d, smin=0):
    """The file name of the output: min4-twoel-<d>.json, or min4-twoel-<d>-s<smin>.json.

    出力のファイル名：min4-twoel-<d>.json か min4-twoel-<d>-s<smin>.json。"""
    return f"min4-twoel-{d}.json" if smin == 0 else f"min4-twoel-{d}-s{smin}.json"


def enumerate_classes(d, smin=0):
    """All classes of rank d down to s = smin (the extension stops below smin), as the dict written to output_name(d, smin).

    階数 d の全部の類（s が smin を下回る拡大はしない）。output_name(d, smin) に書く dict で返す。"""
    t0 = time.time()
    level = []
    for name, U in unimodular_min2(d).items():
        assert U.det() == 1 and int(pari(U).qfminim(1, 0)[0]) == 0, name
        level.append((name, 2 * U))
    classes = []          # classes of all steps / すべてのステージの類
    byinv = {}
    def add(G, origin):
        inv = invariant(G)
        for j in byinv.get(inv, []):
            if forms.qfisom(classes[j]["G"], G) is not None:
                return None
        rec = dict(G=G, origin=origin, s=int(ZZ(G.det()).valuation(2)))
        byinv.setdefault(inv, []).append(len(classes)); classes.append(rec)
        return rec
    cur = []
    for name, G in level:
        r = add(G, f"\u221a2\u00b7{name}")   # the origin "sqrt(2)*name", written as in the lists of components / 出どころ "sqrt(2)*name"。成分の一覧と同じ書き方
        if r: cur.append(r)
    step = 0
    while cur:
        step += 1
        nxt = []
        for r in cur:
            if r["s"] - 2 < smin:
                continue
            for x in extensions(r["G"]):
                Gn = glue(r["G"], [x])
                if int(pari(Gn).qfminim(2, 0)[0]):
                    continue
                rn = add(Gn, r["origin"] + f" +{step}")
                if rn: nxt.append(rn)
        print(f"step {step}: new classes {len(nxt)} (total {len(classes)}, {time.time()-t0:.0f} s)", flush=True)
        cur = nxt
    recs = []
    for r in classes:
        G = r["G"]; D = G.inverse()
        delta = 0 if all(D[i, i] in ZZ for i in range(d)) else 1
        recs.append(dict(gram=[[int(x) for x in row] for row in G.rows()], s=r["s"], delta=delta, origin=r["origin"],
                         aut=int(forms.qfauto(G)[0]), norm4=int(pari(G).qfminim(4, 0)[0])))
        print(f"  s={r['s']} delta={delta} |Aut|={recs[-1]['aut']} norm4={recs[-1]['norm4']}  <- {r['origin']}", flush=True)
    print(f"-> rank {d}: {len(recs)} 2-elementary even lattices of minimum >= 4", flush=True)
    return dict(d=d, smin=smin, n=len(recs), classes=recs)
