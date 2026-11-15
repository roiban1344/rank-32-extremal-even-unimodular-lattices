"""Building the lattice L of rank 32 from a gluing (F, K, phi), and reading phi back from (L, sigma); the functions used by the checks of the involution gluing.

Used by computation C4 of the first paper (Table 11; the checks of the gluing) and by C34 (the
numbers e(sigma), which build lattices from gluings as well).
Paths: LISTS (default computation-outputs/C4 of the release assets) holds the lists <pair>.json and
the checkpoints <pair>.ckpt.json of the gluing; OUT (default results/C4 under the current directory)
receives the outputs of the checks and the temporary files of the bliss tool.
  build_gram(GF, GK, phi): L = F + K + {(x_i, phi(x_i))} and sigma (+1 on F, -1 on K) in an LLL
    basis of L.
  glue_of(G, sigma, GF0, GK0): phi in the bases of the components (GF0, GK0), or None when the fixed
    lattice L^sigma or the anti-fixed lattice L^{-sigma} is not isometric to the component.
    short_profile (numbers of short vectors and a table of inner products) is compared first, since
    qfisom may not return on large groups that are not isomorphic.
  centralizer_order(gens, sigma, G): |Aut L| and |C(sigma)| through a permutation group on
    shortest vectors (GAP).
  aut_tool(G, tag): |Aut L| and generators with rank32.bliss_automorphisms.
  class_table(pair, g, res, kFK): key -> record for all classes of a closed pair.

対合の貼り合わせ (F, K, phi) から階数 32 の格子 L を作り、(L, sigma) から phi を読み戻す。対合の貼り合わせの検算が使う関数。

第 1 の論文の計算 C4（表 11。貼り合わせの検算）と C34（数 e(sigma)。これも貼り合わせから格子を作る）が使う。
パス：LISTS（既定は Release のアセットの computation-outputs/C4）に、貼り合わせの一覧 <組>.json と
チェックポイント <組>.ckpt.json がある。OUT（既定は今のディレクトリの下の results/C4）に、検算の出力と
bliss の道具の一時ファイルを書く。
  build_gram(GF, GK, phi)：L = F + K + {(x_i, phi(x_i))} と sigma（F で +1、K で -1）を、L の LLL 基底で。
  glue_of(G, sigma, GF0, GK0)：成分 (GF0, GK0) の基底での phi。sigma の固定格子 L^sigma か反固定格子
    L^{-sigma} が成分と等長でなければ None。先に short_profile（短いベクトルの本数と内積の表）を比べる。
    群が大きくて同型でない組では、qfisom が返らないことがあるからである。
  centralizer_order(gens, sigma, G)：最短ベクトルの上の置換群（GAP）で |Aut L| と |C(sigma)|。
  aut_tool(G, tag)：rank32.bliss_automorphisms で |Aut L| と生成元。
  class_table(pair, g, res, kFK)：閉じた組の全部の類について、キー -> 記録。
"""
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from sage.all import ZZ, QQ, GF, matrix, vector, pari, libgap

from . import gluing as glue
from . import forms as qfsafe
from . import bliss_automorphisms
from . import involution_gluing as gg

LISTS = gg.DATA / "C4"     # the lists <pair>.json and checkpoints of the gluing (input) / 貼り合わせの一覧とチェックポイント（入力）
OUT = gg.OUT               # outputs of the checks, by default results/C4 / 検算の出力。既定は results/C4


def disc_gens(G):
    """Basis of the discriminant group chosen as in antiisometry_sampler.disc (rational vectors), and the matrix that reads coordinates.

    antiisometry_sampler.disc と同じ基底の取り方の、判別式群の基底の有理ベクトルと、座標を読む行列。"""
    d = G.nrows(); Xi = G.inverse(); gens = []; span = matrix(GF(2), 0, d)
    for j in range(d):
        c = Xi.column(j); v2 = vector(GF(2), [int(2 * x) % 2 for x in c])
        if span.stack(v2).rank() > span.rank():
            gens.append(c); span = span.stack(v2)
    J = span.pivots(); Sinv = span.matrix_from_columns(J).inverse()
    return gens, J, Sinv


def mask_of(v, J, Sinv):
    """The mask of the class of v in F^# (rational coordinates).

    F^# の元 v（有理座標）の類の mask。"""
    w2 = vector(GF(2), [int(2 * v[j]) % 2 for j in J]) * Sinv
    return sum(int(w2[l]) << l for l in range(len(J)))


def build_gram(GF_, GK_, phi):
    """The Gram matrix of L = F + K + {(x_i, phi(x_i))} and the matrix of sigma (+1 on F, -1 on K) in a basis of L.

    L = F + K + {(x_i, phi(x_i))} の Gram と、sigma（F で +1、K で -1）の行列を L の基底で返す。"""
    a, b = GF_.nrows(), GK_.nrows()
    gF, _, _ = disc_gens(GF_); gK, _, _ = disc_gens(GK_)
    rows = []
    for i in range(a + b):
        v = [QQ(0)] * (a + b); v[i] = QQ(1); rows.append(v)
    for i, xi in enumerate(gF):
        m = int(phi[i]); yi = vector(QQ, [0] * b)
        for l, g in enumerate(gK):
            if (m >> l) & 1:
                yi += g
        rows.append(list(xi) + list(yi))
    M = matrix(QQ, rows)
    den = 2
    H = (den * M).change_ring(ZZ).hermite_form(include_zero_rows=False)
    B = H.change_ring(QQ) / den                                             # rows are a basis of L (coordinates of F + K) / 行が L の基底（F + K の座標）
    A = matrix.block_diagonal([GF_, GK_])
    G = (B * A * B.transpose()).change_ring(ZZ)
    assert G.det() == 1
    assert all(G[i, i] % 2 == 0 for i in range(a + b))
    U = G.LLL_gram()
    G2 = U.transpose() * G * U
    B2 = U.transpose() * B                                                  # rows are the new basis / 行が新しい基底
    D = matrix.diagonal(QQ, [1] * a + [-1] * b)
    # sigma (B2 D = S B2 in the basis of row vectors) as the action g on column coordinates (g^T G g = G)
    # sigma（行ベクトルの基底で B2 D = S B2）を、列ベクトルの座標への作用 g（g^T G g = G）にする
    S = (B2 * D * B2.inverse())
    g = S.transpose().change_ring(ZZ)
    assert g.transpose() * G2 * g == G2 and g * g == 1
    return G2, g


def min_info(G):
    """The minimum, and the numbers of vectors of norm 2 and 4 (as counted by qfminim).

    最小ノルムと、ノルム 2, 4 の本数（qfminim の数え方）。"""
    r2 = int(pari(G).qfminim(2, 0, 0)[0])
    if r2:
        return 2, r2, None
    r4 = int(pari(G).qfminim(4, 0, 0)[0])
    return 4, 0, r4


def load(pair):
    """The components of a pair and its list of classes <pair>.json in LISTS.

    組の成分と、LISTS の中の類の一覧 <組>.json。"""
    nF, nK = gg.PAIRS[pair]
    GF_, F, gF = gg.part(nF); GK_, K, gK = gg.part(nK)
    res = json.load(open(Path(LISTS) / f"{pair}.json"))
    return GF_, GK_, F, K, gF, gK, res


def aut_tool(G, tag, threads=2):
    """|Aut L| and generators (integer matrices with g^T G g = G) by rank32.bliss_automorphisms; its temporary files go to OUT/bliss-scratch.

    rank32.bliss_automorphisms で Aut の生成元を得る（g^T G g = G の整数行列）。一時ファイルは OUT/bliss-scratch に置く。"""
    m = bliss_automorphisms
    m.OUT = Path(OUT) / "bliss-scratch"; m.OUT.mkdir(parents=True, exist_ok=True)
    return m.automorphisms(m.Shell(G, tag, threads))


def centralizer_order(gens, sigma, G=None):
    """The order of the matrix group and of the centralizer of sigma in it (through a permutation group on orbits of short vectors, GAP).

    With G, it starts from shortest vectors of full rank instead of the basis (the orbits stay among
    the shortest vectors); basis vectors can be long, and with a large Aut the orbit then has nearly
    |Aut| points (more than 5 GB for |Stab| = 16384).

    行列群とその中の sigma の中心化群の位数（短いベクトルの軌道で置換群に直して GAP で求める）。

    G を渡すと、基底の代わりに最短ベクトルのうち階数が満ちるものから出発する（軌道は最短ベクトルの中に
    収まる）。基底のベクトルは長いことがあり、Aut が大きいと軌道の点の数が |Aut| に近づく
    （|Stab| = 16384 で 5 GB を超えた）。"""
    n = sigma.nrows()
    gl = [matrix(ZZ, g) for g in gens]
    # permutation group on the union of the orbits of the starting points / 出発点の軌道の和集合の上の置換群
    pts = {}; order = []
    if G is None:
        start = [tuple(1 if j == i else 0 for j in range(n)) for i in range(n)]
    else:
        M = qfsafe.short_vectors(G, 4)
        rows = matrix(GF(1000003), M[:20000].tolist()).pivot_rows()
        assert len(rows) == n, ("the shortest vectors do not have full rank", len(rows))
        start = [tuple(int(x) for x in M[i]) for i in rows]
    for e in start:
        for sg in (1, -1):
            k = tuple(sg * x for x in e)
            if k not in pts:
                pts[k] = len(order); order.append(k)
    q = list(order); arrs = [np.array(g, dtype=np.int64) for g in gl + [sigma]]
    while q:
        Q = np.array(q, dtype=np.int64); nx = []
        for A in arrs:
            im = Q @ A.T
            for r in map(tuple, im.tolist()):
                if r not in pts:
                    pts[r] = len(order); order.append(r); nx.append(r)
        q = nx
    V = np.array(order, dtype=np.int64)
    def perm(A):
        return libgap.PermList([pts[r] + 1 for r in map(tuple, (V @ A.T).tolist())])
    H = libgap.Group([perm(a) for a in arrs[:-1]])
    s = perm(arrs[-1])
    assert bool(libgap.IN(s, H)), "sigma is not in Aut"
    return int(libgap.Size(H)), int(libgap.Size(libgap.Centralizer(H, s)))


def short_profile(G):
    """Isometry invariants: the number of vectors of norm <= 4, the table of inner products of the lines of norm 4 (distribution per row), and the sizes of the components of the graph of non-orthogonality.

    等長で変わらない量：ノルム <= 4 の本数と、ノルム 4 の直線の内積の表（行ごとの分布）、非直交のグラフの連結成分の大きさ。"""
    V = qfsafe.short_vectors(G, 4); W = V @ np.array(G.rows(), dtype=np.int64)
    N = (W * V).sum(axis=1); V4 = V[N == 4]; W4 = W[N == 4]
    if len(V4) > 20000:
        return len(V)
    IP = np.abs(W4 @ V4.T)
    rows = Counter(tuple(np.bincount(r, minlength=5).tolist()) for r in IP)
    seen = np.zeros(len(V4), dtype=bool); comps = []
    for i in range(len(V4)):
        if seen[i]:
            continue
        seen[i] = True; front = np.array([i]); size = 1
        while len(front):
            nb = (IP[front] != 0).any(axis=0) & ~seen
            seen |= nb; front = np.nonzero(nb)[0]; size += len(front)
        comps.append(size)
    return len(V), sorted(rows.items()), sorted(comps)


def glue_of(Gm, sg, GF0, GK0):
    """phi in the bases of the components (GF0, GK0) from a lattice Gm and an involution sg (g^T G g = G); None if the components are not isometric.

    格子 Gm と対合 sg（g^T G g = G）から、成分 (GF0, GK0) の基底での phi を取り出す。成分と等長でなければ None。"""
    n = Gm.nrows()
    fb = (sg - 1).right_kernel().basis_matrix().transpose(); kb = (sg + 1).right_kernel().basis_matrix().transpose()
    GFl = (fb.transpose() * Gm * fb).change_ring(ZZ); GKl = (kb.transpose() * Gm * kb).change_ring(ZZ)
    if GFl.nrows() != GF0.nrows():
        return None
    # different invariants mean not isometric; look at them before qfisom (for pairs with huge groups that are not
    # isomorphic, such as sqrt(2) E_8^2 and sqrt(2) D_16^+, qfisom did not return in 40 minutes)
    # 不変量が違えば等長でない。qfisom の前に見る（sqrt(2) E_8^2 と sqrt(2) D_16^+ のように群が巨大で同型でない対は qfisom が 40 分経っても返らなかった）
    if short_profile(GF0) != short_profile(GFl) or short_profile(GK0) != short_profile(GKl):
        return None
    S = qfsafe.qfisom(GF0, GFl)
    if S is None:
        return None
    T = qfsafe.qfisom(GK0, GKl)
    if T is None:
        return None
    _, JF, SF = disc_gens(GF0); _, JK, SK = disc_gens(GK0)
    Si = S.inverse(); Ti = T.inverse(); Fi = GFl.inverse(); Ki = GKl.inverse()
    pairs = []
    for t in range(n):
        v = vector(QQ, [1 if j == t else 0 for j in range(n)])
        vF = (v + sg * v) / 2; vK = (v - sg * v) / 2
        cF = Fi * (fb.transpose() * Gm * vF); cK = Ki * (kb.transpose() * Gm * vK)
        pairs.append((mask_of(Si * cF, JF, SF), mask_of(Ti * cK, JK, SK)))
    s = len(JF); basis = {}
    for x, y in pairs:
        while x:
            h = x.bit_length() - 1
            if h in basis:
                x ^= basis[h][0]; y ^= basis[h][1]
            else:
                basis[h] = (x, y); break
    assert len(basis) == s
    phi = [0] * s
    for i in range(s):
        x, y = 1 << i, 0
        while x:
            h = x.bit_length() - 1
            x ^= basis[h][0]; y ^= basis[h][1]
        phi[i] = y
    return np.array(phi, dtype=np.int64)


def class_table(pair, g, res, kFK):
    """Key -> record. For a pair already rechecked, the normal forms are not recomputed when the set of keys of the checkpoint has the keys_sha256 of the recheck.

    キー -> 記録。recheck の済んだ組は、チェックポイントのキーの集合が recheck の keys_sha256 と同じなら正規形を計算し直さない。"""
    ck = Path(LISTS) / f"{pair}.ckpt.json"; vr = Path(OUT) / f"verify-recheck-{pair}.json"
    if ck.exists() and vr.exists():
        t = gg._load_checkpoint(ck)[0]
        if len(t) == res["n_classes"] and hashlib.sha256(b"".join(sorted(t))).hexdigest() == json.load(open(vr))["keys_sha256"]:
            print(pair, "keys are the set of the recheck (from the checkpoint)", flush=True)
            return {k: dict(stab=c["stab"] * kFK, extremal=c["ext"]) for k, c in t.items()}
    cl = {}
    for c in res["classes"]:
        k, n = g.canon(np.array(c["phi"], dtype=np.int64)); cl[k] = c
        assert n * kFK == c["stab"]
    return cl


def element_order(w):
    """The order of an integer matrix of finite order.

    有限位数の整数行列の位数。"""
    I = matrix.identity(ZZ, w.nrows()); x = w; o = 1
    while x != I:
        x = x * w; o += 1
        assert o < 100000
    return o
