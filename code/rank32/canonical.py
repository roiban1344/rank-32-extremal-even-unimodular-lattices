"""Automorphism groups and isometry tests of lattices (and of Hermitian structures) through a graph of short vectors and bliss.

PARI qfauto / qfisom stall on some lattices of rank 16 to 32; this is the replacement.

S = all vectors of norm <= m, where m is the least value for which S generates the lattice
(an invariant of the lattice). The vertices are S and the colour of the edge (x, y) is the inner
product (x, y) (with omega, the pair ((x, y), (x, omega y))). Since S generates the lattice,
the automorphisms of the graph correspond one to one to the isometries of the lattice (with
omega, to those that commute with omega).
Isometry: if the hashes of the canonical forms differ, the lattices are not isometric. If they
agree, we build the map and check that it is an integral isometry.
For a large S the colours with small absolute inner product are left out (keep_min). The graph
may then have more automorphisms, so every generator is checked.

Row-vector conventions: basis b_i, x = c b, omega acts as c -> c W. An isometry is c -> c R
with R G R^T = G.
The C program _canonical.c is compiled on first use (gcc -O2, linked to bliss) into the
directory $RANK32_BUILD (default $XDG_CACHE_HOME/rank32, that is ~/.cache/rank32). bliss is looked
up under $BLISS_PREFIX, else $CONDA_PREFIX, else the prefix of the running Python. The temporary
graph files go to $RANK32_TMP if it is set, else to the system default.

格子（とエルミート構造）の自己同型群と同型判定を、短いベクトルのグラフと bliss で求める。

PARI の qfauto / qfisom は、階数 16 から 32 の一部の格子で止まる。その代わりである。

S = ノルム <= m のベクトル全部とする。m は S が格子を生成する最小の値で、格子の不変量である。
頂点を S、辺 (x, y) の色を内積 (x, y)（omega があるときは組 ((x, y), (x, omega y))）とする。
S が格子を生成するので、グラフの自己同型は格子の等長写像と 1 対 1 に対応する（omega があれば、
omega と可換なものと）。
同型：標準形のハッシュが違えば同型でない。同じなら写像を作り、整数の等長写像であることを検算する。
大きい S では、内積の絶対値が小さい色を省く（keep_min）。そのときはグラフの自己同型の方が
大きくなりうるので、生成元をすべて検算する。

行ベクトルの約束：基底 b_i、x = c b、omega は c -> c W。等長写像は c -> c R（R G R^T = G）。
C のプログラム _canonical.c は、初めて使うときに $RANK32_BUILD（既定は $XDG_CACHE_HOME/rank32、つまり ~/.cache/rank32）へ
コンパイルする（gcc -O2、bliss をリンクする）。bliss は $BLISS_PREFIX、なければ $CONDA_PREFIX、
なければ実行中の Python の prefix の下で探す。グラフの一時ファイルは、$RANK32_TMP があればその下、
なければシステムの既定の場所に置く。
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from sage.all import ZZ, QQ, matrix, pari, PermutationGroup

from .forms import short_vectors

_SRC = Path(__file__).with_name("_canonical.c")


def _build_dir():
    """The directory for compiled helpers.

    コンパイルしたものを置くディレクトリ。"""
    d = os.environ.get("RANK32_BUILD")
    if d:
        return Path(d)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "rank32"


def _bliss_prefix():
    return Path(os.environ.get("BLISS_PREFIX") or os.environ.get("CONDA_PREFIX") or sys.prefix)


def _binary():
    env = _bliss_prefix()
    build = _build_dir()
    exe = build / ("canonical-" + hashlib.sha256(_SRC.read_bytes()).hexdigest()[:12])
    if not exe.exists():
        build.mkdir(parents=True, exist_ok=True)
        tmp = exe.with_name(exe.name + f".{os.getpid()}.tmp")
        subprocess.run(["gcc", "-O2", str(_SRC), "-o", str(tmp), "-I" + str(env / "include"), "-L" + str(env / "lib"),
                        "-Wl,-rpath," + str(env / "lib"), "-lbliss", "-lstdc++"], check=True)
        os.replace(tmp, exe)
    return exe


def generating_shell(G, maxvec=60000):
    """(m, S): S is all vectors of norm <= m (both signs, coordinates as rows), for the least m such that S generates Z^n; (m, None) if S would have more than maxvec rows.

    (m, S)：S はノルム <= m のベクトル全部（+- 両方、行が座標）で、Z^n を生成する最小の m。S の行が maxvec を超えるなら (m, None)。"""
    n = G.nrows()
    q = pari(G)
    m = 2
    while True:
        # qfminim counts v and -v separately (240 for E_8 and m = 2), so cnt is the number of rows S would have
        # qfminim は v と -v を別に数える（E_8、m = 2 で 240）。だから cnt は S の行の数になる
        cnt = int(q.qfminim(m, 0)[0])
        if cnt * 2 > maxvec * 2:
            return m, None
        if cnt >= n:
            half = short_vectors(G, m)
            if np.linalg.matrix_rank(half.astype(float)) == n:
                H = matrix(ZZ, half.tolist()).hermite_form(include_zero_rows=False)
                if abs(H.det()) == 1:
                    Gn = np.array(G.list(), dtype=np.int64).reshape(n, n)
                    norms = np.einsum("ij,jk,ik->i", half, Gn, half)
                    half = half[norms > 0]
                    return m, np.vstack([half, -half])
        m += 2


def _edges(S, G, W, keep_min):
    n = G.nrows()
    Gn = np.array(G.list(), dtype=np.int64).reshape(n, n)
    A = S @ Gn
    if W is not None:
        Wn = np.array(W.list(), dtype=np.int64).reshape(n, n)
        SW = (S @ Wn).T  # (x, omega y) is built row by row (the full table is large) / (x, omega y) は行ごとに作る（表全体は大きい）
    chunks = []
    for i in range(S.shape[0]):
        a = A[i] @ S.T
        if W is None:
            c = (a + 64) * 256
        else:
            c = (a + 64) * 256 + (A[i] @ SW + 64)
        mask = np.abs(a) >= keep_min
        mask[i] = True  # include the colour of the vertex itself (its norm) / 自分自身の色（ノルム）を含める
        js = np.nonzero(mask)[0]
        chunks.append(np.column_stack([np.full(js.shape, i, dtype=np.int32), js.astype(np.int32), c[js].astype(np.int32)]))
    return np.vstack(chunks)


def canonical(G, W=None, maxvec=60000, keep_min=None, timeout=None):
    """Return dict(m, S, order, hash, gens (permutations of S), lab), or None if the graph cannot be built.

    dict(m, S, order, hash, gens（S の置換の配列）, lab) を返す。作れなければ None。"""
    m, S = generating_shell(G, maxvec)
    if S is None:
        return None
    if keep_min is None:
        # dense graphs take several GB in bliss (more than 5 GB for a lattice of the genus (20, 6), the even rank-20
        # lattices F with F^#/F = F_3^6 of Appendix D.4.6 of the first paper); for a lattice Z of rank z = 13 over
        # Z[omega] (type 3-(13, 6)-5, §3.3 of the first paper; 20000 vectors) keep 3 gave the same order as PARI
        # and keep 4 gave extra automorphisms
        # 密なグラフは bliss が数 GB になる（種 (20, 6) の格子で 5 GB 超。第 1 の論文の付録 D.4.6 の、F^#/F = F_3^6 で
        # 階数 20 の偶格子 F）。Z[omega] 上の階数 z = 13 の格子 Z（型 3-(13, 6)-5、第 1 の論文の §3.3。2 万本）で、
        # keep 3 は位数が PARI と一致し、keep 4 は余分な自己同型が出た
        keep_min = 0 if S.shape[0] <= 2000 else (2 if S.shape[0] <= 8000 else 3)
    E = _edges(S, G, W, keep_min)
    with tempfile.TemporaryDirectory(dir=os.environ.get("RANK32_TMP")) as td:
        inp = Path(td) / "g.bin"
        out = Path(td) / "g.out"
        with open(inp, "wb") as f:
            f.write(np.array([S.shape[0], E.shape[0]], dtype=np.int32).tobytes())
            f.write(E.tobytes())
        r = subprocess.run([str(_binary()), str(inp), str(out)], capture_output=True, text=True, check=True, timeout=timeout)
        info = json.loads(r.stdout)
        gens = np.fromfile(str(out) + ".gens", dtype=np.int32).reshape(-1, S.shape[0])
        lab = np.fromfile(str(out) + ".lab", dtype=np.int32)
    res = dict(m=m, S=S, gens=gens, lab=lab, hash=info["hash"], nodes=info["nodes"], keep_min=keep_min, approx=info["group_size"])
    check_generators(res, G, W)
    res["order"] = group_order(res)
    return res


def shell_hash(G, W=None, bound=4, timeout=1800, keep_min=None):
    """Hash of the canonical form of the configuration of the vectors of norm <= bound; different values mean non-isometric (equal values do not prove isometry).

    S need not generate the lattice. The rule for leaving out colours depends only on |S|, so
    isometric lattices give the same value.

    ノルム <= bound のベクトルの配置の標準形のハッシュ。違えば同型でない（同じでも同型とは限らない）。

    S が格子を生成しなくてもよい。省く色の規則は |S| だけで決まるので、同型な格子では同じ値になる。"""
    half = short_vectors(G, bound)
    S = np.vstack([half, -half])
    if keep_min is None:
        keep_min = 0 if S.shape[0] <= 2000 else (2 if S.shape[0] <= 8000 else 3)
    E = _edges(S, G, W, keep_min)
    with tempfile.TemporaryDirectory(dir=os.environ.get("RANK32_TMP")) as td:
        inp = Path(td) / "g.bin"
        with open(inp, "wb") as f:
            f.write(np.array([S.shape[0], E.shape[0]], dtype=np.int32).tobytes())
            f.write(E.tobytes())
        r = subprocess.run([str(_binary()), str(inp), str(Path(td) / "g.out")], capture_output=True, text=True, check=True, timeout=timeout)
    return json.loads(r.stdout)["hash"]


def _lift(S1, perm_image, G1, G2):
    """The rational matrix R (c -> c R) that sends the rows of S1 to perm_image (rows of S2), or None unless it is an integral isometry.

    S1 の行を perm_image（S2 の行）へ移す有理行列 R（c -> c R）を作り、整数の等長写像でなければ None。"""
    n = G1.nrows()
    idx = []
    rows = []
    # choose n independent rows / 独立な n 本を選ぶ
    for i in range(S1.shape[0]):
        c = matrix(ZZ, rows + [list(map(int, S1[i]))])
        if c.rank() == len(rows) + 1:
            rows.append(list(map(int, S1[i])))
            idx.append(i)
        if len(rows) == n:
            break
    A = matrix(QQ, rows)
    B = matrix(QQ, [list(map(int, perm_image[i])) for i in idx])
    R = A.inverse() * B
    if R.denominator() != 1:
        return None
    R = R.change_ring(ZZ)
    if R * G2 * R.transpose() != G1:
        return None
    if abs(R.det()) != 1:
        return None
    if not np.array_equal(S1 @ np.array(R.list(), dtype=np.int64).reshape(n, n), perm_image):
        return None
    return R


def check_generators(res, G, W=None):
    """Check that every generator is an isometry of the lattice (commuting with omega); raise on failure.

    生成元がすべて格子の等長写像（omega と可換）であることを検算する。失敗は例外。"""
    S = res["S"]
    for p in res["gens"]:
        R = _lift(S, S[p], G, G)
        if R is None:
            raise ValueError("an automorphism of the graph is not an isometry of the lattice (lower keep_min)")
        if W is not None and R * W != W * R:
            raise ValueError("does not commute with omega")


def generators(res, G):
    """The generators as integer matrices g in the column convention (g^T G g = G), as rank32.forms.qfauto returns them.

    生成元を列の約束の整数行列 g（g^T G g = G）で返す（rank32.forms.qfauto と同じ約束）。"""
    S = res["S"]
    return [_lift(S, S[p], G, G).transpose() for p in res["gens"]]


def group_order(res):
    """The order of the group generated by res["gens"]: the value of bliss when it is below 2^52, else computed by Sage.

    res["gens"] の生成する群の位数。bliss の値が 2^52 未満ならそれを使い、それ以上なら Sage で求める。"""
    gens = res["gens"]
    if len(gens) == 0:
        return ZZ(1)
    approx = int(res["approx"])
    if approx < 2 ** 52:
        return ZZ(approx)
    return ZZ(PermutationGroup([[int(x) + 1 for x in p] for p in gens], canonicalize=False).order())


def isometry(res1, G1, res2, G2, W1=None, W2=None):
    """Compare canonical forms: an R with R G2 R^T = G1 (c -> c R maps L1 to L2) if isometric, else None.

    res1 and res2 must come from canonical() with the same keep_min. The default rule chooses keep_min from |S| only,
    so it agrees for isometric lattices; with different values the hashes can differ, so None may be returned
    for isometric lattices.

    標準形で比べる。同型なら R（R G2 R^T = G1 の行列、c -> c R が L1 -> L2）、違えば None。
    res1 と res2 は、同じ keep_min の canonical() の結果でなければならない。既定の規則は keep_min を |S| だけで決めるので、
    同型な格子では一致する。値が違うと、同型な格子でもハッシュが違うことがあり、None を返しうる。"""
    if res1["hash"] != res2["hash"] or res1["S"].shape != res2["S"].shape or res1["m"] != res2["m"]:
        return None
    inv2 = np.argsort(res2["lab"])
    image = res2["S"][inv2[res1["lab"]]]
    R = _lift(res1["S"], image, G1, G2)
    if R is None:
        raise ValueError("the hashes agree but the map is not an isometry (lower keep_min)")
    if W1 is not None and R * W2 != W1 * R:
        raise ValueError("the hashes agree but the map does not commute with omega")
    return R
