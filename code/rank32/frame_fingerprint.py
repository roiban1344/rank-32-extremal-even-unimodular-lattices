"""The frame distribution a_k of an extremal even unimodular lattice of rank 32, as a fingerprint.

Paper: an auxiliary computation, not in Table 11 of the first paper (the paper does not use the frame
distributions). The fingerprints are an invariant for comparing lattices; the checks of the tools (C37)
use them on the lattices of a random walk. The entry point with the command is
code/toolcheck/frame_fingerprint.py.
a_k = number of admissible classes of L/2L (nonzero classes whose norms are divisible by 4) that
contain exactly k pairs +-z of vectors of norm 8. Vectors of norm 8 in the same class are mutually
orthogonal, and two orthogonal ones z, z' in one class give the minimal vectors (z +- z')/2.
So we form the classes mod 2, b_x (32 bits), of the 73440 lines of minimal vectors and count the
XORs b_x ^ b_x' of all pairs {x, x'}: this gives the number m(c) = k(k-1)/2 of orthogonal pairs in
each class (the C program _frame_hash.c). The classes with k = 1 are invisible; a_1 follows from
sum k a_k = number of lines of norm 8, and a_0 from the number of admissible classes.
The fingerprint fp is the first 16 hexadecimal digits of the sha256 of the sorted list of (k, a_k).
Lattices are looked up by name in GRAMS (filled by load_lattices); NAME@uS applies the random
unimodular change of basis with seed S (a positive control).
The C program is compiled on first use ($CC, default gcc: -O3 -march=native -fopenmp -DPL=$FRAME_PL, default 12)
into $RANK32_BUILD (default ~/.cache/rank32); its binary input goes to $FRAME_SCRATCH or to the system
temporary directory. One lattice takes about 2 GB and one minute with 4 threads.

階数 32 の極値偶ユニモジュラー格子の「枠の分布」a_k を、格子の指紋として求める。

論文：第 1 の論文の表 11 に無い補助の計算（論文は枠の分布を使わない）。指紋は格子を比べるための不変量で、
道具の検査（C37）がランダムウォークの格子で使う。コマンドの入口は code/toolcheck/frame_fingerprint.py である。
a_k = L/2L の許容類（ノルムが 4 で割れる 0 でない類）のうち、ノルム 8 のベクトルの対 +-z がちょうど
k 対入るものの個数。同じ類のノルム 8 のベクトルは互いに直交し、同じ類の直交する 2 本 z、z' から
最短ベクトル (z +- z')/2 が出る。そこで最短ベクトルの直線 73440 本の mod 2 の類 b_x（32 ビット）を作り、
全対 {x, x'} の XOR b_x ^ b_x' を数えると、各類に入る直交対の数 m(c) = k(k-1)/2 が出る（C のプログラム
_frame_hash.c）。k = 1 の類は見えないので、a_1 は sum k a_k = ノルム 8 の直線の数から、a_0 は許容類の数から出す。
指紋 fp は、(k, a_k) を並べたリストの sha256 の先頭 16 桁（16 進）である。
格子は名前で GRAMS（load_lattices で読む）から引く。NAME@uS は、種 S の乱択ユニモジュラー基底変換を
かけたもの（陽性対照）である。
C のプログラムは、初めて使うときに $RANK32_BUILD（既定は ~/.cache/rank32）へコンパイルする
（$CC、既定は gcc：-O3 -march=native -fopenmp -DPL=$FRAME_PL、既定 12）。その入力のバイナリは $FRAME_SCRATCH か、
システムの一時ディレクトリに置く。1 つの格子で約 2 GB、4 スレッドで約 1 分かかる。
"""
import hashlib
import json
import os
import random
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

CSRC = Path(__file__).with_name("_frame_hash.c")
A_ADM = 2**31 + 2**15 - 1                       # admissible classes (nonzero, norm = 0 mod 4) / 許容類の数
N8_LINES = 4844836800 // 2                      # lines of norm 8 (half the coefficient of the extremal theta series) / ノルム 8 の直線の数
GRAMS = {}                                      # name -> Gram matrix (list of rows) / 名前 -> Gram 行列


def _build_dir():
    d = os.environ.get("RANK32_BUILD")
    if d:
        return Path(d)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "rank32"


def _scratch():
    return Path(os.environ.get("FRAME_SCRATCH") or tempfile.gettempdir())


def binary():
    """The compiled C program for the split PL = $FRAME_PL (default 12; the table has 2^(32-PL) entries).

    分け方 PL = $FRAME_PL（既定 12。表の大きさは 2^(32-PL)）のための、コンパイルした C のプログラム。"""
    pl = os.environ.get("FRAME_PL", "12")
    build = _build_dir()
    exe = build / f"frame_hash-pl{pl}-{hashlib.sha256(CSRC.read_bytes()).hexdigest()[:12]}"
    if not exe.exists():
        build.mkdir(parents=True, exist_ok=True)
        cc = os.environ.get("CC", "gcc")
        tmp = exe.with_name(exe.name + f".{os.getpid()}.tmp")
        subprocess.run([cc, "-O3", "-march=native", "-fopenmp", f"-DPL={pl}", "-o", str(tmp), str(CSRC)], check=True)
        os.replace(tmp, exe)
    return exe


def load_lattices(path, prefix=""):
    """Add lattices to GRAMS from a JSON Lines file of rows (id, gram), or from a JSON object name -> {gram}.

    格子を GRAMS に足す。(id, gram) の行の JSON Lines か、名前 -> {gram} の JSON オブジェクトから読む。"""
    text = Path(path).read_text()
    try:
        obj = json.loads(text)
    except json.JSONDecodeError:
        obj = None
    if isinstance(obj, dict):
        for k, v in obj.items():
            GRAMS[prefix + k] = v["gram"]
    else:
        for line in text.splitlines():
            if line.strip():
                r = json.loads(line)
                GRAMS[prefix + r["id"]] = r["gram"]


def random_rebase(G, seed):
    """G -> U G U^t for a random unimodular U (400 elementary operations and a permutation of rows, seeded).

    乱択ユニモジュラー U（基本変形 400 回と行の置換、種つき）で G -> U G U^t。"""
    from sage.all import ZZ, matrix
    rng = random.Random(int(seed))
    n = G.nrows()
    U = matrix.identity(ZZ, n)
    for _ in range(400):
        i, j = rng.sample(range(n), 2)
        U.add_multiple_of_row(i, j, rng.choice([-2, -1, 1, 2]))
    p = list(range(n))
    rng.shuffle(p)
    U = matrix(ZZ, [U.row(k) for k in p])
    assert abs(U.det()) == 1
    return U * G * U.transpose()


def gram(name):
    """The Gram matrix of a lattice by name: an entry of GRAMS, or NAME@uS (random change of basis with seed S).

    名前から格子の Gram 行列を返す。GRAMS の項目か、NAME@uS（種 S の乱択基底変換）。"""
    from sage.all import ZZ, matrix
    base, _, u = name.partition("@u")
    G = matrix(ZZ, GRAMS[base])
    if u:
        G = random_rebase(G, u)
    return G


def min_lines(G):
    """All minimal vectors (one of each pair +-x) by qfminim, with the norms checked in integers.

    最小ベクトル（+-x の片方）を qfminim で全部。ノルムは整数で検算する。"""
    from sage.all import pari
    r = pari(G).qfminim(4)  # all of them (flag 0); norms are checked below / 全部（flag 0）。ノルムは下で検算
    assert int(r[0]) == 146880, int(r[0])
    rows = str(r[2].mattranspose())[1:-1].split(";")
    X = np.array([[int(c) for c in x.split(",")] for x in rows], dtype=np.int64)  # (73440, 32)
    Gi = np.array([[int(x) for x in row] for row in G.rows()], dtype=np.int64)
    nn = ((X @ Gi) * X).sum(axis=1)
    assert X.shape == (73440, 32) and (nn == 4).all()
    return X, Gi


def q_tables(Gi):
    """Byte tables of q(c) = (v,v)/2 mod 2 (v = sum c_i e_i).

    q(c) = (v,v)/2 mod 2（v = sum c_i e_i）をバイトごとの表に分ける。"""
    qb = np.zeros((4, 256), dtype=np.uint8)
    bb = np.zeros((6, 256, 256), dtype=np.uint8)
    bits = ((np.arange(256)[:, None] >> np.arange(8)) & 1).astype(np.int64)        # (256, 8)
    for a in range(4):
        Gaa = Gi[8 * a:8 * a + 8, 8 * a:8 * a + 8]
        diag = (np.diag(Gaa) // 2) % 2
        up = np.triu(Gaa, 1) % 2
        qb[a] = ((bits @ diag) + np.einsum("ui,ij,uj->u", bits, up, bits)) % 2
    for t, (a, c) in enumerate([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]):
        Gac = Gi[8 * a:8 * a + 8, 8 * c:8 * c + 8] % 2
        bb[t] = (bits @ Gac @ bits.T) % 2
    return qb, bb


def run_one(name, threads=4):
    """Run the C program on the lattice gram(name); returns its JSON output with the wall time and the name.

    格子 gram(name) で C のプログラムを走らせ、その JSON の出力に経過時間と名前を足して返す。"""
    X, Gi = min_lines(gram(name))
    b = ((X & 1).astype(np.uint64) << np.arange(32, dtype=np.uint64)).sum(axis=1).astype(np.uint32)
    assert np.unique(b).shape[0] == X.shape[0]                       # the lines of minimal vectors are in distinct classes / 最短ベクトルの直線は類が全部違う
    qb, bb = q_tables(Gi)
    scr = _scratch()
    scr.mkdir(parents=True, exist_ok=True)
    path = scr / f"frame-{os.getpid()}-{name.replace('@', '_').replace(':', '_')}.bin"
    with open(path, "wb") as f:
        f.write(np.int32(X.shape[0]).tobytes()); f.write(b.tobytes()); f.write(qb.tobytes()); f.write(bb.tobytes())
    t0 = time.time()
    r = subprocess.run([str(binary()), str(path), "256", str(threads)], capture_output=True, text=True, check=True)
    res = json.loads(r.stdout); res["wall"] = round(time.time() - t0, 1); res["id"] = name
    path.unlink()
    return res


def derive(hist0, n_min):
    """a_k and the fingerprint from hist0. A class with k pairs of lines of norm 8 has m = k(k-1)/2 orthogonal pairs of
    minimal vectors (two orthogonal lines z, z' of norm 8 in one class give the minimal vectors (z +- z')/2, and
    conversely), so the classes with k = 1 are invisible; a_1 follows from sum k a_k = number of lines of norm 8.
    The classes of minimal vectors contain no vectors of norm 8 (they give the extremal neighbours).

    hist0 から a_k と指紋を出す。k 対のノルム 8 の直線を持つ類では、直交する最短ベクトルの対の数は m = k(k-1)/2
    （同じ類の直交するノルム 8 の直線 z、z' から (z +- z')/2 が最短ベクトル。逆も）。だから k = 1 の類は見えない。
    a_1 は sum k a_k = ノルム 8 の直線の数（強制）から逆算する。最短ベクトルの類はノルム 8 を含まない（極値の近傍を与える）。"""
    ak = {}
    for m, c in hist0:
        k = int(round((1 + (1 + 8 * m) ** 0.5) / 2))
        assert k * (k - 1) // 2 == m, m
        ak[k] = c
    ak[1] = N8_LINES - sum(k * c for k, c in ak.items())
    a0 = A_ADM - sum(ak.values())                                   # includes the classes of minimal vectors / 最短ベクトルの類を含む
    mom = [sum(k ** j * c for k, c in ak.items()) for j in range(7)]
    return dict(ak=dict(sorted(ak.items())), a0=a0, p=a0 / A_ADM, n_min=n_min, moments=mom,
                fp=hashlib.sha256(json.dumps(sorted(ak.items())).encode()).hexdigest()[:16])
