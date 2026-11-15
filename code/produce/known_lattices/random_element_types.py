"""Types (p, f, s) of automorphisms of odd prime order of the known lattices, from random elements of their groups.

Paper: computation C32 of the first paper (the known lattices). For each lattice with a generating set
of Aut(L), random products of the generators (random.Random(1), restarted for each lattice) give
elements g; for every odd prime p dividing the order of g, x = g^(o/p) has prime order p and its type
(p, f, s) is recorded (f = rank of the fixed lattice, s = p-rank of its discriminant group, that is its
dimension over F_p). It also counts how often g^(o/4) squares to -1, and for a fixed-point-free element
of order 3 the size of its conjugacy class and its centralizer.
A type that is found is certain (x^T G x = G and the order are checked in integers), but a type that is
not found is not proved absent.
Usage:
  python3 random_element_types.py [TRIALS] [--known FILE] [--code-grams FILE] [--kmax24 FILE] [--reps-dir DIR] [--out FILE]
    TRIALS default 600 (the run of C32 used 800).
Input (by default in computation-outputs/C32 of the release assets):
        --known FILE       JSON Lines: id, order, generators (the groups of the known lattices)
        --code-grams FILE  lattices from codes (JSON object NAME -> {gram})
        --kmax24 FILE      lattices of the stratum k_max = 24 (JSON Lines: idx, gram)
        --reps-dir DIR     representatives reps-*.jsonl (JSON Lines: id, gram)
Output: --out FILE         (default results/C32/known32-types.json) JSON: id -> {id, order,
                           square_minus_one_seen, types [[p, f, s, seen], ...], fpf3}
Needs: Sage and numpy. Time: about 7 minutes and 2.7 GB with 800 trials.

既知の格子の群のランダムな元から、奇素数位数の自己同型の型 (p, f, s) を集める。

論文：第 1 の論文の計算 C32（既知の格子）。Aut(L) の生成元が分かっている各格子で、生成元のランダムな積
（random.Random(1)。格子ごとに初めから）から元 g を作る。g の位数を割る奇素数 p ごとに、
x = g^(o/p) は位数 p で、その型 (p, f, s) を記録する（f = 固定格子の階数、s = その判別式群の
p-ランク（F_p 上の次元））。
g^(o/4) の 2 乗が -1 になった回数も数え、固定点を持たない位数 3 の元については共役類と中心化群の大きさを求める。
見つかった型は確か（整数行列で x^T G x = G と位数を検算する）だが、見つからなかったことは不在の証明ではない。
使い方：
  python3 random_element_types.py [TRIALS] [--known FILE] [--code-grams FILE] [--kmax24 FILE] [--reps-dir DIR] [--out FILE]
    TRIALS の既定は 600（C32 の実行では 800）。
入力（既定は Release のアセットの computation-outputs/C32）：
      --known FILE       JSON Lines：id、order、generators（既知の格子の群）
      --code-grams FILE  符号から作った格子（NAME -> {gram} の JSON オブジェクト）
      --kmax24 FILE      k_max = 24 の層の格子（JSON Lines：idx、gram）
      --reps-dir DIR     代表 reps-*.jsonl（JSON Lines：id、gram）
出力：--out FILE         （既定は results/C32/known32-types.json）JSON：id -> {id, order,
                         square_minus_one_seen, types [[p, f, s, seen], ...], fpf3}
必要なもの：Sage と numpy。時間：800 回で約 7 分・2.7 GB。
"""
import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix  # noqa: E402

from rank32.io import asset  # noqa: E402
from rank32.random_element_types import I, class_size, fixed_type, grams, mpow  # noqa: E402


def main(trials, known, Gs, out_path):
    out = {}
    for line in Path(known).open():
        rec = json.loads(line)
        rid = rec["id"]
        G = matrix(ZZ, Gs[rid])
        gens = [np.array(g, dtype=np.int64) for g in rec["generators"]]
        order = int(rec["order"])
        rng = random.Random(1)
        pool = [g.copy() for g in gens] * 3
        for _ in range(200):
            i, j = rng.randrange(len(pool)), rng.randrange(len(pool))
            if i != j:
                pool[i] = pool[i] @ pool[j]
        types = {}
        fpf = {}
        square = 0
        for t in range(trials):
            i, j = rng.randrange(len(pool)), rng.randrange(len(pool))
            if i != j:
                pool[i] = pool[i] @ pool[j]
            g = pool[i]
            o, h = 1, g
            while not np.array_equal(h, I):
                h = h @ g
                o += 1
            if o % 4 == 0:
                x = mpow(g, o // 4)
                if np.array_equal(x @ x, -I):
                    square += 1
            for p in (3, 5, 7, 11, 13, 17, 31):
                if o % p == 0:
                    x = mpow(g, o // p)
                    tr = int(np.trace(x))
                    key0 = (p, tr)
                    if key0 not in types:
                        f, s = fixed_type(G, x, p)
                        assert (matrix(ZZ, x.tolist()).T * G * matrix(ZZ, x.tolist()) == G)
                        types[key0] = dict(p=p, f=f, s=s, seen=0)
                        if p == 3 and f == 0:
                            fpf[x.tobytes()] = x
                    types[key0]["seen"] += 1
        row = dict(id=rid, order=str(order), square_minus_one_seen=square, types=sorted(([v["p"], v["f"], v["s"], v["seen"]] for v in types.values())))
        if fpf:
            inv = [np.rint(np.linalg.inv(g)).astype(np.int64) for g in gens]
            x = next(iter(fpf.values()))
            c = class_size(x, gens, inv)
            row["fpf3"] = dict(class_size=c, centralizer=(str(order // c) if c else None), inverse_in_class=None)
        out[rid] = row
        print(row, flush=True)
    Path(out_path).write_text(json.dumps(out, indent=0))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Types of elements of odd prime order of known lattices (C32) / 既知の格子の奇素数位数の元の型（C32）")
    ap.add_argument("trials", type=int, nargs="?", default=600)
    ap.add_argument("--known", type=Path, default=asset("computation-outputs/C32/known32.jsonl"), help="groups of the known lattices / 既知の格子の群")
    ap.add_argument("--code-grams", type=Path, default=asset("computation-outputs/C32/code-lattices.json"), help="lattices from codes / 符号から作った格子")
    ap.add_argument("--kmax24", type=Path, default=asset("computation-outputs/C32/kmax24-stratum.jsonl"), help="stratum k_max = 24 / k_max = 24 の層")
    ap.add_argument("--reps-dir", type=Path, default=asset("computation-outputs/C32"), help="directory of reps-*.jsonl / reps-*.jsonl の場所")
    ap.add_argument("--out", type=Path, default=Path("results/C32/known32-types.json"), help="output / 出力")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    main(a.trials, a.known, grams(a.code_grams, a.kmax24, a.reps_dir), a.out)
