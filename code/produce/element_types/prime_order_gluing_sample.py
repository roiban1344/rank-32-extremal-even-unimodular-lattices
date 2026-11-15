"""For the types 3-(z, d)-s with a fixed lattice: the fraction of anti-isometries phi for which the glued lattice is extremal, and one example each.

Paper: computation C2 of the first paper (Table 11), a check of the gluing lemma (rank32.prime_order_gluing).
For each pair (F, Z) of lists of minimum >= 4, A = {classes of F^# with vectors of norm 2/3} and
B = {classes of Z^# with vectors of norm 4/3}; random anti-isometries phi (seed 20261030) are counted with
phi(A) and B disjoint. If A or B is empty, every phi gives an extremal lattice. One glued lattice L is built
for each pair and checked: determinant 1, even, and no vectors of norm 2 (a check of the argument).
For d = 8 and d = 6 the lattice F is sqrt(3) E8 and sqrt(3) E6^#.
Usage:  python3 prime_order_gluing_sample.py [SAMPLES]   (default 3000)
Input:  --input-dir DIR with genus-p3-n<d>-s<s>.json, unimod-z<z>.json and overlat-z<z>.json (outputs of C2)
Output: --out DIR, glue.json (one row per pair) and glue-lattices.jsonl (the lattices built)
Needs:  SageMath (PARI). Time: about 1 minute with 40 samples; the default (3000) was not timed.

固定格子のある型 3-(z, d)-s について、貼り合わせの格子が極値になる反等長写像 phi の割合を測り、
1 個ずつ実例を作る。

論文：第 1 の論文の計算 C2（表 11）。貼り合わせの補題（rank32.prime_order_gluing）の検算。
最小 >= 4 の一覧の組 (F, Z) ごとに、A = {ノルム 2/3 のベクトルを持つ F^# の類}、B = {ノルム 4/3 の
ベクトルを持つ Z^# の類} とし、乱択の反等長写像 phi（種 20261030）のうち phi(A) と B が交わらないものを
数える。A か B が空なら、どの phi でも極値になる。組ごとに貼り合わせの格子 L を 1 個作り、行列式 1・偶・
ノルム 2 のベクトルが無いことを確かめる（論法の検算）。
d = 8 と d = 6 では F は sqrt(3) E8 と sqrt(3) E6^# である。
使い方：python3 prime_order_gluing_sample.py [標本数]（既定 3000）
入力：--input-dir DIR の genus-p3-n<d>-s<s>.json、unimod-z<z>.json、overlat-z<z>.json（C2 の出力）
出力：--out DIR の glue.json（組ごとに 1 行）と glue-lattices.jsonl（作った格子）
必要なもの：SageMath（PARI）。時間：標本数 40 で約 1 分。既定（3000）は計っていない。
"""
import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix, pari  # noqa: E402

from rank32 import prime_order_gluing as pg  # noqa: E402
from rank32.io import asset, write_json_atomic_indent0  # noqa: E402

DEFAULT_INPUT = str(asset("computation-outputs/C2"))
DEFAULT_OUT = "results/C2"


def run(nsample, src_dir, out):
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261030)
    res = []
    lat = open(out / "glue-lattices.jsonl", "w")
    for z, d, s in pg.TYPES3:
        p = 3
        Fs = pg.F_lattices(p, d, s, src_dir)
        if d == 8:
            from sage.all import CartanMatrix
            E8 = matrix(ZZ, CartanMatrix(["E", 8]))
            Fs = [(3 * E8, 696729600)]
        if d == 6:
            from sage.all import CartanMatrix
            E6 = matrix(ZZ, CartanMatrix(["E", 6]))
            Fs = [((3 * E6.inverse()).change_ring(ZZ), 103680)]
        Zs = pg.Z_lattices(z, s, src_dir)
        if not Fs or not Zs:
            print(f"3-({z},{d})-{s}: no list yet (F {None if Fs is None else len(Fs)}, Z {None if Zs is None else len(Zs)})", flush=True)
            continue
        for iF, (F, aF) in enumerate(Fs):
            DF = pg.Disc(F, p)
            A = DF.norm_classes(2)
            for iZ, (Z, aZ) in enumerate(Zs):
                DZ = pg.Disc(Z, p)
                assert DF.s == DZ.s == s
                B = DZ.norm_classes(2 * (p - 1))
                good = 0
                n = 0
                example = None
                for _ in range(nsample if A and B else 1):
                    T = pg.random_phi(DF, DZ, rng)
                    if T is None:
                        continue
                    n += 1
                    if A and B:
                        arr = (np.array(sorted(A), dtype=np.int64) @ T) % p
                        hit = any(tuple(int(x) for x in row) in B for row in arr)
                    else:
                        hit = False
                    if not hit:
                        good += 1
                        if example is None:
                            example = T
                rec = dict(p=p, z=z, d=d, s=s, iF=iF, iZ=iZ, A=len(A), B=len(B), samples=n, good=good, autF=str(aF), autZ=str(aZ))
                if example is not None:
                    L = pg.build_L(F, Z, DF, DZ, example)
                    q = pari(L)
                    nm2 = int(q.qfminim(2, 0)[0])
                    nm4 = int(q.qfminim(4, 0)[0])
                    assert nm2 == 0, "error in the argument: the lattice has a root"
                    rec["built"] = dict(norm2=nm2, norm4=nm4)
                    lat.write(json.dumps(dict(type=f"3-({z},{d})-{s}", iF=iF, iZ=iZ, gram=[list(map(int, r)) for r in L.rows()])) + "\n")
                    lat.flush()
                print(rec, flush=True)
                res.append(rec)
    lat.close()
    write_json_atomic_indent0(out / "glue.json", res)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("samples", nargs="?", type=int, default=3000, help="random anti-isometries per pair / 組ごとの乱択の反等長写像の数")
    ap.add_argument("--input-dir", default=DEFAULT_INPUT, help="directory of the lists of F and Z / F と Z の一覧のディレクトリ")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    run(a.samples, Path(a.input_dir), Path(a.out))


if __name__ == "__main__":
    main()
