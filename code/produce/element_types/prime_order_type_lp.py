"""The joint linear program on the theta series of F and Z for each type p-(z, d)-s in rank 32 (p = 3, 5, 7).

Paper: computation C2 of the first paper (Table 11), the exclusion of types by linear programming.
The program is described in rank32.joint_theta_lp; the theta series of the blocks are in
rank32.theta_lp_basis.
Usage:  python3 prime_order_type_lp.py [p ...]   (default: 3 5 7)
Output: --out DIR (default results/C2), the file lp-p<p>-<p>-....json: one row per type. The field
        verdict is "excluded" (the program is infeasible) or "remains", and for the types that remain
        the ranges of c_F[2], c_Fg[1], c_Fg[2], c_Z[2], c_Zg[p-1] are given (an end without bound is
        "unbounded"). If the basis of theta series for F or Z is not available, verdict is
        "undecided F:<a> Z:<b>", where <a> and <b> are "ok", "no genus",
        "undecided (rank r < dimension m)" (the blocks do not span the space) or "check failed". The constants are those of rank32.joint_theta_lp.
Needs:  SageMath (with PPL, exact rational arithmetic). The address space is limited to 4 GB (RLIMIT_AS).
        Time: about 1 minute.

階数 32 の型 p-(z, d)-s ごとに、F と Z のテータ級数の同時線形計画を解く（p = 3, 5, 7）。

論文：第 1 の論文の計算 C2（表 11）。線形計画による型の除外。線形計画は rank32.joint_theta_lp に、
ブロックのテータ級数は rank32.theta_lp_basis に書いてある。
使い方：python3 prime_order_type_lp.py [p ...]（既定：3 5 7）
出力：--out DIR（既定は results/C2）の lp-p<p>-<p>-....json。型ごとに 1 行。欄 verdict は "excluded"
      （実行不可能）か "remains" で、残る型には c_F[2], c_Fg[1], c_Fg[2], c_Z[2], c_Zg[p-1] の範囲を付ける
      （有界でない端は "unbounded"）。F か Z のテータ級数の基底が得られなければ、verdict は
      "undecided F:<a> Z:<b>" で、<a> と <b> は "ok"、"no genus"、
      "undecided (rank r < dimension m)"（ブロックが空間を張らない）、"check failed" のどれか。定数は rank32.joint_theta_lp のもの。
必要なもの：SageMath（PPL を含む。有理数で厳密に解く）。アドレス空間を 4 GB に制限する（RLIMIT_AS）。
      時間：約 1 分。
"""
import argparse
import resource
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32 import joint_theta_lp as lp  # noqa: E402
from rank32.io import write_json_atomic_indent0  # noqa: E402

tb = lp.tb
resource.setrlimit(resource.RLIMIT_AS, (4 << 30, 4 << 30))

DEFAULT_OUT = "results/C2"
EXCLUDED, REMAINS, UNDECIDED, NO_GENUS, CHECK_FAILED = lp.EXCLUDED, lp.REMAINS, lp.UNDECIDED, lp.NO_GENUS, lp.CHECK_FAILED


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("primes", nargs="*", type=int, help="the primes (default 3 5 7) / 素数（既定 3 5 7）")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ps = a.primes or [3, 5, 7]
    results = []
    for p in ps:
        th = {}
        for name, G, Gd, n, s in tb.blocks(p):
            th[name] = (tb.theta_q(G, 1), tb.theta_q(Gd, 1), n, s)
        cache = {}
        def sp(n, s):
            if (n, s) not in cache:
                pairs, rank, dim = tb.spaces_for(p, n, s, th)
                if pairs is None:
                    cache[(n, s)] = (NO_GENUS, None)
                elif rank < dim:
                    cache[(n, s)] = (f"{UNDECIDED} (rank {rank} < dimension {dim})", None)
                elif not tb.consistent(pairs):
                    cache[(n, s)] = (CHECK_FAILED, None)
                else:
                    cache[(n, s)] = ("ok", lp.basis(pairs))
            return cache[(n, s)]
        for z, d, s in lp.types_for(p):
            n2 = z * (p - 1)
            stZ, bZ = sp(n2, s)
            stF, bF = sp(d, s) if d else ("ok", [])
            if stZ != "ok" or stF != "ok":
                rec = dict(p=p, z=z, d=d, s=s, verdict=f"{UNDECIDED} F:{stF} Z:{stZ}")
            else:
                N = lp.class_counts(p, d, s) if d else {}
                r = lp.solve(p, d, n2, s, bF, bZ, N)
                rec = dict(p=p, z=z, d=d, s=s, verdict=EXCLUDED if r is None else REMAINS, ranges=r)
            print(rec, flush=True)
            results.append(rec)
    write_json_atomic_indent0(out / f"lp-p{'-'.join(map(str, ps))}.json", results)


if __name__ == "__main__":
    main()
