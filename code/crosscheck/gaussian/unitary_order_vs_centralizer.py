"""Check of the unitary group orders |U| of the backtracking search against centralizers in Aut(L).

Paper: computation C12 of the first paper (a check of |U| used for the masses of the Gaussian classes).
For each lattice of a sample with a known generating set of Aut(L) (from the bliss automorphism tool),
|U| = |C_{Aut(L)}(g)| is computed by GAP, through the reduction of the generators mod 3 (faithful on the
finite group Aut(L)), and compared with |U| from rank32.gaussian_keys.unitary_order.
Usage:
  python3 unitary_order_vs_centralizer.py [--sample FILE] [--orders FILE] [--out FILE]
Input (by default in computation-outputs/C12 of the release assets):
        --sample FILE   JSON Lines with id, gram, W (lattices of the walk)
        --orders FILE   JSON Lines with id, order, generators (Aut(L) of the same lattices)
Output: --out FILE      (default results/C12/unitary-order-vs-centralizer.jsonl) JSON Lines: id, aut,
                        centralizer, U, agree (resumable: done ids are skipped)
Needs: Sage with GAP, gcc. Time: a few seconds per lattice.

後戻りの探索によるユニタリ群の位数 |U| を、Aut(L) の中の中心化群と照合する。

論文：第 1 の論文の計算 C12（Gauss の類の質量に使う |U| の検算）。
Aut(L) の生成元が分かっている標本（bliss の自己同型の道具から）の各格子について、
|U| = |C_{Aut(L)}(g)| を GAP で求める（生成元を mod 3 に落とす。有限群 Aut(L) の上で忠実）。
それを rank32.gaussian_keys.unitary_order の |U| と比べる。
使い方：
  python3 unitary_order_vs_centralizer.py [--sample FILE] [--orders FILE] [--out FILE]
入力（既定は Release のアセットの computation-outputs/C12）：
      --sample FILE   id、gram、W の JSON Lines（ランダムウォークの格子）
      --orders FILE   id、order、generators の JSON Lines（同じ格子の Aut(L)）
出力：--out FILE      （既定は results/C12/unitary-order-vs-centralizer.jsonl）JSON Lines：id、aut、
                      centralizer、U、agree（再開できる。済んだ id は飛ばす）
必要なもの：Sage（GAP つき）、gcc。時間：1 格子数秒。
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32.gaussian_keys import read_jsonl, unitary_order  # noqa: E402
from rank32.io import asset  # noqa: E402


def crosscheck(sample, orders, out):
    """Compare |U| of the search with |C_Aut(g)| computed by GAP from the generators of Aut(L).

    探索の |U| と、Aut(L) の生成元から GAP で求めた |C_Aut(g)| を照合する。"""
    from sage.all import GF, ZZ, matrix, pari
    from sage.libs.gap.libgap import libgap
    pari.allocatemem(1 << 28, 1 << 30)
    rows = {r["id"]: r for r in read_jsonl(sample)}
    res = out
    done = {r["id"] for r in read_jsonl(res)}
    for o in read_jsonl(orders):
        if o["id"] in done:
            continue
        r = rows[o["id"]]
        G, W = matrix(ZZ, r["gram"]), matrix(ZZ, r["W"])
        gens = [matrix(ZZ, g) for g in o["generators"]]
        if not all(g * G * g.T == G for g in gens):
            gens = [g.T for g in gens]
            assert all(g * G * g.T == G for g in gens)
        K = libgap.Group([libgap(matrix(GF(3), g)) for g in gens])
        C = int(libgap.Size(libgap.Centralizer(K, libgap(matrix(GF(3), W))))) if len(gens) else None
        a = unitary_order(r["gram"], r["W"])
        rec = dict(id=o["id"], aut=str(o["order"]), centralizer=C, U=a["U"], agree=C == a["U"])
        with res.open("a") as f:
            f.write(json.dumps(rec) + "\n")
        print(rec, flush=True)


def main():
    ap = argparse.ArgumentParser(description="Check |U| against centralizers (C12) / |U| と中心化群の照合（C12）")
    ap.add_argument("--sample", type=Path, default=asset("computation-outputs/C12/gaussian-sample.jsonl"), help="lattices / 格子")
    ap.add_argument("--orders", type=Path, default=asset("computation-outputs/C12/gaussian-sample-orders.jsonl"), help="Aut(L) of the lattices / 格子の Aut(L)")
    ap.add_argument("--out", type=Path, default=Path("results/C12/unitary-order-vs-centralizer.jsonl"), help="output / 出力")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    crosscheck(a.sample, a.orders, a.out)


if __name__ == "__main__":
    main()
