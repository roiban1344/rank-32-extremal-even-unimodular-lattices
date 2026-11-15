"""Rootless mass M_2 of the Gaussian Hermitian lattices of rank 16, and the masses per root system.

Paper: computation C11 of the first paper (the rootless mass M_2 over Z[i]); the target value of the
Gaussian enumeration (computation C12). The functions are in rank32/hermitian_mass_gaussian.py.
Commands (all write to --out, default results/C11):
  solve R [CPU_LIMIT]   m(emptyset) for rank R (R = 4, 8, 12, 16) -> solve-rR.json; resumable from
                        cache-rR.json; with CPU_LIMIT (seconds) it stops with status 3 and resumes later
  table R               the mass m(R') of every root system R' -> table-rR.json
  splitcheck            the check m_r(R + E) = m_{r-k}(R) / |U(E)| for E = E_8^g, E_8[i] on the tables
                        of rank 4, 8, 12, 16 -> splitcheck.json
  target                the target of the enumeration: m(emptyset) and the masses by root system and by
                        number of roots, from table-r16.json and solve-r16.json -> enumeration-target.json
  check                 small checks: dynamic programming against brute force, known values, block equivalences
Inputs: none (only the formula and the tables of root systems).
Time: rank 4 and 8 take seconds; rank 12 about 4 minutes; rank 16 about one CPU hour (solve).

Z[i] 上の階数 16 のエルミート格子のルート無しの質量 M_2 と、ルート系ごとの質量。

論文：第 1 の論文の計算 C11（Z[i] 上のルート無しの質量 M_2）。Gauss の格子の列挙（計算 C12）の目標値。
関数は rank32/hermitian_mass_gaussian.py にある。
コマンド（どれも --out に書く。既定は results/C11）：
  solve R [CPU_LIMIT]   階数 R（R = 4, 8, 12, 16）の m(emptyset) -> solve-rR.json。cache-rR.json から再開できる。
                        CPU_LIMIT（秒）を与えると、そこで終了状態 3 で止まり、後で再開できる
  table R               各ルート系 R' の質量 m(R') -> table-rR.json
  splitcheck            階数 4, 8, 12, 16 の表で、E = E_8^g, E_8[i] について m_r(R + E) = m_{r-k}(R) / |U(E)|
                        を確かめる -> splitcheck.json
  target                列挙の目標値：m(emptyset) と、ルート系ごと・ルートの本数ごとの質量。table-r16.json と
                        solve-r16.json から -> enumeration-target.json
  check                 小さい検算：動的計画と全列挙の照合、既知の値、ブロックの同値
入力：無し（式とルート系の表だけから計算する）。
時間：階数 4 と 8 は数秒。階数 12 は約 4 分。階数 16 は solve が約 1 CPU 時間。
"""
import argparse
import ast
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from math import factorial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32 import hermitian_mass_gaussian as hm  # noqa: E402
from rank32.hermitian_mass_gaussian import (A_value, S_mass, U_order, all_elements, ann_of, brute_stats,  # noqa: E402
                                            direct_gram, herm, roots_of, stats_cached)


def mode_solve(r, limit):
    try:
        print(hm.solve(r, cpu_limit=limit))
    except hm.CpuLimit:
        print("stopped at the CPU limit (rerun the same command to resume)", flush=True)
        sys.exit(3)


def mode_table(r):
    M = hm.table(r)
    nz = {k: v for k, v in M.items() if v != 0}
    neg = [k for k, v in nz.items() if v < 0]
    print(f"r={r}: nonzero {len(nz)}, negative {len(neg)}, sum {sum(M.values())}, mu = {S_mass(r)}")
    rows = []
    for k, v in sorted(nz.items(), key=lambda kv: -roots_of(kv[0])):
        rows.append({"type": str(k), "mass": str(v), "inverse": str(1 / v), "roots": roots_of(k)})
        print(k, v, "1/m =", 1 / v)
    hm._write_json(hm.OUT / f"table-r{r}.json", {"r": r, "nonzero": len(nz), "negative": len(neg), "sum": str(sum(M.values())), "mu": str(S_mass(r)), "rows": rows})


def mode_check():
    # 1. dynamic programming against brute force / 動的計画と全列挙の照合
    d11, d13, d21, d23 = ("d", 1, 1), ("d", 1, 3), ("d", 2, 1), ("d", 2, 3)
    H1, P1, P0 = ("H", 1), ("P", 0, 1), ("P", 0, 0)
    tests = [[d11], [d11] * 2, [d11] * 3, [d11] * 4, [d11, d13], [d11, d11, d13], [d13] * 3, [d21], [d23], [d21, d21], [d21, d23], [d21, d11],
             [d23, d11], [d21, d13], [d21, d21, d11], [("d", 3, 1)], [("d", 3, 3), d11], [("d", 3, 1), d21], [("d", 4, 3)], [H1], [H1, d11], [H1, d13],
             [H1, d21], [H1, H1], [P1], [P0], [P1, P0], [P1, P1], [P0, P0], [P1, P1, P0], [P1, d11], [P0, d11], [P1, d13], [P0, d13], [P1, d11, d11],
             [P0, d21], [P1, d21], [H1, P1], [H1, P0], [H1, P1, P0], [P1, P1, d11]]
    bad = 0
    for b in tests:
        dp = stats_cached(tuple(sorted(b)))
        br = brute_stats(b)
        ok = dp == br
        bad += not ok
        print("OK " if ok else "NG ", b, flush=True)
        if not ok:
            print("  dp", dp, "\n  br", br)
    print("mismatches between dynamic programming and brute force:", bad, "/", len(tests))
    # 2. known values: mean number of roots, unimodular T / 既知の値：ルートの本数の平均、ユニモジュラーの T
    print("mean number of roots r=16:", A_value({("A", 1): 1}, 16) / S_mass(16), Fr(16320, 3617))
    print("mean number of roots r=12:", A_value({("A", 1): 1}, 12) / S_mass(12), Fr(65520, 691))
    print("A_8(E_8^g) = mu_4:", A_value({("Eg", 8): 1}, 8) == S_mass(4), " A_12(E_8[i]) = mu_4:", A_value({("E", 8): 1}, 12) == S_mass(4))

    # 3. equivalence of blocks: the distributions of (ann, q) agree / ブロックの同値：元の (ann, q) の分布が一致するか
    def dist(bl):
        A = direct_gram(bl)
        c = Counter()
        for x in all_elements(tuple(bl)):
            c[(ann_of(x), herm(x, A, x).a % 2)] += 1
        return c

    for e in (1, 2):
        print("equivalence of d^2, e =", e, dist([("d", e, 1)] * 2) == dist([("d", e, 3)] * 2))
    for e in (0, 1):
        print("equivalence of P^2, e =", e, dist([("P", e, 1)] * 2) == dist([("P", e, 0)] * 2))


def mode_splitcheck():
    # An independent check: a lattice whose root system contains E_8^g or E_8[i] has it as an orthogonal summand,
    # so m_r(R + E) = m_{r-k}(R) / |U(E)| (with the multiplicity k!).
    # 独立の検査：ルート系に E_8^g、E_8[i] を含む L はそれを直和因子に持つので、m_r(R + E) = m_{r-k}(R) / |U(E)|（多重度 k! つき）
    tabs = {}
    for r in (4, 8, 12, 16):
        tb = json.load(open(hm.OUT / f"table-r{r}.json"))
        tabs[r] = {ast.literal_eval(x["type"]): Fr(x["mass"]) for x in tb["rows"]}
    tabs[0] = {(): Fr(1)}
    UN = {("Eg", 8): (4, 46080), ("E", 8): (8, U_order(("E", 8)))}
    checked = bad = 0
    for r in (8, 12, 16):
        for R, m in tabs[r].items():
            d = dict(R)
            red = dict(d)
            rr = r
            den = 1
            for Y, (k, u) in UN.items():
                mult = red.pop(Y, 0)
                rr -= k * mult
                den *= u ** mult * factorial(mult)
            if rr == r:
                continue
            key = tuple(sorted(red.items()))
            want = tabs[rr].get(key, Fr(0)) / den
            checked += 1
            if want != m:
                bad += 1
                print("NG", r, R, m, want)
        # converse: R built from a nonzero entry of smaller rank is in the table / 逆向き：小さい階数の非零から作った R が表にあること
        for k, (Y, (rk, u)) in enumerate(UN.items()):
            if r - rk not in tabs:
                continue
            for R0, m0 in tabs[r - rk].items():
                d = dict(R0)
                if Y in d:
                    continue
                d[Y] = 1
                key = tuple(sorted(d.items()))
                checked += 1
                # the denominators were checked above; here only existence / 分母は上で検査済み。ここは存在だけ
                if key not in tabs[r]:
                    bad += 1
                    print("missing", r, key)
    res = {"checked": checked, "bad": bad}
    hm._write_json(hm.OUT / "splitcheck.json", res)
    print(res)


def mode_target():
    # The target of the Gaussian enumeration: the exact fraction, and the masses per root system and per number of roots.
    # Gauss の格子の列挙が使う目標の値：厳密な分数と、ルート系ごと・ルートの本数ごとの質量
    tb = json.load(open(hm.OUT / "table-r16.json"))
    sv = json.load(open(hm.OUT / "solve-r16.json"))

    def nm(Y):
        t, n = Y
        return {"A": f"A_{n}[i]", "D": f"D_{n}[i]", "E": f"E_{n}[i]", "G": f"D_{2 * n}^g", "Eg": "E_8^g"}[t]

    rows = []
    strata = defaultdict(Fr)
    for x in tb["rows"]:
        T = ast.literal_eval(x["type"])
        m = Fr(x["mass"])
        strata[x["roots"]] += m
        rows.append({"root_lattice": " + ".join(f"{nm(Y)}^{k}" if k > 1 else nm(Y) for Y, k in T) or "(rootless)", "key": x["type"], "roots": x["roots"], "mass": x["mass"], "mass_float": float(m)})
    me = Fr(sv["m_empty"])
    mu = S_mass(16)
    out = {
        "purpose": "target value for the completion of the enumeration of the Gaussian lattices of rank 16 (computation C12)",
        "target_mass": {"fraction": str(me), "numerator": me.numerator, "denominator": me.denominator, "denominator_factored": "2^21 * 3^3 * 5^3 * 17", "float": float(me), "fraction_of_mu16": float(me / mu)},
        "mu16": {"fraction": str(mu), "float": float(mu)},
        "how_to_use": "m(empty) = sum 1/|U| is the sum over Hermitian classes, complex conjugates counted separately. Summing over classes up to conjugation, a class without an antilinear element counts 2/|U| and a class with one counts 1/|U|. The enumeration is complete when the sum equals target_mass exactly; a larger sum means an error in the identification or in |U|.",
        "roots_meaning": "roots is the number of vectors of norm 2 of the Z-lattice (4 times the number of Z[i]-lines)",
        "strata_by_root_count": {str(k): {"mass": str(v), "float": float(v), "fraction_of_mu16": float(v / mu)} for k, v in sorted(strata.items())},
        "by_root_lattice": rows,
        "reproduce": "rootless_hermitian_mass_gaussian.py solve 16; table 16; target",
    }
    hm._write_json(hm.OUT / "enumeration-target.json", out)
    print(out["target_mass"], len(rows), len(out["strata_by_root_count"]))


def main():
    ap = argparse.ArgumentParser(description="Rootless mass of Gaussian Hermitian lattices (C11) / Z[i] 上のルート無しの質量（C11）")
    ap.add_argument("--out", type=Path, default=Path("results/C11"), help="output directory / 出力先")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("solve")
    p.add_argument("r", type=int)
    p.add_argument("cpu_limit", type=float, nargs="?")
    p = sub.add_parser("table")
    p.add_argument("r", type=int)
    sub.add_parser("splitcheck")
    sub.add_parser("target")
    sub.add_parser("check")
    a = ap.parse_args()
    hm.OUT = a.out
    hm.OUT.mkdir(parents=True, exist_ok=True)
    if a.cmd == "solve":
        mode_solve(a.r, a.cpu_limit)
    elif a.cmd == "table":
        mode_table(a.r)
    elif a.cmd == "splitcheck":
        mode_splitcheck()
    elif a.cmd == "target":
        mode_target()
    else:
        mode_check()


if __name__ == "__main__":
    main()
