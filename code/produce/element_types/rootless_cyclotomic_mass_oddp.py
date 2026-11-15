"""Rootless masses of the Hermitian genera over Z[zeta_p] for the types of order 5, 7 and 17, and the short classes of F^#.

Paper: computation C9 of the first paper (Table 11). The formulas are in rank32.cyclotomic_rootless_mass.
Modes:
  values   for (p, z, s) = (5,4,4), (5,5,5), (5,6,4), (5,6,6), (5,7,3), (5,8,0), (7,4,4), (17,2,0):
           the rootless mass, mass_H, their ratio, and m(A^a E^b) for every root system; checks that the
           table sums to mass_H. Output: values.json
  fcount   for the lattices F of minimum >= 4 of the genera (p, d, s) = (5,4,3), (5,8,4), (5,8,6),
           (5,12,5), (7,8,4) and one F of (5,16,4): the number of vectors of norm 2/p of F^# and the
           number of their classes. Output: fcount.json
  test     the tables for small genera (p = 5 (2,0), (4,0) with E^2 = 1/720000, s = z for several p, ...)
  std      mass_H against the standard mass of [Jür15, Satz 3.3.15]
  level    the classes of A_4 / 5 A_4^# and the values of Q on them
Inputs (fcount): --genus-dir DIR with genus-p<p>-n<d>-s<s>.json (output of C2), --example FILE (one F of (16, 4)_5)
Output: --out DIR
Needs:  SageMath. Time: a few seconds for each mode.

位数 5, 7, 17 の型について、Z[zeta_p] 上のエルミートの種のルート無しの質量と、F^# の短い類を求める。

論文：第 1 の論文の計算 C9（表 11）。式は rank32.cyclotomic_rootless_mass にある。
モード：
  values   (p, z, s) = (5,4,4), (5,5,5), (5,6,4), (5,6,6), (5,7,3), (5,8,0), (7,4,4), (17,2,0) について、
           ルート無しの質量、mass_H、その比、全てのルート系の m(A^a E^b)。表の和が mass_H になることを
           確かめる。出力：values.json
  fcount   種 (p, d, s) = (5,4,3), (5,8,4), (5,8,6), (5,12,5), (7,8,4) の最小 >= 4 の格子 F と、(5,16,4) の
           F の 1 つについて、F^# のノルム 2/p のベクトルの数と、その類の数。出力：fcount.json
  test     小さい種の表（p = 5 の (2,0)、E^2 = 1/720000 となる (4,0)、いくつかの p の s = z など）
  std      mass_H を [Jür15, Satz 3.3.15] の標準質量と照合する
  level    A_4 / 5 A_4^# の類と、その上の Q の値
入力（fcount）：--genus-dir DIR の genus-p<p>-n<d>-s<s>.json（C2 の出力）、--example FILE（(16, 4)_5 の F の 1 つ）
出力：--out DIR
必要なもの：SageMath。時間：どのモードも数秒。
"""
import argparse
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32 import cyclotomic_rootless_mass as cm  # noqa: E402
from rank32.io import asset, write_json_atomic_indent1  # noqa: E402

DEFAULT_OUT = "results/C9"
DEFAULT_GENUS_DIR = str(asset("computation-outputs/C2"))
DEFAULT_EXAMPLE = str(asset("computation-outputs/C17/example-p5-16-4.json"))


def test():
    table, rootless = cm.table, cm.rootless
    print("p=5 (2,0) full table", {k: str(v) for k, v in table(5, 2, 0).items() if v}, "m(empty)", rootless(5, 2, 0))
    print("p=5 (4,0) full table", {k: str(v) for k, v in table(5, 4, 0).items() if v})
    print("  expected E^2 = 1/720000:", Fr(1, 720000))
    for p, zs in [(5, range(1, 9)), (7, range(1, 6)), (11, range(1, 3)), (13, range(1, 3)), (17, range(1, 3))]:
        print(p, "s = z:", [(z, str(rootless(p, z, z))) for z in zs])
    print("p=5 (6,0)", str(rootless(5, 6, 0)), {k: str(v) for k, v in table(5, 6, 0).items() if v})
    print("p=7 (4,0)", str(rootless(7, 4, 0)), {k: str(v) for k, v in table(7, 4, 0).items() if v})
    print("p=13 (2,0)", str(rootless(13, 2, 0)), {k: str(v) for k, v in table(13, 2, 0).items() if v})


def values(out):
    rows = []
    for p, z, s in [(5, 4, 4), (5, 5, 5), (5, 6, 4), (5, 6, 6), (5, 7, 3), (5, 8, 0), (7, 4, 4), (17, 2, 0)]:
        M = cm.table(p, z, s)
        r = cm.rootless(p, z, s)
        assert M[(0, 0)] == r
        mh = cm.mass_H(p, z, s)
        assert sum(M.values()) == mh, (p, z, s)
        rows.append({"p": p, "z": z, "s": s, "rootless": str(r), "float": float(r), "mass_H": str(mh), "fraction": float(r / mh), "negative": sum(1 for v in M.values() if v < 0),
                     "table": {f"A^{k[0]} E^{k[1]}": str(v) for k, v in M.items() if v}})
        print(p, z, s, r, float(r), "ratio", float(r / mh), "negative", rows[-1]["negative"], rows[-1]["table"], flush=True)
    out.mkdir(parents=True, exist_ok=True)
    write_json_atomic_indent1(out / "values.json", rows)


def std():
    # [Jür15, Satz 3.3.15]: std_N = 2^(1 - N[F:Q]) prod |L(1 - j, chi^j)| * (1 / 2^(-t)); t = 1 (p ramifies)
    # [Jür15, Satz 3.3.15]：std_N = 2^(1 - N[F:Q]) prod |L(1 - j, chi^j)| * (1 / 2^(-t))。t = 1（p で分岐）
    bad = 0
    for p in (3, 5, 7, 11, 13, 17):
        d = (p - 1) // 2
        for z in range(1, 9 if p <= 7 else 4):
            L = Fr(1)
            for j in range(1, z + 1):
                L *= cm.F(cm.ptmass.l_value(p, j))
            std = Fr(1, 2 ** (z * d - 1)) * L * (Fr(1, 2) if z % 2 else 1)
            s = z % 2
            n = z * (p - 1)  # noqa: F841
            m = cm.mass_H(p, z, s)
            ok = (m == std) or (m == 0)
            bad += not ok
            print(p, z, s, "std", std, "mass_H", m, "OK" if m == std else ("no genus" if m == 0 else "NG"), flush=True)
    print("disagreements", bad)


def level():
    # the F of 5-(7,4)-3 and 5-(6,8)-6: the number of classes with q = 2/5
    # 5-(7,4)-3 と 5-(6,8)-6 の F：q = 2/5 の類の数
    for n, s, name in [(3, 3, "D(sqrt5 A_4^#)")]:
        pass
    import itertools

    # the Gram matrix of A_4, D = A_4 / 5 A_4^#, Q(a) = (a, a)/2 mod 5 (the norms of F^# are (a, a)/5)
    # A_4 の Gram、D = A_4 / 5 A_4^#、Q(a) = (a, a)/2 mod 5（F^# のノルムは (a, a)/5）
    c = [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]]
    seen = {}
    for co in itertools.product(range(5), repeat=4):
        v = sum(co[i] * c[i][j] * co[j] for i in range(4) for j in range(4))
        # the class of a mod 5 A_4^#: in coordinates of A_4^# it is the inverse of c; a ~ a + 5 w (w in A_4^#) <=> c a mod 5 agree
        # 類は a mod 5 A_4^#：A_4^# 座標では c の逆。a ~ a + 5 w（w in A_4^#）<=> c a mod 5 が同じ
        key = tuple(sum(c[i][j] * co[j] for j in range(4)) % 5 for i in range(4))
        seen.setdefault(key, set()).add(v // 2 % 5 if v % 2 == 0 else None)
    print("number of classes", len(seen), "Q takes 1 value on each class", all(len(x) == 1 for x in seen.values()))
    print("classes with Q = 1", sum(1 for x in seen.values() if 1 in x), "(the classes of roots are 20)")


def fcount(out, genus_dir, example):
    # stage 1: count the short classes of F^# (the inputs are the genus lists of C2)
    # ステージ 1：F^# の短い類を数える（入力は C2 の種の一覧）
    from sage.all import ZZ, matrix, pari

    rows = []
    for (p, n, s) in [(5, 4, 3), (5, 8, 4), (5, 8, 6), (5, 12, 5), (7, 8, 4)]:
        g = json.load(open(genus_dir / f"genus-p{p}-n{n}-s{s}.json"))
        for c in g["classes"]:
            if c["min"] < 4:
                continue
            G = matrix(ZZ, c["gram"])
            Dm = (p * G.inverse()).change_ring(ZZ)
            cnt = int(pari(Dm).qfminim(2, 0, 0)[0])
            classes = set()
            if cnt:
                Vm = matrix(ZZ, pari(Dm).qfminim(2, None, 0)[2].sage())
                for col in Vm.columns():
                    for sg in (1, -1):
                        v = G.inverse() * (sg * col)
                        classes.add(tuple(str(x - x.floor()) for x in v))
            rows.append({"p": p, "d": n, "s": s, "aut": c["aut"], "min": c["min"], "norm_2_over_p_vectors": cnt, "classes": len(classes)})
            print(rows[-1], flush=True)
    # the F of (16, 4): one class given by its Gram matrix, not a list / (16, 4) の F：一覧ではなく 1 類の Gram
    ex = json.load(open(example))["hurwitz:s1-21"]
    G = matrix(ZZ, ex["gram"])
    Dm = (5 * G.inverse()).change_ring(ZZ)
    rows.append({"p": 5, "d": 16, "s": 4, "aut": ex["aut"], "min": 4, "norm_2_over_p_vectors": int(pari(Dm).qfminim(2, 0, 0)[0]), "classes": 0})
    print(rows[-1], flush=True)
    out.mkdir(parents=True, exist_ok=True)
    write_json_atomic_indent1(out / "fcount.json", rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mode", choices=["test", "values", "std", "level", "fcount"], help="the mode / モード")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    ap.add_argument("--genus-dir", default=DEFAULT_GENUS_DIR, help="directory of the genus lists (fcount) / 種の一覧のディレクトリ（fcount）")
    ap.add_argument("--example", default=DEFAULT_EXAMPLE, help="the F of (16, 4)_5 (fcount) / (16, 4)_5 の F（fcount）")
    a = ap.parse_args()
    if a.mode == "test":
        test()
    elif a.mode == "values":
        values(Path(a.out))
    elif a.mode == "std":
        std()
    elif a.mode == "level":
        level()
    else:
        fcount(Path(a.out), Path(a.genus_dir), Path(a.example))


if __name__ == "__main__":
    main()
