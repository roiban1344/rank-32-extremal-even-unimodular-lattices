"""Masses of the rootless Hermitian unimodular lattices over Z[zeta_p], p <= 31, rank up to 32/(p - 1).

Paper: computation C2 of the first paper (Table 11); the types whose cyclic part would need a rootless
unimodular lattice of mass 0 are excluded. The formula is in rank32.rootless_unimodular_mass.
Output: --out DIR, the file noroot.json ({p: [{z, mass, mass_float}, ...]})
Needs:  SageMath. Time: a few seconds.

Z[zeta_p] 上のルートを持たないエルミート・ユニモジュラー格子の質量（p <= 31、階数は 32/(p - 1) まで）。

論文：第 1 の論文の計算 C2（表 11）。巡回部分に質量 0 のルート無しのユニモジュラー格子が要る型を除く。
式は rank32.rootless_unimodular_mass にある。
出力：--out DIR の noroot.json（{p: [{z, mass, mass_float}, ...]}）
必要なもの：SageMath。時間：数秒。
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32.io import write_json_atomic_indent0  # noqa: E402
from rank32.rootless_unimodular_mass import noroot_masses  # noqa: E402

DEFAULT_OUT = "results/C2"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    rows = {}
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        zmax = 32 // (p - 1)
        ms = noroot_masses(p, zmax)
        rows[str(p)] = [dict(z=z, mass=str(m), mass_float=float(m)) for z, m in enumerate(ms)]
        print(f"p = {p}")
        for z, m in enumerate(ms):
            assert m >= 0, (p, z, m)
            extra = f"  1/m = {float(1 / m):.6g}" if m else ""
            print(f"  z = {z}: {float(m):.6g} {'= ' + str(m) if m.denominator() < 10**15 else ''}{extra}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    write_json_atomic_indent0(out / "noroot.json", rows)


if __name__ == "__main__":
    main()
