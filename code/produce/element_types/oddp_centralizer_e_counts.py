"""The numbers e(h) and |C(h)| for elements of order 5, 7, 17, 31, and E = sum e(h)/|C(h)| per type.

Paper: computation C9 of the first paper (Table 11). The method is in rank32.oddp_e_counts.
Usage (python3 oddp_centralizer_e_counts.py MODE ...):
  check ID...                       direct count against sum |Aut L| e(h)/|C(h)| for known lattices (check.jsonl)
  eknown ID@P...                    (|C(h)|, e(h), fixed points) per class of order P (e-known.jsonl)
  pairs TYPE [all]                  TYPE = 5-7-4-3 or 5-6-8-6 (E-TYPE.json); `all` builds the group for every pair
  generic P Z S D ZLIST TAG         a type whose F is one class (E-TAG.json); ZLIST is the list of Z
  phis P Z S D ZLIST TAG            the same computation, writing phi and L of every pair (phis-TAG.jsonl.gz, phis-E-TAG.json)
  summary                           E per type and per p (E-summary.json)
Paths: --out DIR (outputs, default results/C9; summary reads the outputs of the other modes from there),
--known FILE (the 26 known lattices, output of C32), --pairs-dir DIR (lists of pairs and known-types.jsonl),
--genus-dir DIR (genus lists, output of C2), --example FILE (example-p5-16-4.json).
Needs: SageMath and the GAP executable `gap`. Time: seconds to minutes per command (check: up to 8 minutes per lattice).

位数 5, 7, 17, 31 の元の e(h) と |C(h)|、型ごとの E = sum e(h)/|C(h)|。

論文：第 1 の論文の計算 C9（表 11）。方法は rank32.oddp_e_counts にある。
使い方（python3 oddp_centralizer_e_counts.py MODE ...）：
  check ID...                       既知の格子で、直接の数え上げと sum |Aut L| e(h)/|C(h)| を比べる（check.jsonl）
  eknown ID@P...                    位数 P の類ごとの (|C(h)|, e(h), 不動点の数)（e-known.jsonl）
  pairs TYPE [all]                  TYPE = 5-7-4-3 か 5-6-8-6（E-TYPE.json）。`all` なら全部の対で群を作る
  generic P Z S D ZLIST TAG         F が 1 類の型（E-TAG.json）。ZLIST は Z の一覧
  phis P Z S D ZLIST TAG            同じ計算で、対ごとの phi と L を書く（phis-TAG.jsonl.gz、phis-E-TAG.json）
  summary                           型ごと・p ごとの E（E-summary.json）
パス：--out DIR（出力。既定は results/C9。summary は、ほかのモードの出力をそこから読む）、
--known FILE（既知の 26 格子。C32 の出力）、--pairs-dir DIR（対の一覧と known-types.jsonl）、
--genus-dir DIR（種の一覧。C2 の出力）、--example FILE（example-p5-16-4.json）。
必要なもの：SageMath と GAP の実行ファイル `gap`。時間：コマンドごとに数秒から数分（check は格子 1 つで最大 8 分）。
"""
import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32 import oddp_e_counts as k  # noqa: E402

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mode", choices=["check", "eknown", "pairs", "generic", "phis", "summary"], help="mode / モード")
    ap.add_argument("args", nargs="*", help="arguments of the mode / モードの引数")
    ap.add_argument("--out", type=Path, default=k.OUT, help="output directory (default results/C9) / 出力先（既定は results/C9）")
    ap.add_argument("--known", type=Path, default=k.KNOWN, help="the 26 known lattices / 既知の 26 格子")
    ap.add_argument("--pairs-dir", type=Path, default=k.PAIRS_DIR, help="lists of pairs / 対の一覧")
    ap.add_argument("--genus-dir", type=Path, default=k.GENUS_DIR, help="genus lists / 種の一覧")
    ap.add_argument("--example", type=Path, default=k.EXAMPLE_P5_16_4, help="example-p5-16-4.json")
    a = ap.parse_args()
    k.OUT, k.KNOWN, k.PAIRS_DIR, k.GENUS_DIR, k.EXAMPLE_P5_16_4 = a.out, a.known, a.pairs_dir, a.genus_dir, a.example
    k.OUT.mkdir(parents=True, exist_ok=True)
    x = a.args
    if a.mode == "check":
        k.check(x)
    elif a.mode == "summary":
        k.summary()
    elif a.mode == "eknown":
        for spec in x:
            i_, p_ = spec.rsplit("@", 1)
            k.eknown(i_, int(p_))
    elif a.mode == "pairs":
        k.group_pairs(x[0], nverify=10 ** 9 if x[1:] == ["all"] else 2)
    elif a.mode == "phis":
        with gzip.open(k.OUT / f"phis-{x[5]}.jsonl.gz", "wt") as fh_:
            k.generic(int(x[0]), int(x[1]), int(x[2]), int(x[3]), x[4], x[5], phi_out=fh_)
    elif a.mode == "generic":
        k.generic(int(x[0]), int(x[1]), int(x[2]), int(x[3]), x[4], x[5])
