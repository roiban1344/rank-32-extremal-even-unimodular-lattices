"""All double cosets of the gluing of one pair of components (F, K) of an involution, closed by the mass.

Paper: computation C4 of the first paper (Table 11). The method is in rank32.involution_gluing.
Usage:  python3 involution_gluing_classes.py PAIR [LIMIT_SECONDS]
  PAIR: 14-18a ... 14-18e, 15-17, 16-16s14, 16-16d8, 16-16-s<s>d<delta>..., 12-20-c<i>
Input (by default in computation-outputs/C4 and C7 of the release assets):
        --groups FILE (groups.json, from discriminant_action.py), --groups-12-20 FILE, --scale FILE
        (|O(q)|, involution-scale.json), --parts DIR (the lists of components, output of C7)
Output: --out DIR (default results/C4)   <PAIR>.json (the list of classes: phi, |Stab|, extremal) when the
        mass is complete, <PAIR>.progress.json, <PAIR>.runs.jsonl (seed of every start), <PAIR>.ckpt.json
        (checkpoint; the same command resumes).
Time: seconds to minutes for most pairs; 16-16d8 took several days.

対合の成分の組 (F, K) の 1 つについて、貼り合わせの両側剰余類の全部を、質量で閉じて求める。

論文：第 1 の論文の計算 C4（表 11）。方法は rank32.involution_gluing にある。
使い方：python3 involution_gluing_classes.py PAIR [LIMIT_SECONDS]
  PAIR：14-18a ... 14-18e、15-17、16-16s14、16-16d8、16-16-s<s>d<delta>...、12-20-c<i>
入力（既定は Release のアセットの computation-outputs/C4 と C7）：
      --groups FILE（groups.json。discriminant_action.py の出力）、--groups-12-20 FILE、--scale FILE
      （|O(q)|、involution-scale.json）、--parts DIR（成分の一覧。C7 の出力）
出力：--out DIR（既定は results/C4）   質量が満ちたら <PAIR>.json（類の一覧：phi、|Stab|、極値か）、
      <PAIR>.progress.json、<PAIR>.runs.jsonl（起動ごとの乱数の種）、<PAIR>.ckpt.json（チェックポイント。同じコマンドで再開する）。
時間：多くの組は数秒から数分。16-16d8 は数日かかった。
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32 import involution_gluing as gl  # noqa: E402

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pair", help="name of the pair / 組の名前")
    ap.add_argument("limit", nargs="?", type=float, default=600, help="time limit in seconds / 制限秒")
    ap.add_argument("--groups", type=Path, default=gl.GROUPS, help="groups.json")
    ap.add_argument("--groups-12-20", type=Path, default=gl.GROUPS_12_20, help="table of the pairs (12, 20) / (12, 20) の組の表")
    ap.add_argument("--scale", type=Path, default=gl.SCALE_FILE, help="orders |O(q)| / |O(q)| の表")
    ap.add_argument("--parts", type=Path, default=gl.PARTS_DIR, help="directory of min4-twoel-*.json / min4-twoel-*.json のディレクトリ")
    ap.add_argument("--out", type=Path, default=gl.OUT, help="output directory / 出力先")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    gl.load_tables(a.groups, a.groups_12_20, a.scale, a.parts)
    gl.search(a.pair, a.limit, a.out)
