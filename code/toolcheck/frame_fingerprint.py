"""Frame distributions a_k (fingerprints) of extremal even unimodular lattices of rank 32.

Paper: an auxiliary computation, not in Table 11 of the first paper (the paper does not use the frame
distributions). The fingerprints are an invariant for comparing lattices; the checks of the tools (C37)
use them on the lattices of a random walk. The functions are in rank32/frame_fingerprint.py.
Usage:
  python3 frame_fingerprint.py run --lattices FILE [--lattices FILE ...] [--code-grams FILE] [--out FILE] ID ...
    ID is the id of a lattice in a --lattices file (JSON Lines of rows id, gram), code:NAME for the
    lattice NAME of --code-grams (a JSON object NAME -> {gram}), or ID@uS (ID after the random change of
    basis with seed S; a positive control).
Output: --out FILE   (default results/frame-fingerprints.jsonl) JSON Lines, one row per lattice: id, seconds,
                     wall, threads, saturated, min_q1, mincls, hist0, hist1, top (rows already in the file are
                     skipped, so the same command resumes)
Needs: Sage (PARI), gcc with OpenMP. Time: about one minute and 2 GB per lattice (4 threads).

階数 32 の極値偶ユニモジュラー格子の枠の分布 a_k（指紋）。

論文：第 1 の論文の表 11 に無い補助の計算（論文は枠の分布を使わない）。指紋は格子を比べるための不変量で、
道具の検査（C37）がランダムウォークの格子で使う。関数は rank32/frame_fingerprint.py にある。
使い方：
  python3 frame_fingerprint.py run --lattices FILE [--lattices FILE ...] [--code-grams FILE] [--out FILE] ID ...
    ID は --lattices のファイル（id、gram の行の JSON Lines）の格子の id、--code-grams（NAME -> {gram} の
    JSON オブジェクト）の格子 NAME を指す code:NAME、または ID@uS（ID に種 S の乱択基底変換をかけたもの。陽性対照）。
出力：--out FILE   （既定は results/frame-fingerprints.jsonl）JSON Lines。1 行 1 格子：id、seconds、wall、
                   threads、saturated、min_q1、mincls、hist0、hist1、top（ファイルにある行は飛ばすので、同じコマンドで
                   再開できる）
必要なもの：Sage（PARI）、OpenMP つきの gcc。時間：1 格子約 1 分・2 GB（4 スレッド）。
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

from rank32 import frame_fingerprint as fh  # noqa: E402


def run(ids, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {json.loads(l)["id"] for l in open(path)} if path.exists() else set()
    for name in ids:
        if name in done:
            continue
        res = fh.run_one(name)
        rec = dict(id=name, seconds=res["seconds"], wall=res["wall"], threads=res["threads"], saturated=res["saturated"],
                   min_q1=res["min_q1"], mincls=res["mincls"], hist0=res["hist0"], hist1=res["hist1"], top=res["top"])
        d = fh.derive(rec["hist0"], sum(c for _, c in rec["mincls"]))
        with open(path, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(name, res["seconds"], "s  p =", round(d["p"], 7), "fp", d["fp"], flush=True)


def main():
    ap = argparse.ArgumentParser(description="Frame distributions a_k (auxiliary) / 枠の分布 a_k（補助）")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("ids", nargs="+")
    p.add_argument("--lattices", type=Path, action="append", default=[], help="JSON Lines of (id, gram) / (id, gram) の JSON Lines")
    p.add_argument("--code-grams", type=Path, help="JSON object NAME -> {gram}, used as code:NAME / code:NAME として使う")
    p.add_argument("--out", type=Path, default=Path("results/frame-fingerprints.jsonl"), help="output / 出力")
    a = ap.parse_args()
    for f in a.lattices:
        fh.load_lattices(f)
    if a.code_grams:
        fh.load_lattices(a.code_grams, prefix="code:")
    run(a.ids, a.out)


if __name__ == "__main__":
    main()
