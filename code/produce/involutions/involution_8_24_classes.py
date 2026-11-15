"""The type (8, 24) of involutions: from the orbits W of a Niemeier lattice, the lattices M and their invariants, one line per W.

Paper: computation C6 of the first paper (Table 11). The method is in rank32.involution_8_24_classes.
Usage:  python3 involution_8_24_classes.py P [PROCESSES]
  P: Leech, A1 or A2 (the Niemeier lattice: entries 179, 177, 178 of rank 24 of the list), or an entry number
Input:  --unimodular FILE   the unimodular lattices of rank <= 24 without norm 1 (key "24"; default
                            computation-outputs/C6/unimodular-min2.json in the release assets)
        --orbits ORBITS     <ORBITS>/<P>/W.jsonl, the orbits W (output of niemeier_good_subspace_orbits.py,
                            added in a later version); default computation-outputs/C6 in the release assets
                            (orbits computed here are in results/C6: give --orbits results/C6)
Output: --out OUT (default results/C6): <OUT>/<P>/classes-<shard>.jsonl, one line per W of the shard, in the order of W
        (resumable: a shard continues after its last line), and progress-classes.json.
        Each line: cert (complete invariant), index, method, orbit, |Aut M|, census, np, i, stab, gram.
Time: Leech about 11 minutes; A1 about 9 hours with 4 processes (with the orbits W).

対合の (8, 24) 型：Niemeier 格子の軌道 W から格子 M と不変量を求める。W ごとに 1 行。

論文：第 1 の論文の計算 C6（表 11）。方法は rank32.involution_8_24_classes にある。
使い方：python3 involution_8_24_classes.py P [PROCESSES]
  P：Leech、A1、A2（Niemeier 格子。一覧の階数 24 の 179、177、178 番）、または一覧の番号
入力：--unimodular FILE   ノルム 1 を持たない階数 24 以下のユニモジュラー格子（キー "24"。既定は
                          Release のアセットの computation-outputs/C6/unimodular-min2.json）
      --orbits ORBITS     <ORBITS>/<P>/W.jsonl。軌道 W（niemeier_good_subspace_orbits.py の出力。この
                          プログラムは後の版で加わる）。既定は Release のアセットの
                          computation-outputs/C6（自分で求めた軌道は results/C6 にある。--orbits results/C6 と渡す）
出力：--out OUT（既定は results/C6）。<OUT>/<P>/classes-<shard>.jsonl。シャードの W ごとに 1 行、W の順（再開できる：シャードは最後の行の
      続きから動く）と progress-classes.json。各行：cert（完全不変量）、index、method、orbit、|Aut M|、
      census、np、i、stab、gram。
時間：Leech は約 11 分。A1 は 4 プロセスで約 9 時間（軌道 W と合わせて）。
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix  # noqa: E402

from rank32 import involution_8_24_classes as cl  # noqa: E402


def worker(args):
    tag, shard, nshards, unimod, orbits, outdir = args
    lib = cl.load_canon()
    idx = cl.NAMES.get(tag)
    idx = int(tag) if idx is None else idx
    G = matrix(ZZ, json.load(open(unimod))["24"][idx]["gram"])
    out = Path(outdir) / tag
    out.mkdir(parents=True, exist_ok=True)
    recs = [json.loads(l) for l in open(Path(orbits) / tag / "W.jsonl")]
    recs = [r for r in recs if "W" in r]
    fn = out / f"classes-{shard}.jsonl"
    done = 0
    if fn.exists():
        done = sum(1 for _ in open(fn))
    mine = [i for i in range(len(recs)) if i % nshards == shard][done:]
    t0 = time.time()
    with open(fn, "a") as f:
        for n, i in enumerate(mine):
            r = recs[i]
            a, Mg = cl.analyse(lib, cl.M_from_W(G, r["W"]))
            a.update(i=i, stab=r["stab"], gram=[[int(x) for x in row] for row in Mg.rows()])
            f.write(json.dumps(a) + "\n")
            f.flush()
            if n % 50 == 0 and shard == 0:
                tmp = str(out / "progress-classes.json") + ".tmp"
                with open(tmp, "w") as g:
                    json.dump(dict(job=f"type-8-24-classes-{tag}", stage="isometry classes", state="running", done=done + n, total=len(mine) + done, unit="W (shard 0)",
                                   updated=time.strftime("%Y-%m-%dT%H:%M:%S%z"), elapsed_s=round(time.time() - t0, 1), shards=nshards), g)
                os.replace(tmp, out / "progress-classes.json")
    return shard


def main(tag, k, unimod, orbits, outdir):
    cl.load_canon()
    with mp.Pool(k) as pool:
        for s in pool.imap_unordered(worker, [(tag, i, k, str(unimod), str(orbits), str(outdir)) for i in range(k)]):
            print("shard", s, "done", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("tag", help="Leech, A1, A2 or an entry number / Leech、A1、A2 または番号")
    ap.add_argument("processes", nargs="?", type=int, default=1, help="number of processes (= shards) / プロセス数（= シャード数）")
    ap.add_argument("--unimodular", type=Path, default=cl.UNIMOD, help="list of unimodular lattices / ユニモジュラー格子の一覧")
    ap.add_argument("--orbits", type=Path, default=cl.ORBITS, help="directory with <P>/W.jsonl / <P>/W.jsonl のあるディレクトリ")
    ap.add_argument("--out", type=Path, default=cl.OUTDIR, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    main(a.tag, a.processes, a.unimodular, a.orbits, a.out)
