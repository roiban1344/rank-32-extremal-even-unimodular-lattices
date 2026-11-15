"""Census of the pi-neighbour walk on Gaussian lattices of rank 16: invariants at every step.

Paper: computation C12 of the first paper. The walk continues the two earlier chains of the uniform
pi-neighbour walk from their final states, and records at each step the invariants and the canonical
form of rank32.gaussian_keys.analyse. The final states of these census walks are the starting points
of the two chains of the enumeration (code/produce/gaussian/gaussian_pi_neighbour_chain.py, added in a
later version).
Lattices with roots are kept with their Gram matrix every 30 steps; rootless lattices without a
canonical form are always kept.
Usage:
  python3 gaussian_census.py SEED MAX_SECONDS [--walk-dir DIR] [--out DIR]
    SEED 11 and 12 were used (they continue the chains 2 and 1: chain 1 + SEED mod 2).
Input:  --walk-dir DIR   gaussian-walk-r16-s1.json, gaussian-walk-r16-s2.json (final states G, W of the chains;
                         default computation-outputs/C12 in the release assets)
Output: --out DIR        (default results/C12) census-seedSEED.jsonl (one row per step), census-seedSEED-grams.jsonl (kept Gram
                         matrices), census-seedSEED.state.json (state to resume; the same command resumes),
                         census-seedSEED.progress.json
Needs: Sage, gcc. Time: the run stops after MAX_SECONDS (about 0.8 s per step).

Z[i] 上の階数 16 の格子の pi 近傍のランダムウォークの調査：各歩で不変量を取る。

論文：第 1 の論文の計算 C12。一様な pi 近傍のランダムウォークの、先の 2 本の鎖を、終わりの状態から続ける。各歩で
rank32.gaussian_keys.analyse の不変量と標準形を記録する。この調査のランダムウォークの終わりの状態が、列挙の
2 本の鎖（code/produce/gaussian/gaussian_pi_neighbour_chain.py。後の版で加わる）の出発点である。
ルートのある格子は 30 歩に 1 個、標準形の無いルート無しの格子は全部、Gram 行列も残す。
使い方：
  python3 gaussian_census.py SEED MAX_SECONDS [--walk-dir DIR] [--out DIR]
    SEED は 11 と 12 を使った（鎖 2 と鎖 1 の続き。鎖 1 + SEED mod 2）。
入力：--walk-dir DIR   gaussian-walk-r16-s1.json、gaussian-walk-r16-s2.json（鎖の終わりの状態 G、W。
                       既定は Release のアセットの computation-outputs/C12）
出力：--out DIR        （既定は results/C12）census-seedSEED.jsonl（1 歩 1 行）、census-seedSEED-grams.jsonl（残した Gram 行列）、
                       census-seedSEED.state.json（再開用。同じコマンドで再開する）、census-seedSEED.progress.json
必要なもの：Sage、gcc。時間：MAX_SECONDS で止まる（1 歩約 0.8 秒）。
"""
import argparse
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from rank32.gaussian_keys import NON_DISCRETE, NOT_GENERATING, analyse, step  # noqa: E402
from rank32.io import asset, write_text_atomic as atomic  # noqa: E402


def census(seed, maxsec, walk_dir, out):
    """Continue the chain 1 + seed % 2 from its final state and record the invariants at every step.

    鎖 1 + seed % 2 を終わりの状態から続け、各歩の不変量を記録する。"""
    from sage.all import ZZ, matrix, pari
    pari.allocatemem(1 << 28, 1 << 30)
    out.mkdir(parents=True, exist_ok=True)
    res, prog, state, keep = (out / f"census-seed{seed}{x}" for x in (".jsonl", ".progress.json", ".state.json", "-grams.jsonl"))
    rng = random.Random(seed)
    if state.exists():
        s = json.loads(state.read_text())
        G, W = matrix(ZZ, s["G"]), matrix(ZZ, s["W"])
        rng.setstate((s["rng"][0], tuple(s["rng"][1]), s["rng"][2]))
    else:
        s0 = json.loads((walk_dir / f"gaussian-walk-r16-s{1 + seed % 2}.json").read_text())
        G, W = matrix(ZZ, s0["G"]), matrix(ZZ, s0["W"])
        s = dict(seed=seed, steps=0, elapsed=0.0)
    t0 = time.time()
    e0 = s["elapsed"]
    last = 0

    def save(st):
        s.update(G=[list(map(int, x)) for x in G.rows()], W=[list(map(int, x)) for x in W.rows()], elapsed=e0 + time.time() - t0)
        x = rng.getstate()
        s["rng"] = [x[0], list(x[1]), x[2]]
        atomic(state, json.dumps(s))
        atomic(prog, json.dumps(dict(job=f"census-seed{seed}", stage="walk and invariants", state=st, done=s["steps"], total=None, unit="steps", failures=0, elapsed_s=round(s["elapsed"], 1),
                                     updated=time.strftime("%Y-%m-%d %H:%M:%S %z")), ensure_ascii=False))

    while time.time() - t0 < maxsec:
        tick = time.process_time()
        G, W = step(G, W, rng)
        walk = time.process_time() - tick
        s["steps"] += 1
        Gl, Wl = [list(map(int, x)) for x in G.rows()], [list(map(int, x)) for x in W.rows()]
        a = analyse(Gl, Wl)
        a["sec"]["step"] = round(walk, 3)
        a["id"] = f"c{seed}-{s['steps']}"
        for k in ("hist", "rhist"):
            a[k] = a[k][:32]
        with res.open("a") as f:
            f.write(json.dumps(a, ensure_ascii=False) + "\n")
        if (a["roots"] and s["steps"] % 30 == 0) or (a["status"] in (NON_DISCRETE, NOT_GENERATING)):
            with keep.open("a") as f:
                f.write(json.dumps(dict(id=a["id"], roots=a["roots"], status=a["status"], gram=Gl, W=Wl)) + "\n")
        if time.time() - last > 60:
            save("running")
            last = time.time()
    save("time limit")
    print("steps", s["steps"], flush=True)


def main():
    ap = argparse.ArgumentParser(description="Census of the pi-neighbour walk (C12) / pi 近傍のランダムウォークの調査（C12）")
    ap.add_argument("seed", type=int)
    ap.add_argument("maxsec", type=float)
    ap.add_argument("--walk-dir", type=Path, default=asset("computation-outputs/C12"), help="final states of the chains / 鎖の終わりの状態")
    ap.add_argument("--out", type=Path, default=Path("results/C12"), help="output directory / 出力先")
    a = ap.parse_args()
    census(a.seed, a.maxsec, a.walk_dir, a.out)


if __name__ == "__main__":
    main()
