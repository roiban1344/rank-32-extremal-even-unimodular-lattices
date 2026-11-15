"""Generators of Aut F of the components of the involution gluing, their action on the discriminant group, and the orders of the image and the kernel.

Paper: computation C4 of the first paper (Table 11); the input table of the gluing of involutions.
For every component (rank32.discriminant_action.PARTS) it computes Aut F with qfauto (rank32.forms),
checks the order against the list of components, writes the action of the generators on D_F = F^#/F,
and computes with GAP the order of the image bar G_F, the kernel |Aut F| / |bar G_F| and the
sizes of the orbits of bar G_F on D_F. The generators depend on the basis chosen by qfauto; the
orders and the orbits do not.
Input:  --parts DIR   the lists of 2-elementary lattices of minimum 4 (output of C7; default
                      computation-outputs/C7 in the release assets)
Output: --out FILE    groups.json (default results/C4/groups.json), one entry per component (resumable:
                      components already present are skipped)
Needs: SageMath (PARI, GAP). Time: about a minute.

対合の貼り合わせの成分について、Aut F の生成元、その判別式群への作用、像とカーネルの位数を求める。

論文：第 1 の論文の計算 C4（表 11）。対合の貼り合わせの入力の表。
成分（rank32.discriminant_action.PARTS）ごとに、qfauto（rank32.forms）で Aut F を求めて位数を
成分の一覧と照らし、生成元の D_F = F^#/F への作用を書き、GAP で像 bar G_F の位数、カーネル
|Aut F| / |bar G_F|、bar G_F の D_F 上の軌道の大きさを求める。生成元は qfauto の選ぶ基底に依るが、
位数と軌道は依らない。
入力：--parts DIR   最小 4 の 2-初等的な格子の一覧（C7 の出力。既定は Release のアセットの
                    computation-outputs/C7）
出力：--out FILE    groups.json（既定は results/C4/groups.json）。成分ごとに 1 項目（再開できる：既にある
                    成分は飛ばす）
必要なもの：SageMath（PARI、GAP）。時間：約 1 分。
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix, libgap  # noqa: E402

from rank32 import forms  # noqa: E402
from rank32.discriminant_action import PARTS, action, disc_basis  # noqa: E402
from rank32.io import asset  # noqa: E402


def main(parts_dir, path):
    res = json.load(open(path)) if path.exists() else {}
    for name, (fn, i) in PARTS.items():
        if name in res:
            continue
        c = json.load(open(parts_dir / fn))["classes"][i]; G = matrix(ZZ, c["gram"]); t0 = time.time()
        order, g = forms.qfauto(G)
        assert order == c["aut"], (name, order, c["aut"])
        gens, J, Sinv = disc_basis(G); k = len(gens)
        imgs = [action(G, x, gens, J, Sinv) for x in g]
        # permutations of the whole discriminant group (0 is point 1, mask m is point m + 1)
        # 判別式全体への置換（0 を点 1、mask m を点 m + 1）
        n = 1 << k
        perms = []
        for im in imgs:
            t = np.zeros(n, dtype=np.int64)
            for l, y in enumerate(im):
                t[1 << l: 2 << l] = t[0: 1 << l] ^ y
            perms.append(t)
        P = [libgap.PermList([int(x) + 1 for x in t]) for t in perms]
        H = libgap.Group(P); img = int(H.Size())
        orbs = sorted(int(libgap.Size(o)) for o in libgap.Orbits(H, list(range(1, n + 1))))
        res[name] = dict(file=fn, index=i, s=k, aut=str(order), image=str(img), kernel=order // img, ngens=len(g), orbits=orbs, seconds=time.time() - t0,
                         gens_disc=imgs)
        print(name, "s", k, "aut", order, "image", img, "kernel", order // img, "orbits", orbs[:30], len(orbs), f"{time.time()-t0:.1f}s", flush=True)
        tmp = str(path) + ".tmp"; json.dump(res, open(tmp, "w")); Path(tmp).replace(path)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--parts", type=Path, default=asset("computation-outputs/C7"), help="directory of min4-twoel-*.json / min4-twoel-*.json のディレクトリ")
    ap.add_argument("--out", type=Path, default=Path("results/C4/groups.json"), help="output groups.json / 出力の groups.json")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    main(a.parts, a.out)
