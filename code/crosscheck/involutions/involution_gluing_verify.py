"""Checks of the results of the involution gluing.

Paper: computation C4 of the first paper (Table 11); the mode known uses the known lattices of C32.
Modes:
  recheck PAIR      pass every saved representative through the normal form again: phi is an
                    anti-isometry, the keys are distinct, |Stab| and the extremal flag agree with
                    the record, and the masses are M_up and M_ext. No search. Also writes the
                    sha256 of the sorted keys.
  gram PAIR N       build the Gram matrix of rank 32 of N extremal and N non-extremal classes
                    (random, seed 1): unimodular, even, minimum 4 exactly for the extremal ones
                    (then 146880 vectors of norm 4).
  stab PAIR N       for up to N extremal classes with distinct |Stab| (seed 2): |Aut L| with
                    rank32.bliss_automorphisms, and |C(sigma)| = |Stab| (GAP). Resumable.
  known PAIR...     take involutions of the 26 known lattices (random words in the generators,
                    seed 3), read (F, K, phi) back, and check that the class is in the list and
                    is extremal, with |C(sigma)| recomputed.
Input (by default in computation-outputs/C4, C7 and C32 of the release assets):
        --lists DIR (the lists <PAIR>.json and the checkpoints <PAIR>.ckpt.json of the gluing), the
        tables of rank32.involution_gluing (--groups, --groups-12-20, --scale, --parts),
        --known FILE (the 26 known lattices with generators, output of C32)
Output: verify-recheck-<PAIR>.json, verify-gram-<PAIR>.json, verify-stab-<PAIR>.json,
        verify-known-<first PAIR>.json in --out DIR (default results/C4)
Time: recheck seconds (6.7 hours for 16-16d8), known about 13 minutes.

対合の貼り合わせの結果の検算。

論文：第 1 の論文の計算 C4（表 11）。モード known は C32 の既知の格子を使う。
モード：
  recheck PAIR      保存された全部の代表を、もう一度正規形に通す。phi が反等長であること、キーが互いに
                    異なること、|Stab| と極値の判定が記録と同じこと、質量が M_up と M_ext であることを
                    確かめる。探索はしない。ソートしたキーの sha256 も書く。
  gram PAIR N       極値の類 N 個と極値でない類 N 個（乱択、種 1）から、階数 32 の Gram 行列を作る。
                    ユニモジュラー・偶で、最小 4 になるのがちょうど極値の類であること（そのとき
                    ノルム 4 は 146880 本）を確かめる。
  stab PAIR N       |Stab| の違う極値の類を N 個まで（種 2）。rank32.bliss_automorphisms で |Aut L| を
                    求め、|C(sigma)| = |Stab| を確かめる（GAP）。再開できる。
  known PAIR...     既知の 26 格子の対合（生成元の乱択の語、種 3）から (F, K, phi) を読み戻し、その類が
                    一覧にあって極値であることを、|C(sigma)| を計算し直して確かめる。
入力（既定は Release のアセットの computation-outputs/C4、C7、C32）：
      --lists DIR（貼り合わせの一覧 <PAIR>.json とチェックポイント <PAIR>.ckpt.json）、
      rank32.involution_gluing の表（--groups、--groups-12-20、--scale、--parts）、
      --known FILE（生成元つきの既知の 26 格子。C32 の出力）
出力：--out DIR（既定は results/C4）の中の verify-recheck-<PAIR>.json、verify-gram-<PAIR>.json、
      verify-stab-<PAIR>.json、verify-known-<最初の PAIR>.json
時間：recheck は数秒（16-16d8 は 6.7 時間）、known は約 13 分。
"""
import argparse
import hashlib
import json
import random
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix  # noqa: E402

from rank32 import gluing as glue  # noqa: E402
from rank32 import involution_gluing as gg  # noqa: E402
from rank32 import involution_gluing_checks as vv  # noqa: E402
from rank32.involution_gluing_checks import (aut_tool, build_gram, centralizer_order, class_table, element_order,  # noqa: E402
                                             glue_of, load, min_info)
from rank32.io import asset  # noqa: E402

KNOWN = asset("computation-outputs/C32/known-s-26.jsonl")


def code_commit():
    """The commit of the checked-out code, or "" when git is not installed or this is not a git checkout.

    取り出したコードのコミット。git が無いか、git の作業ツリーでなければ ""。"""
    try:
        return subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except FileNotFoundError:
        return ""


def mode_recheck(pair):
    _, _, F, K, gF, gK, res = load(pair)
    t = time.time()
    g = glue.Glue(F, K, gF["gens_disc"], gK["gens_disc"], int(gF["image"]), int(gK["image"]))
    kFK = gF["kernel"] * gK["kernel"]; keys = set(); mass = Fraction(0); ext = Fraction(0); bad = []
    for i, c in enumerate(res["classes"]):
        phi = np.array(c["phi"], dtype=np.int64)
        assert g.is_antiisometry(phi), ("not an anti-isometry", i)
        k, n = g.canon(phi)
        e = g.is_extremal(phi)
        if n * kFK != c["stab"] or e != c["extremal"] or k in keys:
            bad.append(i)
        keys.add(k); mass += Fraction(1, n * kFK); ext += Fraction(1, n * kFK) if e else 0
    assert not bad, ("classes that disagree with the record", bad[:10], len(bad))
    assert len(keys) == res["n_classes"] and str(mass) == res["M_up"] and str(ext) == res["M_ext"]
    r = dict(pair=pair, classes=len(keys), M_up=str(mass), M_ext=str(ext), keys_sha256=hashlib.sha256(b"".join(sorted(keys))).hexdigest(),
             commit=code_commit(), seconds=time.time() - t)
    gg.wr(Path(vv.OUT) / f"verify-recheck-{pair}.json", r)
    print(r, flush=True)
    return g


def mode_gram(pair, N):
    GF_, GK_, F, K, gF, gK, res = load(pair)
    rng = random.Random(1); cl = res["classes"]
    ext = [c for c in cl if c["extremal"]]; non = [c for c in cl if not c["extremal"]]
    out = []
    for c in rng.sample(ext, min(N, len(ext))) + rng.sample(non, min(N, len(non))):
        t = time.time(); G, s = build_gram(GF_, GK_, c["phi"]); m, r2, r4 = min_info(G)
        assert (m == 4) == c["extremal"], ("the extremal flag differs", c)
        if m == 4:
            assert r4 == 146880, r4                     # vectors of norm 4 (qfminim counts twice the number of pairs) / ノルム 4 の本数（qfminim の 2 倍の本数）
        out.append(dict(stab=c["stab"], extremal=c["extremal"], min=m, roots=r2, norm4=r4, seconds=time.time() - t))
        print(out[-1], flush=True)
    gg.wr(Path(vv.OUT) / f"verify-gram-{pair}.json", dict(pair=pair, checks=out))
    print(pair, "Gram check OK", len(out))


def mode_stab(pair, N):
    GF_, GK_, F, K, gF, gK, res = load(pair)
    rng = random.Random(2); ext = [c for c in res["classes"] if c["extremal"]]
    out = []
    # choose so that the values of |Stab| are spread / |Stab| の値がばらけるように選ぶ
    by = {}
    for c in ext:
        by.setdefault(c["stab"], []).append(c)
    keys = sorted(by); sel = []
    for k in keys:
        sel.append(rng.choice(by[k]))
    rng.shuffle(sel); sel = sorted(sel[:N], key=lambda c: c["stab"])
    # holding the list of 220000 classes while bliss runs took more than 4 GB (16-16d8); drop it after the choice.
    # A stopped run continues with the same choice.
    # 22 万類の一覧を持ったまま bliss を回すと 4 GB を超えた（16-16d8）。選んだあとは捨てる。止まったら同じ種の選び方で続きから
    del res, ext, by
    path = Path(vv.OUT) / f"verify-stab-{pair}.json"
    if path.exists():
        out = json.load(open(path))["checks"]
        assert [x["stab"] for x in out] == [c["stab"] for c in sel[:len(out)]], "the choice differs from the previous run"
    for c in sel[len(out):]:
        t = time.time(); G, s = build_gram(GF_, GK_, c["phi"])
        rec = aut_tool(G, f"{pair}-{c['stab']}")
        ord_, cen = centralizer_order(rec["generators"], s, G)
        assert str(ord_) == str(rec["order"]), (ord_, rec["order"])
        out.append(dict(stab=c["stab"], aut=str(ord_), centralizer=cen, ok=cen == c["stab"], seconds=time.time() - t))
        print(out[-1], flush=True)
        assert cen == c["stab"], ("|Stab| differs from the centralizer", cen, c["stab"])
        gg.wr(path, dict(pair=pair, checks=out))
    print(pair, "|Stab| = |C(sigma)| check OK", len(out))


def mode_known(pairs, known=KNOWN, per=12):
    """Collect involutions of the 26 known lattices at random (representatives of types), and check that those that fall in the given pairs appear among the classes, with |C(sigma)| = |Stab|.

    unmatched counts the involutions whose fixed lattice has the rank of F of a pair but which are
    isometric to no (F, K) of the given pairs. Given all pairs of a type, this checks that the list of components is
    complete (components with roots do not occur in extremal lattices).

    既知の 26 格子の対合を乱択で集め（型の代表）、渡した組に当たるものが結果の類に現れることと、中心化群の位数 = |Stab| を確かめる。

    unmatched は、固定格子の階数が組の F の階数なのに、渡した組のどの (F, K) にも等長でない対合の数。型の全部の組を渡せば
    成分の一覧の漏れの検査になる（ルートを持つ成分は極値の格子に現れない）。"""
    data = {}
    for pair in pairs:
        GF_, GK_, F, K, gF, gK, res = load(pair)
        g = glue.Glue(F, K, gF["gens_disc"], gK["gens_disc"], int(gF["image"]), int(gK["image"]))
        cl = class_table(pair, g, res, gF["kernel"] * gK["kernel"])
        del res
        data[pair] = (GF_, GK_, g, cl, gF["kernel"] * gK["kernel"])
    rng = random.Random(3); out = []
    path = Path(vv.OUT) / f"verify-known-{pairs[0]}.json"
    for line in open(known):
        r = json.loads(line); G = matrix(ZZ, r["gram"]); gens = [matrix(ZZ, x) for x in r["generators"]]
        order = int(r["order"])
        rows = {}
        for _ in range(60):
            w = matrix.identity(ZZ, 32)
            for _ in range(rng.randrange(10, 40)):
                w = w * gens[rng.randrange(len(gens))]
            o = element_order(w)
            if o % 2:
                continue
            x = w ** (o // 2)
            if x == -matrix.identity(ZZ, 32):
                continue
            d = int((x - 1).right_kernel().dimension())
            for xx, dd in ((x, d), (-x, 32 - d)):
                l = rows.setdefault(dd, [])
                if len(l) < 6 and xx not in l:
                    l.append(xx)
        found = []; unmatched = {}
        for d in sorted({v[0].nrows() for v in data.values()}):
            for x_ in rows.get(d, []):
                hit = False
                for pair, (GF_, GK_, g, cl, kFK) in data.items():
                    if GF_.nrows() != d:
                        continue
                    phi = glue_of(G, x_, GF_, GK_)
                    if phi is None:
                        continue
                    hit = True
                    assert g.is_antiisometry(phi), ("not an anti-isometry", r["id"], pair)
                    k, n = g.canon(phi)
                    c = cl.get(k)
                    ok = c is not None and bool(c["extremal"])
                    found.append(dict(pair=pair, in_list=c is not None, extremal=ok, stab=None if c is None else int(c["stab"]), stab_recomputed=int(n) * kFK))
                    print(r["id"], pair, found[-1], flush=True)
                    assert ok, ("an involution of a known lattice is not in the list", r["id"], pair)
                # for pairs with F != K, an involution of the opposite orientation (-x matches) is not counted as missing
                # F != K の組は向きが逆の対合（-x が当たる）を「漏れ」に数えない
                if not hit and not any(F_.nrows() == d and glue_of(G, x_, K_, F_) is not None for (F_, K_, *_) in data.values() if F_ != K_):
                    unmatched[str(d)] = unmatched.get(str(d), 0) + 1
        out.append(dict(id=r["id"], order=r["order"], dims=sorted(rows), found=found, unmatched=unmatched))
        print(r["id"], "order", r["order"], "fixed ranks of involutions", sorted(rows), "hits", len(found), "involutions matching no pair", unmatched, flush=True)
        gg.wr(path, dict(pairs=pairs, lattices=out))
    print("check of the known lattices OK,", sum(1 for x in out if x["found"]), "with hits")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mode", choices=["recheck", "gram", "stab", "known"], help="mode / モード")
    ap.add_argument("args", nargs="+", help="PAIR [N] or PAIR... / PAIR [N] または PAIR...")
    ap.add_argument("--lists", type=Path, default=vv.LISTS, help="directory of the lists <PAIR>.json / 一覧 <PAIR>.json のディレクトリ")
    ap.add_argument("--out", type=Path, default=gg.OUT, help="directory of the outputs / 出力のディレクトリ")
    ap.add_argument("--groups", type=Path, default=gg.GROUPS, help="groups.json")
    ap.add_argument("--groups-12-20", type=Path, default=gg.GROUPS_12_20, help="table of the pairs (12, 20) / (12, 20) の組の表")
    ap.add_argument("--scale", type=Path, default=gg.SCALE_FILE, help="orders |O(q)| / |O(q)| の表")
    ap.add_argument("--parts", type=Path, default=gg.PARTS_DIR, help="directory of min4-twoel-*.json / min4-twoel-*.json のディレクトリ")
    ap.add_argument("--known", type=Path, default=KNOWN, help="the 26 known lattices / 既知の 26 格子")
    a = ap.parse_args()
    gg.load_tables(a.groups, a.groups_12_20, a.scale, a.parts)
    vv.LISTS, vv.OUT = a.lists, a.out
    a.out.mkdir(parents=True, exist_ok=True)
    if a.mode == "recheck":
        mode_recheck(a.args[0])
    elif a.mode == "gram":
        mode_gram(a.args[0], int(a.args[1]))
    elif a.mode == "stab":
        mode_stab(a.args[0], int(a.args[1]))
    elif a.mode == "known":
        mode_known(a.args, a.known)
