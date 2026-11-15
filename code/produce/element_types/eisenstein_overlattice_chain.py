"""For p = 3: from the rootless unimodular Z[omega]-lattices U of rank z, build the Z[omega]-lattices Z of minimum >= 4, for s = z, z - 2, ... in turn.

Paper: computation C2 of the first paper (Table 11), the lists of the cyclic parts Z of the types 3-(z, d)-s.
Method (rank32.eisenstein_overlattices): from each class of determinant 3^s, the isotropic lines whose class has no
vector of norm 2 are split into orbits of U(Z) (qfauto on the two forms); each representative gives an
overlattice of determinant 3^(s - 2); the new lattices are kept up to isometry (same group order, then an isometry
from the canonical forms of bliss when the shortest vectors generate, otherwise qfisom).
Check (the mass identity), at every step:
    sum_{Z'} N_s / |U(Z')| = sum_Z g(Z) / |U(Z)|,
where g(Z) is the number of isotropic lines without vectors of norm 2, and N_s = (3^(z - s + 2) - 1)/2 is
the number of lines of Z'/Z'^* (a symplectic space of dimension z - s + 2; Z^* as in Lemma 4.4 of the first paper).
Usage:  python3 eisenstein_overlattice_chain.py Z SMIN
Input:  --input-dir DIR with unimod-z<Z>.json (the list of unimodular lattices U; only the rootless ones are used)
Output: --out DIR, overlat-z<Z>.json (the classes for each s, in the list levels, and the mass checks), and
        overlat-z<Z>-s<s>.ckpt.json (the checkpoint of the s being processed; the same command resumes)
Needs:  SageMath (PARI), bliss (rank32.canonical). Time: z = 9, 10, 11 take 2 to 25 seconds; the larger z were not timed.

p = 3：階数 z のルートを持たないユニモジュラー Z[omega]-格子 U から、最小 >= 4 の Z[omega]-格子 Z を
s ごとに順に作る。

論文：第 1 の論文の計算 C2（表 11）。型 3-(z, d)-s の巡回部分 Z の一覧。
方法（rank32.eisenstein_overlattices）：行列式 3^s の各類について、類にノルム 2 のベクトルを持たない等方直線を
U(Z) の軌道に分ける（2 つの形式での qfauto）。各代表から行列式 3^(s - 2) の上格子を作り、等長を除いて残す
（群の位数が同じものについて、最短ベクトルが生成するなら bliss の標準形からの等長写像、そうでなければ qfisom）。
検算（質量の恒等式）を各段階で行う：
    sum_{Z'} N_s / |U(Z')| = sum_Z g(Z) / |U(Z)|
g(Z) はノルム 2 のベクトルを持たない等方直線の数、N_s = (3^(z - s + 2) - 1)/2 は Z'/Z'^*（Z^* は第 1 の論文の補題 4.4 のもの。z - s + 2 次元の
斜交空間）の直線の数。
使い方：python3 eisenstein_overlattice_chain.py Z SMIN
入力：--input-dir DIR の unimod-z<Z>.json（ユニモジュラー格子 U の一覧。ルートを持たないものだけを使う）
出力：--out DIR の overlat-z<Z>.json（各 s の類（リスト levels）と質量の検算）と overlat-z<Z>-s<s>.ckpt.json（処理中の s の
      チェックポイント。同じコマンドで再開する）
必要なもの：SageMath（PARI）、bliss（rank32.canonical）。時間：z = 9, 10, 11 は 2〜25 秒。それより大きい z は計っていない。
"""
import argparse
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix, pari  # noqa: E402

from rank32 import canonical as autgraph  # noqa: E402
from rank32.eisenstein_overlattices import forms, good_lines, orbits, overlattice, qfsafe  # noqa: E402
from rank32.io import asset, write_json_atomic_compact  # noqa: E402

DEFAULT_INPUT = str(asset("computation-outputs/C2"))
DEFAULT_OUT = "results/C2"


def run(z, smin, src_dir, out):
    """Build the lists for s = z, z - 2, ..., down to smin (or until a list is empty); write overlat-z<z>.json in out.

    s = z, z - 2, ... の一覧を smin まで（または一覧が空になるまで）作り、out に overlat-z<z>.json を書く。"""
    out.mkdir(parents=True, exist_ok=True)
    src = json.loads((src_dir / f"unimod-z{z}.json").read_text())
    assert src["done"], "the classification of U is not complete"
    cur = []
    assert not src.get("undecided"), "U has undecided cases"
    for r in src["classes"]:
        if r["inv"][0] == 0:
            cur.append(dict(G=matrix(ZZ, r["G"]), W=matrix(ZZ, r["W"]), aut=int(r["aut"])))
    s = z
    result = {"z": z, "levels": []}
    t0 = time.time()
    while True:
        rec = dict(s=s, classes=[])
        next_cl = []
        sum_g = Fraction(0)
        # resume: within one s the state is saved after each U / 途中再開：1 つの s の中で U ごとに保存する
        ck = out / f"overlat-z{z}-s{s}.ckpt.json"
        start = 0
        if ck.exists():
            d = json.loads(ck.read_text())
            start, sum_g, rec = d["done"], Fraction(d["sum_g"]), d["rec"]
            next_cl = [dict(G=matrix(ZZ, e["G"]), W=matrix(ZZ, e["W"]), aut=int(e["aut"])) for e in d["next"]]
            print(f"resume: s={s} from item {start}, next level {len(next_cl)} classes", flush=True)
        for ci, c in enumerate(cur):
            if ci < start:
                continue
            G, W = c["G"], c["W"]
            q = pari(G)
            n4 = int(q.qfminim(4, 0)[0])
            info = dict(aut=str(c["aut"]), norm4=n4, min_dual3=int(pari((3 * G.inverse()).change_ring(ZZ)).qfminim()[1]),
                        gram=[list(map(int, r)) for r in G.rows()], omega=[list(map(int, r)) for r in W.rows()])
            if s - 2 >= smin:
                D, iso, good = good_lines(G)
                info.update(iso_lines=len(iso), good_lines=len(good))
                sum_g += Fraction(len(good), c["aut"])
                if good:
                    _, gens = qfsafe.qfauto(forms(G, W), tries=30)
                    for T in gens:
                        assert T.transpose() * G * T == G and T * W.transpose() == W.transpose() * T
                    reps = orbits(D, gens, good)
                    info["orbits"] = [int(x[1]) for x in reps]
                    for key, size in reps:
                        Gn, Wn = overlattice(G, W, D, key)
                        assert int(pari(Gn).qfminim()[1]) >= 4
                        new = True
                        aut = qfsafe.qfauto(forms(Gn, Wn), tries=30)[0]
                        cn = None
                        for e in next_cl:
                            if e["aut"] != aut:
                                continue
                            # if the shortest vectors generate, an isometry from the canonical forms decides (without PARI)
                            # 最短ベクトルが生成するなら標準形から写像を作って確定する（PARI を使わない）
                            if cn is None:
                                cn = autgraph.canonical(Gn, Wn, maxvec=4000) or False
                            if e.get("canon") is None:
                                e["canon"] = autgraph.canonical(e["G"], e["W"], maxvec=4000) or False
                            if cn and e["canon"]:
                                if autgraph.isometry(e["canon"], e["G"], cn, Gn, e["W"], Wn) is not None:
                                    new = False
                                    break
                                continue
                            if qfsafe.qfisom(forms(e["G"], e["W"]), forms(Gn, Wn), tries=12, cap=120) is not None:
                                new = False
                                break
                        if new:
                            next_cl.append(dict(G=Gn, W=Wn, aut=aut))
            rec["classes"].append(info)
            write_json_atomic_compact(ck, dict(done=ci + 1, sum_g=str(sum_g), rec=rec, next=[dict(G=[list(map(int, r)) for r in e["G"].rows()], W=[list(map(int, r)) for r in e["W"].rows()], aut=str(e["aut"])) for e in next_cl]))
            print(f"z={z} s={s}: |U|={c['aut']} norm4={n4} min sqrt3 Z^#={info['min_dual3']} isotropic lines {info.get('iso_lines')} good lines {info.get('good_lines')} orbits {len(info.get('orbits', []))} ({time.time() - t0:.0f} s)", flush=True)
        # check / 検算
        if s - 2 >= smin:
            Ns = (3 ** (z - (s - 2)) - 1) // 2
            lhs = sum((Fraction(Ns, e["aut"]) for e in next_cl), Fraction(0))
            ok = lhs == sum_g
            rec["mass_check_next"] = dict(lhs=str(lhs), rhs=str(sum_g), ok=ok)
            print(f"  mass check s={s} -> {s - 2}: {lhs} = {sum_g}: {'agrees' if ok else 'differs'}", flush=True)
            assert ok
        result["levels"].append(rec)
        write_json_atomic_compact(out / f"overlat-z{z}.json", result)
        if s - 2 < smin or not next_cl:
            if not next_cl and s - 2 >= smin:
                print(f"z={z}: no lattice of minimum >= 4 at s = {s - 2}", flush=True)
            break
        cur = next_cl
        s -= 2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("z", type=int, help="the rank over Z[omega] / Z[omega] 上の階数")
    ap.add_argument("smin", type=int, help="the smallest s to build / 作る最小の s")
    ap.add_argument("--input-dir", default=DEFAULT_INPUT, help="directory with unimod-z<z>.json / unimod-z<z>.json のディレクトリ")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    run(a.z, a.smin, Path(a.input_dir), Path(a.out))


if __name__ == "__main__":
    main()
