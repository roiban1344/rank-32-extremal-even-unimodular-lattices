"""The 2-neighbour method that walks only orbit representatives; completeness is decided by the mass.

Paper: computation C2 of the first paper (Table 11), lists of genera for the types of order 3; the mode
hermr gives the rootless classes of the Hermitian genus (13, 5) used later (C16). The functions are in
rank32.orbit_walk.
Modes:
  real P N S [SECONDS]   the genus of even P-elementary Z-lattices of rank N and determinant P^S (odd
                         discriminant): one neighbour for each v in L/2L. Done when the mass of the genus is reached.
  herm Z [SECONDS]       the rootless part of the unimodular Z[omega]-lattices of rank Z: two neighbours for each
                         isotropic line for the inert prime 2. Done when the rootless mass is reached.
  hermr Z S [SECONDS]    only the rootless classes of the Hermitian genus (Z, S). The rootless mass is not known,
                         so completeness cannot be decided (a search for a lower bound and for examples).
  hermg Z S [SECONDS]    the whole Hermitian genus (Z, S) (also the classes with roots). Done when mass_H(Z, S) is reached.
SECONDS (default 3000) is the time limit: the walk stops when it is reached, even if the mass is not.
Isometry: first the group order, then a short qfisom, then (if the shortest vectors generate) an isometry
from the canonical forms of bliss, then a longer qfisom, finally the coarse canonical form of the vectors of
norm <= 6. What remains undecided is recorded as undecided (the mass then does not close).
The random seed is the start time, so the order of the classes and the Gram matrices differ from run to run;
the list (group orders, invariants, mass) does not.
Input:  --input-dir DIR (hermr: unimod-z<Z>.json; real with P = 3: unimod-z6.json for K_12 as a block, if present)
Output: --out DIR, genus-p<P>-n<N>-s<S>.json, unimod-z<Z>.json, hermr-z<Z>-s<S>.json or hermg-z<Z>-s<S>.json,
        updated atomically during the run (the same command resumes), and the same name with
        .progress.json (every minute)
Needs:  SageMath (PARI), bliss (rank32.canonical). The memory limit is left to the caller (the soft limit
        of RLIMIT_AS is raised to the hard limit). Time: seconds to hours depending on the genus.

軌道の代表だけを歩く 2-近傍法。完了は質量の一致で判定する。

論文：第 1 の論文の計算 C2（表 11）。位数 3 の型の種の一覧。モード hermr は、後の計算（C16）が使う
エルミートの種 (13, 5) のルート無しの類を与える。関数は rank32.orbit_walk にある。
モード：
  real P N S [秒]   偶・P-基本・階数 N・行列式 P^S の Z-格子の種（判別式は奇数）。v in L/2L ごとに近傍が 1 個。
                    種の質量に届けば完了。
  herm Z [秒]       階数 Z のユニモジュラー Z[omega]-格子のルート無しの部分。不活性素数 2 の等方直線ごとに
                    近傍が 2 個。ルート無しの質量に届けば完了。
  hermr Z S [秒]    エルミート種 (Z, S) のルート無しの類だけを歩く。ルート無しの部分の質量は分からないので、
                    完了の判定はできない（下界と実例を得るための探索）。
  hermg Z S [秒]    エルミート種 (Z, S) の全体（ルートのある類も含む）。mass_H(Z, S) に届けば完了。
秒（既定 3000）は時間の上限である。質量に届いていなくても、そこで歩くのをやめる。
同型の判定：群の位数で絞る、短い qfisom、最短ベクトルが生成するなら bliss の標準形からの等長写像、長い
qfisom、最後にノルム 6 以下のベクトルの粗い標準形。それでも決まらないものは未判定として記録する
（そのとき質量は閉じない）。
乱数の種は開始時刻なので、類の順と Gram 行列は回ごとに違う。一覧（群の位数・不変量・質量）は変わらない。
入力：--input-dir DIR（hermr：unimod-z<Z>.json。real で P = 3 のとき、あれば K_12 をブロックにするための unimod-z6.json）
出力：--out DIR の genus-p<P>-n<N>-s<S>.json、unimod-z<Z>.json、hermr-z<Z>-s<S>.json、hermg-z<Z>-s<S>.json の
      どれか。実行中にアトミックに更新する（同じコマンドで再開する）。同名の .progress.json（1 分ごと）
必要なもの：SageMath（PARI）、bliss（rank32.canonical）。メモリの上限は呼び出し側に任せる（RLIMIT_AS の
ソフトリミットをハードリミットまで上げる）。時間：種によって数秒から数時間。
"""
import argparse
import gc
import json
import random
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

from sage.all import ZZ, matrix, pari, block_diagonal_matrix  # noqa: E402

from rank32 import eisenstein_overlattices as ptov  # noqa: E402
from rank32 import orbit_walk as w  # noqa: E402
from rank32.io import asset, write_json_atomic_compact  # noqa: E402

_, _hard = resource.getrlimit(resource.RLIMIT_AS)
resource.setrlimit(resource.RLIMIT_AS, (_hard, _hard))  # the memory limit is left to the caller / メモリの上限は呼び出し側に任せる

DEFAULT_INPUT = str(asset("computation-outputs/C2"))
DEFAULT_OUT = "results/C2"
qfsafe, forms, counts, autgraph = w.qfsafe, w.forms, w.counts, w.autgraph
A2, W2 = w.A2, w.W2


def run(mode, args, src_dir, out):
    # a cyclic GC that frees an old Gen of PARI inside sig_on gives a SystemError; GC runs only between steps
    # 循環 GC が PARI の sig_on の中で古い Gen を解放すると SystemError になる。GC は区切りでだけ走らせる
    gc.disable()
    out.mkdir(parents=True, exist_ok=True)
    ptmass, ptroot = w.ptmass, w.ptroot
    if mode == "real":
        p, n, s = map(int, args[0:3])
        maxsec = float(args[3]) if len(args) > 3 else 3000
        target = Fraction(str(ptmass.mass_Z_genus(p, n, s)))
        path = out / f"genus-p{p}-n{n}-s{s}.json"
    elif mode in ("hermg", "hermr"):
        p, z, s = 3, int(args[0]), int(args[1])
        n = 2 * z
        maxsec = float(args[2]) if len(args) > 2 else 3000
        target = Fraction(str(ptmass.mass_H(3, z, s)))
        path = out / f"{mode}-z{z}-s{s}.json"
    else:
        p, z = 3, int(args[0])
        n = 2 * z
        maxsec = float(args[1]) if len(args) > 1 else 3000
        target = Fraction(str(ptroot.noroot_masses(3, z)[z]))
        path = out / f"unimod-z{z}.json"
    prog = path.with_suffix(".progress.json")
    t0 = time.time()
    rng = random.Random(int(t0))
    classes, byinv, undecided = [], {}, []
    herm = mode in ("herm", "hermg", "hermr")
    rootless_only = mode in ("herm", "hermr")

    def mass():
        return sum((Fraction(1, c["aut"]) for c in classes if (not rootless_only or c["inv"][0] == 0)), Fraction(0))

    def add(G, W):
        inv = counts(G)
        if rootless_only and inv[0] > 0 and any(c["inv"][0] == 0 for c in classes):
            return False
        G, W = w.short_basis(G, W, rng)
        # order of the tests: group order -> short qfisom -> isometry from the canonical forms if the shortest vectors
        # generate -> long qfisom -> coarse canonical form. The dense graph of the shell of norm 6 is not built (more than 8 GB).
        # The classes with roots have large groups on which PARI tends to stall, so the canonical form (bliss) is tried first.
        # 判定の順序：群の位数で絞る -> 短い qfisom -> 最短ベクトルが生成するなら標準形の写像で確定 -> 長い qfisom -> 粗い標準形。
        # ノルム 6 の殻の密なグラフは作らない（8 GB を超えた）。ルートのある類は群が大きく PARI が止まりやすいので、標準形（bliss）を先に試す
        canon = w.safe_canonical(G, W) if inv[0] > 0 else None
        try:
            order = int(canon["order"]) if canon else qfsafe.qfauto(forms(G, W), tries=20)[0]
        except qfsafe.QfTimeout:
            undecided.append(dict(inv=list(inv), why="qfauto", G=[list(map(int, x)) for x in G.rows()], W=None if W is None else [list(map(int, x)) for x in W.rows()]))
            print(f"undecided (qfauto): {inv}", flush=True)
            return False
        for j in byinv.get(inv, []):
            c = classes[j]
            if c["aut"] != order:
                continue
            if not canon:
                r = w.quick_isom(c, G, W, 8)
                if r:
                    return False
                if r is False:
                    continue
            if canon is None:
                canon = w.safe_canonical(G, W)
            if c.get("canon") is None:
                c["canon"] = w.safe_canonical(c["G"], c["W"])
            if canon and c["canon"]:
                if autgraph.isometry(c["canon"], c["G"], canon, G, c["W"], W) is not None:
                    return False
                print(f"separated by the canonical form (same order {order}): {inv}", flush=True)
                continue
            try:
                if qfsafe.qfisom(forms(c["G"], c["W"]), forms(G, W), tries=4, cap=60) is not None:
                    return False
                continue
            except qfsafe.QfTimeout:
                pass
            if autgraph.shell_hash(G, W, bound=6, keep_min=4) != autgraph.shell_hash(c["G"], c["W"], bound=6, keep_min=4):
                print(f"separated by the coarse canonical form up to norm 6: {inv}", flush=True)
                continue
            undecided.append(dict(inv=list(inv), why="qfisom", G=[list(map(int, x)) for x in G.rows()], W=None if W is None else [list(map(int, x)) for x in W.rows()]))
            print(f"undecided (qfisom): {inv}", flush=True)
            return False
        byinv.setdefault(inv, []).append(len(classes))
        classes.append(dict(G=G, W=W, inv=inv, aut=order, done=False, canon=canon))
        return True

    def save(state):
        rows = []
        for c in classes:
            r = dict(aut=str(c["aut"]), inv=list(c["inv"]), walked=bool(c["done"]))
            if herm:
                r.update(G=[list(map(int, x)) for x in c["G"].rows()], W=[list(map(int, x)) for x in c["W"].rows()])
            else:
                G = c["G"]
                q = pari(G)
                mn = int(q.qfminim()[1])
                r.update(gram=[list(map(int, x)) for x in G.rows()], min=mn, min_dual=int(pari((p * G.inverse()).change_ring(ZZ)).qfminim()[1]),
                         norm4=c["inv"][1] - c["inv"][0])
            rows.append(r)
        m = mass()
        d = dict(target=str(target), mass=str(m), done=m == target, classes=rows, undecided=undecided)
        d.update(dict(z=z) if rootless_only else dict(z=z, s=s) if herm else dict(p=p, n=n, s=s))
        write_json_atomic_compact(path, d)
        write_json_atomic_compact(prog, dict(job=path.stem, stage="neighbours from orbit representatives", state=state, processed=len(classes), unit="classes",
                                             walked=sum(1 for c in classes if c["done"]), undecided=len(undecided), mass=float(m), target=float(target),
                                             elapsed_s=time.time() - t0, updated=time.strftime("%Y-%m-%d %H:%M:%S %z")))

    if path.exists():
        data = json.loads(path.read_text())
        for r in data["classes"]:
            G = matrix(ZZ, r["G"] if herm else r["gram"])
            W = matrix(ZZ, r["W"]) if herm else None
            inv = tuple(r["inv"]) if "inv" in r else counts(G)
            byinv.setdefault(inv, []).append(len(classes))
            classes.append(dict(G=G, W=W, inv=inv, aut=int(r["aut"]), done=bool(r.get("walked"))))
    if mode == "hermr":
        if not classes:
            # start: from a rootless U, add isotropic lines that meet no vector of norm 2, one after another
            # 出発点：ルート無しの U から、ノルム 2 に当たらない等方直線を順に足す
            src = [c for c in json.loads((src_dir / f"unimod-z{z}.json").read_text())["classes"] if c["inv"][0] == 0]
            while not classes:
                c0 = rng.choice(src)
                G0, W0 = matrix(ZZ, c0["G"]), matrix(ZZ, c0["W"])
                for _ in range((z - s) // 2):
                    D, iso, good = ptov.good_lines(G0)
                    if not good:
                        G0 = None
                        break
                    G0, W0 = ptov.overlattice(G0, W0, D, rng.choice(sorted(good)))
                if G0 is not None:
                    assert ZZ(G0.det()) == 3 ** s and counts(G0)[0] == 0
                    add(G0, W0)
    elif rootless_only:
        if not classes:
            add(block_diagonal_matrix([A2] * z, subdivide=False), block_diagonal_matrix([W2] * z, subdivide=False))
    elif herm:
        if not classes:
            # start: a sum of Eisenstein structures of E8 and A_2^s (the roots generate, so bliss is fast); only when z - s is a multiple of 4
            # 出発点：E8 の Eisenstein 構造の直和と A_2^s（ルートが格子を生成するので bliss が速い）。(z - s) が 4 の倍数のときだけ
            assert (z - s) % 4 == 0, "cannot build a starting lattice"
            Ge, We = block_diagonal_matrix([A2] * 4, subdivide=False), block_diagonal_matrix([W2] * 4, subdivide=False)
            for _ in range(2):
                D, iso, _g = ptov.good_lines(Ge)
                Ge, We = ptov.overlattice(Ge, We, D, sorted(iso)[0])
            assert Ge.det() == 1
            k = (z - s) // 4
            G0 = block_diagonal_matrix([Ge] * k + [A2] * s, subdivide=False)
            W0 = block_diagonal_matrix([We] * k + [W2] * s, subdivide=False)
            assert ZZ(G0.det()) == 3 ** s, "the starting lattice has the wrong s"
            add(G0, W0)
    else:
        for G in w.real_seeds(p, n, s, src_dir):
            add(G, None)
    print(f"start: {len(classes)} classes, mass {float(mass()):.8g} / {float(target):.8g}", flush=True)
    save("start")
    if mode == "herm":
        # fast stage: look only at the invariants of the neighbours and pick up lattices with new invariants first
        # (no isometry tests; fast when what remains is one class of high symmetry)
        # 高速なステージ：近傍の不変量だけを見て、未知の不変量の格子を先に拾う（同型判定を飛ばす。残りが対称性の高い 1 類のときに速い）
        for c in sorted([c for c in classes if c["inv"][0] == 0], key=lambda t: -t["aut"]):
            if mass() >= target or time.time() - t0 >= maxsec:
                break
            try:
                _, gens = qfsafe.qfauto(forms(c["G"], c["W"]), tries=20)
            except qfsafe.QfTimeout:
                continue
            reps, total = w.isotropic_reps(c["G"], gens + [c["W"].transpose()], rng, 30000)
            seen = 0
            for x in reps:
                gc.collect()
                for Gn, Wn in w.herm_neighbours(c["G"], c["W"], w.to_vec(x, n)):
                    seen += 1
                    inv = counts(Gn)
                    if inv[0] == 0 and inv not in byinv and add(Gn, Wn):
                        d = classes[-1]
                        print(f"class {len(classes)} (new invariants): {d['inv']}, |Aut| {d['aut']}, mass {float(mass()):.8g} / {float(target):.8g} ({time.time() - t0:.0f} s)", flush=True)
                if mass() >= target or time.time() - t0 >= maxsec:
                    break
            print(f"fast stage: |Aut| {c['aut']}, representatives {len(reps)} / {total}, neighbours {seen} ({time.time() - t0:.0f} s)", flush=True)
            save("fast stage")
    last = time.time()
    cap = 3000  # the largest number of neighbours built from one class in one round / 1 類から 1 周で作る近傍の上限
    sweep = 0
    turn = 0
    while mass() < target and time.time() - t0 < maxsec:
        pool = [c for c in classes if not c["done"] and (not rootless_only or c["inv"][0] == 0 or not any(d["inv"][0] == 0 for d in classes))]
        if not pool:
            sweep += 1
            if all(c["done"] == "all" for c in classes):
                break
            for c in classes:
                if c["done"] != "all":
                    c["done"] = False
            continue
        # alternate: odd turns exhaust the class with the largest group (so that rare classes are not missed),
        # even turns take a few from the classes with small groups (most of the mass is there)
        # 交互に歩く：奇数回は群の大きい類を尽くす（珍しい類を漏らさない）、偶数回は群の小さい類から少しずつ（質量の大半はこちら）
        turn += 1
        small_turn = mode != "herm" and turn % 2 == 0
        if small_turn:
            c = rng.choice(sorted(classes, key=lambda t: t["aut"])[:max(5, len(classes) // 3)])
        else:
            c = max(pool, key=lambda t: t["aut"])
        if c["inv"][0] > 0 and c.get("canon") is None:
            c["canon"] = w.safe_canonical(c["G"], c["W"])
        if c["inv"][0] > 0 and c.get("canon"):
            order, gens = int(c["canon"]["order"]), autgraph.generators(c["canon"], c["G"])
        else:
            try:
                order, gens = qfsafe.qfauto(forms(c["G"], c["W"]), tries=20)
            except qfsafe.QfTimeout:
                print(f"qfauto ran out of time: |Aut| {c['aut']}; this class is skipped", flush=True)
                c["done"] = True
                continue
        assert order == c["aut"], "the group order differs from the record"
        for T in gens:
            assert T.transpose() * c["G"] * T == c["G"]
        if herm:
            gens.append(c["W"].transpose())
        reps, total = w.isotropic_reps(c["G"], gens, rng, 200 if small_turn else cap)
        new = 0
        print(f"|Aut| {c['aut']}: after the orbits {w.rss()}, representatives {total} ({time.time() - t0:.0f} s)", flush=True)
        for x in reps:
            gc.collect()
            v = w.to_vec(x, n)
            for Gn, Wn in (w.herm_neighbours(c["G"], c["W"], v) if herm else w.real_neighbour(c["G"], v)):
                if add(Gn, Wn):
                    new += 1
                    d = classes[-1]
                    print(f"class {len(classes)}: invariants {d['inv']}, |Aut| {d['aut']}, mass {float(mass()):.8g} / {float(target):.8g} ({time.time() - t0:.0f} s, {w.rss()})", flush=True)
            if time.time() - last > 60:
                save("running")
                last = time.time()
            if mass() >= target or time.time() - t0 >= maxsec:
                break
        if not small_turn or len(reps) == total:
            c["done"] = "all" if total <= cap and len(reps) == total else True
        print(f"class with |Aut| {c['aut']}: {len(reps)} of {total} orbit representatives, {new} new classes ({time.time() - t0:.0f} s)", flush=True)
    assert mass() <= target, "the mass exceeds the target: an error in the isometry test or in a group order"
    save("complete" if mass() == target else "incomplete")
    m = mass()
    print(f"RESULT {path.stem}: {len(classes)} classes, mass {'agrees' if m == target else 'short by ' + str(target - m)}, undecided {len(undecided)}", flush=True)
    if herm:
        for c in classes:
            if c["inv"][0] == 0:
                print(f"  |U| = {c['aut']}, norm <= 4: {c['inv'][1]}, norm <= 6: {c['inv'][2]}", flush=True)
    else:
        for c in classes:
            if c["inv"][0] == 0:
                print(f"  |Aut| = {c['aut']}, norm 4: {c['inv'][1]}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mode", choices=["real", "herm", "hermr", "hermg"], help="the mode / モード")
    ap.add_argument("args", nargs="+", help="the arguments of the mode (see above) / モードの引数（上を参照）")
    ap.add_argument("--input-dir", default=DEFAULT_INPUT, help="input directory / 入力のディレクトリ")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    run(a.mode, a.args, Path(a.input_dir), Path(a.out))


if __name__ == "__main__":
    main()
