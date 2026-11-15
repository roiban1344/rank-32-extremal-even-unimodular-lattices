"""The gluing of involutions: all double cosets of anti-isometries phi of a pair of components, found by walking through reflection neighbours until the mass is complete.

Used by computation C4 of the first paper (Table 11).
An involution sigma of an extremal lattice L of rank 32 gives the components F = L^sigma (the fixed
lattice) and K = L^{-sigma} (the anti-fixed lattice), 2-elementary lattices of minimum >= 4, and an
anti-isometry phi: D_F -> D_K. The lattices L with given (F, K) correspond to the double cosets
bar G_K \\ O / bar G_F of the anti-isometries. The key of a double coset is its normal form canon()
(rank32.gluing). New double cosets are found by phi -> phi o t_v (t_v a reflection of D_F). The search
stops when sum 1/|Stab| (in the image) reaches |O(q)| / (|bar G_F| |bar G_K|): then every double coset has
been found, whatever the search did.
The search: an exhaustive test of the neighbours (classes with a large |Stab| first, one
representative per orbit of Stab), an ordinary random walk with three walkers (Metropolis with
beta = 0, 1/2, 1), and an equivariant walk (moves by reflections that commute with an element of
prime order of the stabilizer), in time shares 1 : 1 : 2. The random seed is
20261031 + (time mod 100000); the representatives and the order of discovery change from run to
run, the number of classes, |Stab| and the masses do not.
Tables: load_tables() reads groups.json (the output of produce/involutions/discriminant_action.py),
optionally the table of the pairs (12, 20), the orders |O(q)| (involution-scale.json) and the lists
of components. part() and order_O() load the default paths if load_tables() was not called.
Paths: the tables are read by default from computation-outputs/C4 and C7 of the release assets;
search() writes to OUT, by default results/C4 under the current directory.

対合の貼り合わせ：成分の組の反等長 phi の両側剰余類の全部を、質量が満ちるまで反射の近傍を渡り歩いて求める。

第 1 の論文の計算 C4（表 11）が使う。
階数 32 の極値格子 L の対合 sigma から、成分 F = L^sigma（固定格子）と K = L^{-sigma}（反固定格子）と、
反等長 phi: D_F -> D_K が決まる。F と K は、どちらも最小 4 以上の 2-初等的な格子である。
(F, K) を与えたときの L は、反等長の両側剰余類 bar G_K \\ O / bar G_F に対応する。両側剰余類のキーは正規形 canon()（rank32.gluing）である。
phi -> phi o t_v（t_v は D_F の反射）で新しい両側剰余類を探し、sum 1/|Stab|（像の中で）が
|O(q)| / (|bar G_F| |bar G_K|) に達したら止める。そのとき、探索の方法によらず、全部の両側剰余類が見つかっている。
探索：近傍の全数試験（|Stab| の大きい類から。Stab の軌道ごとに代表 1 つ）、普通のランダムウォーク（歩き手 3 つ。
beta = 0, 1/2, 1 の Metropolis）、等変のランダムウォーク（固定部分群の素数位数の元と可換な反射で動く）を、
1 : 1 : 2 の時間で進める。乱数の種は 20261031 + (時刻 mod 100000) である。代表と見つかる順序は回ごとに
変わるが、類の数、|Stab|、質量は変わらない。
表：load_tables() が groups.json（produce/involutions/discriminant_action.py の出力）、必要なら
(12, 20) の組の表、|O(q)| の表（involution-scale.json）、成分の一覧を読む。load_tables() を呼んで
いなければ、part() と order_O() が既定のパスで読む。
パス：表は、既定では Release のアセットの computation-outputs/C4 と C7 から読む。search() は OUT に書く。
OUT の既定は今のディレクトリの下の results/C4 である。
"""
import hashlib
import json
import os
import random
import resource
import subprocess
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from sage.all import ZZ, matrix

from . import gluing as glue
from . import antiisometry_sampler as gs
from .io import asset

DATA = asset("computation-outputs")
PARTS_DIR = DATA / "C7"                                   # min4-twoel-*.json
GROUPS = DATA / "C4/groups.json"
GROUPS_12_20 = DATA / "C4/groups-12-20.json"             # optional / 無くてもよい
SCALE_FILE = DATA / "C4/involution-scale.json"
OUT = Path("results/C4")                                  # output of search() / search() の出力先

GR = None
SCALE = None
PAIRS = {"14-18b": ("14", "18b"), "14-18e": ("14", "18e"), "14-18a": ("14", "18a"), "14-18c": ("14", "18c"), "14-18d": ("14", "18d"),
         "16-16s14": ("16s14", "16s14"), "15-17": ("15", "17s15"), "16-16d8": ("16d8", "16d8"),
         "16-16-s16d0a": ("16c0", "16c0"), "16-16-s16d0b": ("16c0", "16c1"), "16-16-s16d0c": ("16c1", "16c1"), "16-16-s14d0": ("16c3", "16c3"),
         "16-16-s12d0": ("16c5", "16c5"), "16-16-s12d1": ("16c6", "16c6"), "16-16-s10d0": ("16c7", "16c7"), "16-16-s10d1": ("16c8", "16c8"), "16-16-s8d0": ("16c9", "16c9")}


def load_tables(groups=None, groups_12_20=None, scale=None, parts_dir=None):
    """Read the tables of the components; the pairs (12, 20) are added to PAIRS when their table exists.

    成分の表を読む。(12, 20) の組の表があれば、その組を PAIRS に足す。"""
    global GR, SCALE, PARTS_DIR
    if parts_dir is not None:
        PARTS_DIR = Path(parts_dir)
    GR = json.load(open(groups or GROUPS))
    g1220 = Path(groups_12_20 or GROUPS_12_20)
    if g1220.exists():
        GR.update(json.load(open(g1220)))
    sc = Path(scale or SCALE_FILE)
    SCALE = json.load(open(sc)) if sc.exists() else []
    PAIRS.update({f"12-20-c{v['cid']}": ("12", k) for k, v in GR.items() if "cid" in v})


def _tables():
    if GR is None:
        load_tables()


def part(name):
    """(Gram matrix, table of the discriminant form, entry of groups.json) of a component.

    成分の (Gram 行列、判別式形式の表、groups.json の項目)。"""
    _tables()
    g = GR[name]
    cl = json.load(open(PARTS_DIR / g["file"]))["classes"]
    G = matrix(ZZ, next(c for c in cl if c["cid"] == g["cid"])["gram"] if "cid" in g else cl[g["index"]]["gram"])
    return G, gs.disc(G), g


def order_O(gF, gK):
    """|O(q)|. For (16, 16) there are components with the same |Aut| and different (s, delta), so s and delta are matched as well.

    |O(q)|。(16,16) には |Aut| が同じで (s, delta) の違う成分があるので、s と delta でも照合する。"""
    _tables()
    if "O_q" in gF:
        return int(gF["O_q"])
    cF = json.load(open(PARTS_DIR / gF["file"]))["classes"][gF["index"]]
    for r in SCALE:
        if r["aut_F"] == str(gF["aut"]) and r["aut_K"] == str(gK["aut"]) and r["s"] == cF["s"] and r["delta"] == cF["delta"]:
            return int(r["O_q"])
    raise KeyError((gF["aut"], gK["aut"], cF["s"], cF["delta"]))


def first_phi(F, K, rng):
    """A first anti-isometry D_F -> D_K (random, retried at dead ends).

    最初の反等長 D_F -> D_K（乱択。行き止まりならやり直す）。"""
    for _ in range(100000):
        ys, _c = gs.random_antiisometry(F, K, rng)
        if ys is not None:
            return np.array(ys, dtype=np.int64)
    raise RuntimeError("no anti-isometry found")


def corehash():
    """Hash of the C part of rank32.gluing and of SMALL (the keys depend on both); used to decide whether a checkpoint can be resumed as it is.

    rank32.gluing の C の部分と SMALL のハッシュ（キーは両方に依る）。チェックポイントをそのまま再開できるかの判定に使う。"""
    return hashlib.sha256(Path(glue.__file__).with_name("_gluing.c").read_bytes() + str(glue.SMALL).encode()).hexdigest()


def record_run(pair, seed, limit, out):
    """Append a record of the start (arguments, seed, hash of the C part, commit of the code) to <pair>.runs.jsonl.

    起動の記録（引数、乱数の種、C の部分のハッシュ、コードのコミット）を <組>.runs.jsonl に追記する。"""
    code_dir = Path(__file__).resolve().parents[1]

    def git(*a):
        try:
            return subprocess.run(["git", "-C", str(code_dir), *a], capture_output=True, text=True).stdout.strip()
        except FileNotFoundError:   # git is not installed: the commit is recorded as unknown / git が無い：コミットは unknown と記録する
            return ""
    r = dict(pair=pair, started=time.strftime("%Y-%m-%d %H:%M:%S %z"), argv=sys.argv, limit_sec=limit, out=str(out), commit=git("rev-parse", "HEAD") or "unknown",
             dirty=git("status", "--porcelain", "--", ".").splitlines(), seed=seed, corehash=corehash(), job=os.environ.get("INVOCATION_ID", "interactive"))
    with open(Path(out) / f"{pair}.runs.jsonl", "a") as f:
        f.write(json.dumps(r, ensure_ascii=False) + "\n")


def wr(path, obj):
    """Write JSON through a temporary file and rename.

    一時ファイルに書いてから rename する。"""
    tmp = str(path) + ".tmp"
    json.dump(obj, open(tmp, "w"))
    os.replace(tmp, path)


def orbit_reps(g, els):
    """Representatives of the orbits of <els> on Q1F (phi o t_{av} and phi o t_v are in the same class, so only representatives are tried).

    Q1F の <els> の軌道の代表（phi o t_{av} と phi o t_v は同じ類なので代表だけ試す）。"""
    n = 1 << g.s
    label = np.arange(n)
    tabs = [glue.table(e) for e in els]
    while True:
        new = label.copy()
        for t in tabs:
            np.minimum.at(new, t, label)
            new = np.minimum(new, new[t])
        new = new[new]
        if (new == label).all():
            break
        label = new
    q = g.Q1F
    return np.unique(label[q])


def _save_checkpoint(path, classes, steps, expanded, elapsed):
    """The checkpoint as JSON (keys in hex), written through a temporary file.

    チェックポイントを JSON で（キーは 16 進）。一時ファイルを経由して書く。"""
    wr(path, dict(classes={k.hex(): c for k, c in classes.items()}, steps=steps, expanded=expanded, elapsed=elapsed))


def _load_checkpoint(path):
    d = json.load(open(path))
    return {bytes.fromhex(k): c for k, c in d["classes"].items()}, d["steps"], d["expanded"], d["elapsed"]


def search(pair, limit=600, out=None):
    """Find all double cosets of the pair by the search above, within limit seconds; resumable from <pair>.ckpt.json. Writes <pair>.json when the mass is complete.

    上の探索で組の両側剰余類を全部求める（制限 limit 秒）。<組>.ckpt.json から再開できる。質量が満ちたら <組>.json を書く。"""
    out = Path(out if out is not None else OUT)
    nF, nK = PAIRS[pair]
    GF_, F, gF = part(nF); GK_, K, gK = part(nK)
    t0 = time.time()
    g = glue.Glue(F, K, gF["gens_disc"], gK["gens_disc"], int(gF["image"]), int(gK["image"]))
    O = order_O(gF, gK); imF, imK = int(gF["image"]), int(gK["image"]); kF, kK = gF["kernel"], gK["kernel"]
    target = Fraction(O, imF * imK)
    print(f"{pair}: s={g.s}, stage sizes={[len(x) for x in g.stages]}, length of the chain of F={len(g.chainF)}, |O|={O}, target (mass of the image)={target} ~ {float(target):.4f}, M_up={float(target)/(kF*kK):.4f}", flush=True)
    seed = 20261031 + int(time.time()) % 100000          # a different random walk at every restart / 再開のたびに別の乱数で歩く
    rng = random.Random(seed); nprng = np.random.default_rng(seed)
    record_run(pair, seed, limit, out)
    ckpt = out / f"{pair}.ckpt.json"; prog = out / f"{pair}.progress.json"
    classes = {}; mass = Fraction(0); steps = 0; expanded = 0; prior = 0.0
    if ckpt.exists():
        classes, steps, expanded, prior = _load_checkpoint(ckpt)
        # the keys do not depend on the order of the search; recompute them only if the C part changed
        # キーは探索の順序に依らないので、C の部分を変えていなければそのまま使う
        if not (out / f"{pair}.corehash").exists() or (out / f"{pair}.corehash").read_text() != corehash():
            old = classes; classes = {}
            for c in old.values():
                k, n = g.canon(np.array(c["phi"], dtype=np.int64))
                assert n == c["stab"], ("|Stab| differs on resume", n, c["stab"])
                classes.setdefault(k, c)
        mass = sum((Fraction(1, c["stab"]) for c in classes.values()), Fraction(0))
        print(f"resume: classes {len(classes)}, mass {float(mass):.6f}", flush=True)
    last = 0.0; lastck = time.time() - 240; job = os.environ.get("INVOCATION_ID", "interactive")

    def report(state):
        ext = [c for c in classes.values() if c["ext"]]
        wr(prog, dict(job=job, pair=pair, stage="neighbour search", state=state, classes=len(classes), extremal_classes=len(ext), samples=steps, expanded=expanded,
                      unexpanded=sum(1 for c in classes.values() if not c["done"]), mass_float=float(mass), target_float=float(target), remaining_float=float(target - mass),
                      unit="double cosets (closed by the mass of the image)", updated=time.strftime("%Y-%m-%d %H:%M:%S %z"), elapsed_sec=prior + time.time() - t0, knodes=g.TK.nnodes, failed=0))

    stabbed = [k for k, c in classes.items() if c["stab"] >= 2]

    def add(phi):
        nonlocal mass, steps
        key, n = g.canon(phi); steps += 1
        c = classes.get(key)
        if c is None:
            classes[key] = dict(stab=n, phi=[int(x) for x in phi], ext=g.is_extremal(phi), hits=1, done=False)
            mass += Fraction(1, n)
            if n >= 2:
                stabbed.append(key)
        else:
            assert c["stab"] == n
            c["hits"] += 1
        return n

    def tick():
        nonlocal last, lastck
        if time.time() - last > 60:
            last = time.time(); report("running")
            print(f"  {prior + time.time()-t0:.0f} s: canon {steps} times, classes {len(classes)}, expanded {expanded}, mass {float(mass):.6f} / {float(target):.6f}, nodes K {g.TK.nnodes}, {steps/max(time.time()-t0,1):.0f} per second, memory {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss//1024} MB", flush=True)
        if time.time() - lastck > 300:
            lastck = time.time()
            _save_checkpoint(ckpt, classes, steps, expanded, prior + time.time() - t0)
            (out / f"{pair}.corehash").write_text(corehash())

    if not classes:
        add(first_phi(F, K, rng))
    t_bfs = t_walk = 0.0
    # three walkers: beta = 0 is the uniform distribution on phi (a class has probability ~ 1/|Stab|), beta = 1 is Metropolis
    # with the uniform distribution on classes, beta = 1/2 is in between. The proposal phi -> phi o t_v is accepted with
    # probability min(1, (|Stab'|/|Stab|)^beta). A new class is recorded even if it is not accepted (the mass does not depend on it).
    # ランダムウォークは 3 本：beta = 0 は phi の一様分布（類の確率は 1/|Stab| に比例）、beta = 1 は Metropolis で類の一様分布、beta = 1/2 はその中間。
    # 提案 phi -> phi o t_v を確率 min(1, (|Stab'|/|Stab|)^beta) で受理する。受理されなくても新しい類なら記録する（質量の和に確率は関わらない）。
    start = next(iter(classes.values()))
    walkers = [dict(phi=np.array(start["phi"], dtype=np.int64), n=start["stab"], beta=b_) for b_ in (0.0, 0.5, 1.0)]
    wi = 0

    def walk_step():
        nonlocal wi
        w = walkers[wi]; wi = (wi + 1) % len(walkers)
        cand = glue.compose(w["phi"], g.reflection(int(g.Q1F[nprng.integers(len(g.Q1F))])))
        n2 = add(cand)
        if w["beta"] == 0 or n2 >= w["n"] or rng.random() < (n2 / w["n"]) ** w["beta"]:
            w["phi"] = cand; w["n"] = n2

    # equivariant walk: take an element (a, b) of prime order of the Stab of a known class and move only by products of reflections
    # commuting with a. It walks among the phi fixed by (a, b), so it meets classes with a large |Stab| more often.
    # 等変のランダムウォーク：既知の類の Stab の素数位数の元 (a, b) を選び、a と可換な反射の積だけで動く。(a, b) を保つ phi の中を歩くので、|Stab| の大きい類に当たりやすい。
    chain = None

    def equiv_step():
        nonlocal chain
        if chain is None or chain["left"] == 0:
            phi0 = np.array(classes[stabbed[rng.randrange(len(stabbed))]]["phi"], dtype=np.int64)
            chain = dict(phi=phi0, a=g.stab_prime_element(phi0, rng), left=32)
        if chain["a"] is None:
            chain = None
            return
        chain["left"] -= 1
        for _ in range(64):
            tv = g.orbit_reflection(chain["a"], int(g.Q1F[nprng.integers(len(g.Q1F))]))
            if tv is not None:
                chain["phi"] = glue.compose(chain["phi"], tv); add(chain["phi"])
                return

    # exhaustive test of the neighbours (classes with a large |Stab| first, representatives of the orbits of Stab only) :
    # ordinary walk (meets a class with probability ~ 1/|Stab|) : equivariant walk = 1 : 1 : 2 in time
    # 近傍の全数試験（|Stab| の大きい類から、Stab の軌道の代表だけ）: 普通のランダムウォーク（出会う確率は 1/|Stab| に比例）: 等変のランダムウォーク = 1 : 1 : 2 の時間で進める
    t_eq = 0.0
    while mass < target and time.time() - t0 < limit:
        todo = [(c["stab"], k) for k, c in classes.items() if not c["done"]]
        if not todo:
            break
        k = max(todo, key=lambda z: z[0])[1]
        c = classes[k]
        phi = np.array(c["phi"], dtype=np.int64)
        _, n, els = g.canon(phi, want_stab=True)
        reps = orbit_reps(g, els)
        finished = True
        for v in reps:
            if mass >= target or time.time() - t0 >= limit:
                finished = False; break
            t1 = time.time(); add(glue.compose(phi, g.reflection(int(v)))); t_bfs += time.time() - t1
            while t_walk < t_bfs:
                t1 = time.time(); walk_step(); t_walk += time.time() - t1
            while stabbed and mass < target and t_eq < t_bfs + t_walk:
                t1 = time.time(); equiv_step(); t_eq += time.time() - t1
            tick()
        if finished:
            c["done"] = True; expanded += 1
    report("done" if mass == target else "stopped")
    _save_checkpoint(ckpt, classes, steps, expanded, prior + time.time() - t0)
    assert mass <= target, ("the mass exceeds the target", mass, target)
    print(f"canon {steps} times, classes {len(classes)}, mass {mass} {'= target' if mass == target else '!= target ' + str(target)}, {time.time()-t0:.0f} s, nodes K {g.TK.nnodes}", flush=True)
    print("distribution of |Stab|", sorted(Counter(c["stab"] for c in classes.values()).items()))
    ext = [c for c in classes.values() if c["ext"]]
    mext = sum((Fraction(1, c["stab"]) for c in ext), Fraction(0))
    print(f"extremal classes {len(ext)}, extremal mass of the image {mext} ~ {float(mext):.6f}, p_phi = {float(mext/target):.6f}, M_ext = {mext/(kF*kK)} ~ {float(mext)/(kF*kK):.6f}")
    if mass == target:
        rows = [dict(phi=c["phi"], stab_image=c["stab"], stab=c["stab"] * kF * kK, extremal=c["ext"], hits=c["hits"]) for c in classes.values()]
        wr(out / f"{pair}.json", dict(pair=pair, F=[gF["file"], gF["index"]], K=[gK["file"], gK["index"]], K_cid=gK.get("cid"), s=g.s, O_q=str(O), aut_F=gF["aut"], aut_K=gK["aut"], kernels=[kF, kK],
                                      target_image=str(target), M_up=str(target / (kF * kK)), sum_all=str(mass / (kF * kK)), M_ext=str(mext / (kF * kK)), n_classes=len(rows), n_extremal=len(ext), closed=True, classes=rows))
