"""Tests of rank32.gluing: the normal form of the double cosets, and the mass of two closed pairs.

Input: data/gluing_pairs.json, two pairs of 2-elementary parts (F, K) of the involution gluing
with the tables of their discriminant forms, generators of bar G_F and bar G_K on the discriminant
groups, and the full list of double cosets (phi, |Stab|, extremal) of each pair:
  16-16-s12d0  involution type (16, 16) of §4.1 of the first paper, the components with s = 12, delta = 0;
               F = K, 122 classes;
  12-20-c14    involution type (12, 20) of §4.1: F = sqrt(2) D_12^+ and K one of the 66 lattices of rank 20 with
               s = 12, delta = 1 and minimum at least 4, the one with index 14 (counting from 0) in
               certificates/involution-gluing-masses.json; F != K, 1605 classes.
The names are the identifiers of the pairs in the certificates (certificates/README.md, "Identifiers").
keys_sha256 is the sha256 of the sorted, concatenated keys canon(phi) of all classes.
Needs gcc and GAP (through Sage). Time: about 30 seconds.

rank32.gluing の試験：両側剰余類の正規形と、閉じた 2 組の質量。

入力：data/gluing_pairs.json。対合の貼り合わせの 2-初等的な成分の組 (F, K) を 2 つ収める。
判別式形式の表、判別式群の上の bar G_F と bar G_K の生成元、各組の両側剰余類の全部の一覧
（phi、|Stab|、極値か）である。
  16-16-s12d0  第 1 の論文の §4.1 の対合の型 (16, 16)。s = 12、delta = 0 の成分。F = K、122 類
  12-20-c14    §4.1 の対合の型 (12, 20)。F = sqrt(2) D_12^+、K は s = 12、delta = 1、最小値 4 以上の階数 20 の格子 66 個の
               うち、certificates/involution-gluing-masses.json で（0 から数えて）14 番の格子。F != K、1605 類
名前は certificates の組の識別子である（certificates/README.md の「Identifiers」）。
keys_sha256 は、全部の類のキー canon(phi) をソートして連結したものの sha256 である。
gcc と GAP（Sage 経由）が要る。時間：約 30 秒。
"""
import hashlib
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

pytest.importorskip("sage.all")

from rank32 import gluing  # noqa: E402

DATA = Path(__file__).resolve().parent / "data/gluing_pairs.json"


def load_pairs():
    """The pairs of the input file, with the tables as numpy arrays.

    入力のファイルの組。表は numpy の配列にする。"""
    pairs = json.loads(DATA.read_text())
    for p in pairs.values():
        for side in ("F", "K"):
            d = p[side]
            k = d["k"]
            masks = np.arange(1 << k, dtype=np.int64)
            d["q2"] = np.array([int(c) for c in d["q2"]], dtype=np.int64)
            d["mu2"] = np.array([int(c) for c in d["mu2"]], dtype=np.int64)
            d["B4"] = np.array(d["B4"], dtype=np.int64)
            d["bits"] = ((masks[:, None] >> np.arange(k)) & 1).astype(np.int64)
    return pairs


PAIRS = load_pairs()


def build(name, cls=None):
    """(Glue of the pair, |kernel of Aut F -> O(q_F)| * |kernel of Aut K -> O(q_K)|).

    (組の Glue, |Aut F -> O(q_F) の核| * |Aut K -> O(q_K) の核|)。"""
    p = PAIRS[name]
    F, K = p["F"], p["K"]
    cls = cls or gluing.Glue
    return cls(F, K, F["gens"], K["gens"], int(F["image"]), int(K["image"])), F["kernel"] * K["kernel"]


@pytest.mark.parametrize("name", sorted(PAIRS))
def test_mass_closes(name):
    """The list closes: sum of 1/|Stab| equals |O(q)|/(|Aut F| |Aut K|), and the extremal part equals M_ext.

    一覧が閉じている：1/|Stab| の和が |O(q)|/(|Aut F| |Aut K|) に、極値の部分が M_ext に等しい。"""
    p = PAIRS[name]
    total = sum(Fraction(1, c["stab"]) for c in p["classes"])
    assert total == Fraction(int(p["O_q"]), int(p["F"]["aut"]) * int(p["K"]["aut"])) == Fraction(p["M_up"])
    assert sum(Fraction(1, c["stab"]) for c in p["classes"] if c["extremal"]) == Fraction(p["M_ext"])


@pytest.mark.parametrize("name", sorted(PAIRS))
def test_keys_match_list(name):
    """The keys of all classes are distinct, |Stab| agrees with the list, and the sha256 of the keys agrees.

    全部の類のキーが互いに異なり、|Stab| が一覧と一致し、キーの sha256 も一致する。"""
    g, kFK = build(name)
    p = PAIRS[name]
    keys = []
    for c in p["classes"]:
        k, n = g.canon(np.array(c["phi"], dtype=np.int64))
        assert n * kFK == c["stab"] and n == c["stab_image"]
        keys.append(k)
    assert len(set(keys)) == len(keys)
    assert hashlib.sha256(b"".join(sorted(keys))).hexdigest() == p["keys_sha256"]


def test_canon_invariant():
    """The normal form does not change under bar G_F x bar G_K, and extremality is kept.

    正規形は bar G_F x bar G_K の元で変わらず、極値かどうかも保たれる。"""
    name = "12-20-c14"
    g, _ = build(name)
    p = PAIRS[name]
    rng = random.Random(4)
    gensF = [np.array(x, dtype=np.int64) for x in p["F"]["gens"]]
    gensK = [np.array(x, dtype=np.int64) for x in p["K"]["gens"]]

    def word(gens):
        A = 1 << np.arange(g.s, dtype=np.int64)
        for _ in range(25):
            A = gluing.compose(gens[rng.randrange(len(gens))], A)
        return A

    for c in rng.sample(p["classes"], 40):
        phi = np.array(c["phi"], dtype=np.int64)
        moved = gluing.compose(word(gensK), gluing.compose(phi, gluing.inverse(word(gensF))))
        assert g.is_antiisometry(moved)
        assert g.canon(moved) == g.canon(phi)
        assert g.is_extremal(moved) == c["extremal"]


def test_orbit_reflection_keeps_stabilizer():
    """Equivariant move: moving by an isometry that commutes with a keeps (a, b) in Stab (the reason for the move that finds classes with large |Stab|).

    等変の動き：a と可換な等長写像で動かすと (a, b) は Stab に残る（|Stab| の大きい類を探す動きの根拠）。"""
    name = "16-16-s12d0"
    g, _ = build(name)
    p = PAIRS[name]
    rng = random.Random(3); nprng = np.random.default_rng(3); tried = 0
    cand = [c for c in p["classes"] if c["stab_image"] >= 2]
    for c in rng.sample(cand, min(12, len(cand))):
        phi = np.array(c["phi"], dtype=np.int64)
        a = g.stab_prime_element(phi, rng)
        b = gluing.compose(phi, gluing.compose(a, gluing.inverse(phi)))
        assert np.array_equal(gluing.compose(b, phi), gluing.compose(phi, a))
        for _ in range(64):
            t = g.orbit_reflection(a, int(g.Q1F[nprng.integers(len(g.Q1F))]))
            if t is None:
                continue
            assert np.array_equal(gluing.compose(a, t), gluing.compose(t, a))
            phi2 = gluing.compose(phi, t)
            assert g.is_antiisometry(phi2)
            assert np.array_equal(gluing.compose(b, phi2), gluing.compose(phi2, a))
            _, n, els = g.canon(phi2, want_stab=True)
            assert n >= 2
            tried += 1
            break
    assert tried >= 6


def test_matrix_helpers():
    """The C helpers: apply agrees with the naive computation for 0-, 1- and higher-dimensional input, and inverse is the inverse matrix.

    C のヘルパー：apply は 0 次元・1 次元・多次元で素朴な計算と同じ、inverse は逆行列。"""
    rng = np.random.default_rng(1)
    s = 12
    while True:
        A = rng.integers(0, 1 << s, s, dtype=np.int64)
        try:
            R = gluing.inverse(A)
            break
        except AssertionError:
            continue
    naive = lambda x: int(np.bitwise_xor.reduce([A[i] for i in range(s) if (int(x) >> i) & 1], initial=0))  # noqa: E731
    m = rng.integers(0, 1 << s, 40, dtype=np.int64)
    assert np.array_equal(gluing.apply(A, m), [naive(x) for x in m])
    assert np.array_equal(gluing.apply(A, m[:5]), [naive(x) for x in m[:5]])
    assert np.array_equal(gluing.apply(A, m.reshape(8, 5)).ravel(), [naive(x) for x in m])
    assert int(gluing.apply(A, 5)) == naive(5) and gluing.apply(A, 5).shape == ()
    assert np.array_equal(gluing.compose(A, R), 1 << np.arange(s))


def test_verify_order():
    """The order check in C: 1 for true generators, never 1 when twice the order is asked; _enumerate lists exactly order elements.

    C の位数の確認：本物の生成元では 1、2 倍の位数を要求しても 1 を返さない。_enumerate は位数どおり。"""
    g, _ = build("16-16-s12d0")
    T = g.TK
    node = T.root
    seen = 0
    for _ in range(4):
        label, size, *_ = node.tables(0)
        reps = np.nonzero((label == np.arange(len(label))) & (size > 1))[0]
        rep = int(T.domains[0][reps[0]])
        ch = node.child(rep, 0)
        if ch.order <= 1 or ch.elems is not None:
            break
        G = np.ascontiguousarray(np.stack(ch.gens), dtype=np.int64)
        L = gluing.core()
        assert L.verify_order(G.ctypes.data, len(G), g.s, ch.order, 7, 40) == 1
        assert L.verify_order(G.ctypes.data, len(G), g.s, 2 * ch.order, 7, 40) == 0
        node = ch
        seen += 1
    assert seen >= 1, "no large node was reached (the test checked nothing)"
    H = T._enumerate([np.array([2, 1] + [1 << i for i in range(2, g.s)], dtype=np.int64)], 2)
    assert H.shape == (2, g.s)
