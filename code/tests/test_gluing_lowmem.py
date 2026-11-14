"""Tests of rank32.gluing_lowmem: building the chain on the F side in slices gives the same chain and the same keys as rank32.gluing.

Input: data/gluing_pairs.json (see test_gluing.py). Needs gcc and GAP (through Sage).
Time: about 30 seconds.

rank32.gluing_lowmem の試験：F 側の鎖を区切って作っても、鎖と正規形のキーが rank32.gluing と同じである。

入力：data/gluing_pairs.json（test_gluing.py を参照）。gcc と GAP（Sage 経由）が要る。
時間：約 30 秒。
"""
import hashlib
import random
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする
sys.path.insert(0, str(Path(__file__).resolve().parent))

pytest.importorskip("sage.all")

from rank32 import gluing, gluing_lowmem  # noqa: E402
from test_gluing import PAIRS, build  # noqa: E402


def test_chunked_orbits_equal():
    """The sliced orbits_of returns the same arrays as orbits_of (on nodes that keep all elements, with small slices).

    区切った orbits_of が元の orbits_of と同じ配列を返す（元を全部持つノードで、区切りを小さくして）。"""
    g, _ = build("16-16-s12d0", gluing.Glue)
    node = g.TF.root
    pts = np.arange(1 << g.s)
    seen_elems = False
    while node.order > 1:
        a = node.orbits_of(pts, 0)
        b = gluing_lowmem.orbits_of_chunked(node, pts, 0, chunk_bytes=1 << 16)
        assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])
        seen_elems |= node.elems is not None
        l, sz = a
        reps = pts[(l == pts) & (sz > 1)]
        node = node.child(int(reps[np.argmin(sz[reps])]), 0)
    assert seen_elems, "no node keeping all elements was reached (the test checked nothing)"


@pytest.mark.parametrize("name", sorted(PAIRS))
def test_same_chain_and_keys_as_gluing(name):
    """The subclass has the same chain and the same keys as rank32.gluing.Glue.

    部分クラスの鎖とキーが rank32.gluing.Glue と同じ。"""
    g0, _ = build(name, gluing.Glue)
    g1, _ = build(name, gluing_lowmem.Glue)
    assert len(g0.chainF) == len(g1.chainF)
    for (o0, u0), (o1, u1) in zip(g0.chainF, g1.chainF):
        assert np.array_equal(o0, o1) and np.array_equal(u0, u1)
    rng = random.Random(31)
    classes = PAIRS[name]["classes"]
    for c in rng.sample(classes, min(25, len(classes))):
        phi = np.array(c["phi"], dtype=np.int64)
        assert g0.canon(phi) == g1.canon(phi)


@pytest.mark.parametrize("name", sorted(PAIRS))
def test_keys_match_list(name):
    """All keys of the pair, sorted and concatenated, have the sha256 of the list, and |Stab| agrees.

    組の全部の類のキーをソートして連結した sha256 が一覧のものと同じで、|Stab| も一致する。"""
    g, kFK = build(name, gluing_lowmem.Glue)
    p = PAIRS[name]
    keys = []
    for c in p["classes"]:
        k, n = g.canon(np.array(c["phi"], dtype=np.int64))
        assert n * kFK == c["stab"]
        keys.append(k)
    assert len(set(keys)) == len(keys)
    assert hashlib.sha256(b"".join(sorted(keys))).hexdigest() == p["keys_sha256"]
