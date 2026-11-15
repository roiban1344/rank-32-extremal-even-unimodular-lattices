"""Tests of rank32.gaussian_keys: invariants, canonical forms, |U| and a step of the pi-neighbour walk.

Two rootless lattices of the walk (tests/data/gaussian_lattices.json): one with a canonical form
(status discrete, |U| = 4) and one whose unique-colour lines do not generate L (|U| = 8).
Needs SageMath and gcc. Time: a few seconds.

rank32.gaussian_keys の試験：不変量、標準形、|U|、pi 近傍のランダムウォークの 1 歩。

ランダムウォークのルート無しの格子 2 つ（tests/data/gaussian_lattices.json）：標準形を持つもの（状態 discrete、|U| = 4）と、
色が唯一の直線が L を生成しないもの（|U| = 8）。
SageMath と gcc が要る。時間：数秒。
"""
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, matrix, identity_matrix  # noqa: E402

from rank32 import gaussian_keys as gk  # noqa: E402

DATA = json.loads((Path(__file__).with_name("data") / "gaussian_lattices.json").read_text())["lattices"]


@pytest.mark.parametrize("k", range(len(DATA)))
def test_analyse_and_unitary_order(k):
    r = DATA[k]
    a = gk.analyse(r["gram"], r["W"])
    assert a["roots"] == 0 and (a["A"], a["C"], a["D"], a["status"]) == (r["A"], r["C"], r["D"], r["status"])
    assert a.get("key") == r["key"]
    assert gk.unitary_order(r["gram"], r["W"])["U"] == r["U"]


def test_iso_and_step():
    a, b = DATA
    assert gk.iso(a, a)["linear"]
    assert gk.iso(a, b) == {"linear": False, "antilinear": False}
    G, W = matrix(ZZ, a["gram"]), matrix(ZZ, a["W"])
    G1, W1 = gk.step(G, W, random.Random(1))
    G2, W2 = gk.step(G, W, random.Random(1))
    assert (G1, W1) == (G2, W2)
    assert G1.det() == 1 and all(G1[i, i] % 2 == 0 for i in range(32))
    assert W1 * W1 == -identity_matrix(ZZ, 32) and W1 * G1 * W1.transpose() == G1


def test_enum_short_counts():
    import numpy as np
    V = gk.enum_short(np.array(DATA[0]["gram"]))
    G = np.array(DATA[0]["gram"])
    norms = np.einsum("ij,ij->i", V @ G, V)
    assert len(V) == 73440 and (norms == 4).all()      # rootless extremal: 146880 minimal vectors / ルート無しの極値格子
