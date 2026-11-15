"""Tests of rank32.antiisometry_sampler and rank32.discriminant_action on sqrt(2) E_8.

D = (sqrt(2) E_8)^# / sqrt(2) E_8 is F_2^8 with the quadratic form of plus type; the classes of norm
1 (2q = 2 mod 4) are the 120 of the roots / 2. Needs SageMath (PARI). Time: a few seconds.

rank32.antiisometry_sampler と rank32.discriminant_action の試験（sqrt(2) E_8 で）。

D = (sqrt(2) E_8)^# / sqrt(2) E_8 は + 型の 2 次形式をもつ F_2^8 で、ノルム 1（2q = 2 mod 4）の類はルートの半分の
120 個である。SageMath（PARI）が要る。時間：数秒。
"""
import random
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

pytest.importorskip("sage.all")

from sage.all import ZZ, CartanMatrix, matrix  # noqa: E402

from rank32 import antiisometry_sampler as gs  # noqa: E402
from rank32 import discriminant_action as da  # noqa: E402
from rank32 import forms  # noqa: E402

G = 2 * matrix(ZZ, CartanMatrix(["E", 8]))


def test_disc_table():
    t = gs.disc(G)
    assert t["k"] == 8 and len(t["q2"]) == 256
    assert int((t["q2"] == 2).sum()) == 120 and int((t["q2"] == 0).sum()) == 136
    assert t["mu2"][0] == 0 and set(t["mu2"][1:].tolist()) <= {2, 4, 8}
    assert ((t["mu2"] == 2) == (t["q2"] == 2)).all()           # norm 1 classes contain the roots / 2 / ノルム 1 の類はルートの半分を含む


def test_random_antiisometry_is_antiisometry():
    t = gs.disc(G)
    found = 0
    for seed in range(10):
        ys, counts = gs.random_antiisometry(t, t, random.Random(seed))
        assert len(counts) >= 1
        if ys is None:
            continue
        found += 1
        phi = gs.phi_table(ys, 8)
        assert len(set(phi.tolist())) == 256
        assert ((t["q2"][phi] + t["q2"]) % 4 == 0).all()
    assert found


def test_action_preserves_the_form():
    order, gens = forms.qfauto(G)
    assert order == 696729600
    t = gs.disc(G)
    b = da.disc_basis(G)
    assert len(b[0]) == 8
    for g in gens:
        im = da.action(G, g, *b)
        tab = np.zeros(256, dtype=np.int64)
        for l, y in enumerate(im):
            tab[1 << l: 2 << l] = tab[0: 1 << l] ^ y
        assert len(set(tab.tolist())) == 256 and (t["q2"][tab] == t["q2"]).all()
