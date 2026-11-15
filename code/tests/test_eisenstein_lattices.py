"""Tests of rank32.eisenstein_overlattices, rank32.orbit_walk and rank32.prime_order_gluing on small lattices.

E8 is built twice: as an overlattice of A2^4 (with its Eisenstein structure) and as A2 glued to E6.
Time: about 20 seconds.

rank32.eisenstein_overlattices、rank32.orbit_walk、rank32.prime_order_gluing の、小さい格子での試験。

E8 を 2 通りに作る：A2^4 の上格子として（Eisenstein 構造つきで）と、A2 と E6 の貼り合わせとして。
時間：約 20 秒。
"""
import random
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, CartanMatrix, GF, block_diagonal_matrix, matrix, pari  # noqa: E402

from rank32 import eisenstein_overlattices as ov  # noqa: E402
from rank32 import orbit_walk as w  # noqa: E402
from rank32 import prime_order_gluing as pg  # noqa: E402

A2 = matrix(ZZ, [[2, -1], [-1, 2]])
W2 = matrix(ZZ, [[1, 1], [-1, 0]])


def eisenstein_e8():
    G, W = block_diagonal_matrix([A2] * 4, subdivide=False), block_diagonal_matrix([W2] * 4, subdivide=False)
    for _ in range(2):
        D, iso, _good = ov.good_lines(G)
        G, W = ov.overlattice(G, W, D, sorted(iso)[0])
    return G, W


def test_overlattices_give_e8():
    G, W = eisenstein_e8()
    assert G.det() == 1 and W * G * W.transpose() == G
    assert w.counts(G) == (240, 2400, 9120)   # 240, 2160, 6720 vectors of norm 2, 4, 6 / ノルム 2, 4, 6 のベクトルが 240, 2160, 6720 本


def test_lines_and_orbits():
    G, W = block_diagonal_matrix([A2] * 4, subdivide=False), block_diagonal_matrix([W2] * 4, subdivide=False)
    D, iso, good = ov.good_lines(G)
    assert D.s == 4 and len(iso) == 16 and not good   # 16 of the 40 lines are isotropic, and each meets a root / 直線 40 本のうち等方は 16 本で、どれもルートに当たる
    order, gens = ov.qfsafe.qfauto(ov.forms(G, W), tries=30)
    assert order == 6 ** 4 * 24
    reps = ov.orbits(D, gens, iso)
    assert sum(size for _, size in reps) == len(iso)
    assert all(ov.line_key(ov.key_to_vec(k, 4)) == k for k in iso)


def test_neighbours_of_e8():
    G, W = eisenstein_e8()
    order, gens = w.qfsafe.qfauto(w.forms(G, W), tries=20)
    reps, total = w.isotropic_reps(G, gens + [W.transpose()], random.Random(1), 0)
    assert total == len(reps) and total >= 1
    for x in reps:
        for Gn, Wn in w.herm_neighbours(G, W, w.to_vec(x, 8)):
            assert Gn.det() == 1 and Wn * Gn * Wn.transpose() == Gn
    E8 = matrix(ZZ, CartanMatrix(["E", 8]))
    D16 = block_diagonal_matrix([E8, E8], subdivide=False)
    v = w.to_vec(1 | 1 << 8, 16)
    nb = w.real_neighbour(D16, v)
    assert all(Gn.det() == 1 for Gn, _ in nb)


def test_real_seeds_and_short_basis():
    seeds = w.real_seeds(3, 12, 6)
    assert seeds and all(ZZ(G.det()) == 3 ** 6 and G.nrows() == 12 for G in seeds)
    G = 3 * matrix(ZZ, CartanMatrix(["E", 8]))
    G2, _ = w.short_basis(G, None, random.Random(2), bound=6)
    assert G2.det() == G.det()


def test_gluing_a2_e6_gives_e8():
    F, Z = A2, matrix(ZZ, CartanMatrix(["E", 6]))
    DF, DZ = pg.Disc(F, 3), pg.Disc(Z, 3)
    assert DF.s == DZ.s == 1
    assert len(DF.norm_classes(2)) == 2 and len(DZ.norm_classes(4)) == 2   # norm 2/3 in A2^#, norm 4/3 in E6^# / A2^# のノルム 2/3、E6^# のノルム 4/3
    T = pg.random_phi(DF, DZ, random.Random(0))
    L = pg.build_L(F, Z, DF, DZ, T)
    assert L.det() == 1 and int(pari(L).qfminim(2, 0)[0]) == 240


def test_nullspace():
    rng = np.random.default_rng(0)
    for _ in range(20):
        M = rng.integers(0, 3, size=(3, 6))
        B = pg.nullspace(M, 3)
        assert len(B) == 6 - matrix(GF(3), M.tolist()).rank()
        assert all(not ((M @ b) % 3).any() for b in B)
