"""Tests of rank32.canonical: known group orders and isometry tests (the bliss route, independent of PARI).

Needs gcc and bliss (the C program is compiled on first use). Time: a few seconds.

rank32.canonical の試験：既知の群の位数と同型判定（PARI と独立な bliss の経路）。

gcc と bliss が要る（C のプログラムを初めて使うときにコンパイルする）。時間：数秒。
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, matrix, CartanMatrix, block_diagonal_matrix  # noqa: E402

from rank32 import canonical  # noqa: E402

A2 = matrix(ZZ, [[2, -1], [-1, 2]])
W2 = matrix(ZZ, [[1, 1], [-1, 0]])


def test_weyl_orders():
    assert canonical.canonical(matrix(ZZ, CartanMatrix(["E", 8])))["order"] == 696729600
    assert canonical.canonical(matrix(ZZ, CartanMatrix(["E", 6])))["order"] == 103680
    assert canonical.canonical(matrix(ZZ, CartanMatrix(["D", 4])))["order"] == 1152


def test_hermitian_orders():
    """A_2^k with the Eisenstein structure omega: |U| = 6^k k!; without omega: 12^k k!.

    Eisenstein の構造 omega つきの A_2^k は |U| = 6^k k!、omega なしでは 12^k k!。"""
    for k, order in ((1, 6), (2, 72), (3, 1296)):
        G = block_diagonal_matrix([A2] * k, subdivide=False)
        W = block_diagonal_matrix([W2] * k, subdivide=False)
        assert canonical.canonical(G, W)["order"] == order
        assert canonical.canonical(G)["order"] == 12 ** k * ZZ(k).factorial()


def test_isometry_changes_basis():
    G = matrix(ZZ, CartanMatrix(["E", 6]))
    T = matrix(ZZ, 6, 6, lambda i, j: 1 if i == j else (1 if j == i + 1 else 0))
    H = T * G * T.transpose()
    a, b = canonical.canonical(G), canonical.canonical(H)
    R = canonical.isometry(a, G, b, H)
    assert R is not None and R * H * R.transpose() == G


def test_non_isometric_same_theta_prefix():
    """E8 and D4 + D4 (both of rank 8) are found not isometric.

    E8 と D4 + D4（どちらも階数 8）が、同型でないと判定される。"""
    G = matrix(ZZ, CartanMatrix(["E", 8]))
    G2 = block_diagonal_matrix([matrix(ZZ, CartanMatrix(["D", 4]))] * 2, subdivide=False)
    a, b = canonical.canonical(G), canonical.canonical(G2)
    assert canonical.isometry(a, G, b, G2) is None


def test_shell_hash_invariant():
    G = matrix(ZZ, CartanMatrix(["E", 6]))
    T = matrix(ZZ, 6, 6, lambda i, j: 1 if i == j else (1 if j == i + 1 else 0))
    H = T * G * T.transpose()
    assert canonical.shell_hash(G) == canonical.shell_hash(H)
    A = matrix(ZZ, CartanMatrix(["D", 4]))
    B = block_diagonal_matrix([matrix(ZZ, CartanMatrix(["A", 1]))] * 4, subdivide=False)
    assert canonical.shell_hash(A) != canonical.shell_hash(B)


def test_generators_convention():
    G = matrix(ZZ, CartanMatrix(["E", 6]))
    r = canonical.canonical(G)
    for g in canonical.generators(r, G):
        assert g.transpose() * G * g == G
