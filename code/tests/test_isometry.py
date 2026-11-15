"""Tests of rank32.isometry on small lattices.

rank32.isometry の小さい格子での試験。
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, matrix, identity_matrix, CartanMatrix  # noqa: E402

from rank32.isometry import is_isometry, is_automorphism, all_automorphisms  # noqa: E402

D4 = matrix(ZZ, CartanMatrix(["D", 4]))


def test_change_of_basis():
    U = identity_matrix(ZZ, 4)
    U[0, 1] = 1
    assert is_isometry(U, D4, U.T * D4 * U)
    assert not is_isometry(U, D4, D4)


def test_determinant_is_checked():
    T = 2 * identity_matrix(ZZ, 4)
    assert T.T * D4 * T == 4 * D4 and not is_isometry(T, D4, 4 * D4)


def test_reflections_are_automorphisms():
    refl = []
    for i in range(4):
        S = identity_matrix(ZZ, 4)
        for j in range(4):
            S[i, j] -= D4[i, j]          # column j: e_j - C_ij e_i / 列 j：e_j - C_ij e_i
        refl.append(S)
    assert all_automorphisms(refl + [-identity_matrix(ZZ, 4)], D4)
    P = matrix(ZZ, 4, 4, lambda i, j: 1 if j == (i + 1) % 4 else 0)
    assert not is_automorphism(P, D4)
