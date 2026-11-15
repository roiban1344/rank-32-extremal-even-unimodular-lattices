"""Tests of rank32.permutation_model on (2E_8)^4 (960 vectors of norm 4, Aut = W(E_8) wr S_4).

The group order is |W(E_8)|^4 * 24, and perm_to_row turns the permutations of the generators
back into their row actions. Needs SageMath (with GAP). Time: a few seconds.

rank32.permutation_model の (2E_8)^4（ノルム 4 のベクトル 960 本、Aut = W(E_8) wr S_4）での試験。

群の位数が |W(E_8)|^4 * 24 であること、perm_to_row が生成元の置換を行の作用に戻すことを確かめる。
SageMath（GAP を含む）が要る。時間：数秒。
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import CartanMatrix, ZZ, block_diagonal_matrix, identity_matrix, matrix, libgap  # noqa: E402

from rank32.permutation_model import perm_model, perm_to_row  # noqa: E402


def e8_4():
    """(2E_8)^4 and generators of its group (simple reflections, a 4-cycle and a transposition of the components); g^T G g = G.

    (2E_8)^4 とその群の生成元（各成分の単純鏡映、成分の巡回と互換）。g^T G g = G。"""
    C = matrix(ZZ, CartanMatrix(["E", 8]))
    refl = []
    for i in range(8):
        S = identity_matrix(ZZ, 8)
        for j in range(8):
            S[i, j] -= C[i, j]          # column j: e_j - C_ij e_i / 列 j：e_j - C_ij e_i
        refl.append(S)
    G = block_diagonal_matrix([2 * C] * 4)
    I8 = identity_matrix(ZZ, 8)
    Z = matrix(ZZ, 8, 8)
    gens = [block_diagonal_matrix([S if c == b else I8 for c in range(4)]) for b in range(4) for S in refl]
    def perm_blocks(p):
        return matrix.block([[I8 if p[c] == r else Z for c in range(4)] for r in range(4)])
    gens += [perm_blocks([1, 2, 3, 0]), perm_blocks([1, 0, 2, 3])]
    assert all(g.transpose() * G * g == G for g in gens)
    return G, gens


def test_order_and_rows():
    G, gens = e8_4()
    K, pts = perm_model(G, gens, "points")
    assert len(pts) == 960 and int(libgap.Size(K)) == 696729600 ** 4 * 24
    K2, xs = perm_model(G, gens, [g.transpose() for g in gens[-3:]])
    assert int(libgap.Size(K2)) == int(libgap.Size(K))
    for x, g in zip(xs, gens[-3:]):
        assert perm_to_row(x, pts) == g.transpose()


def test_minus_one():
    G, _ = e8_4()
    K, pts = perm_model(G, [-identity_matrix(ZZ, 32)], "points")
    assert int(libgap.Size(K)) == 2
