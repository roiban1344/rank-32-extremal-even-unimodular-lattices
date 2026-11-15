"""Tests of rank32.theta_lp_basis and rank32.joint_theta_lp: theta series of the blocks, the spaces, and one exclusion.

Time: about 30 seconds.

rank32.theta_lp_basis と rank32.joint_theta_lp の試験：ブロックのテータ級数、空間、除外の 1 例。

時間：約 30 秒。
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")

from rank32 import joint_theta_lp as lp  # noqa: E402
from rank32 import theta_lp_basis as tb  # noqa: E402


def thetas(p):
    return {name: (tb.theta_q(G, 1), tb.theta_q(Gd, 1), n, s) for name, G, Gd, n, s in tb.blocks(p)}


def test_theta_e8():
    th = thetas(3)
    f = th["E8"][0]
    assert [f[j] for j in range(4)] == [1, 240, 2160, 6720]
    assert [f[j] for j in range(3)] == [th["rt3E8"][1][j] for j in range(3)]   # sqrt(3) (sqrt(3) E8)^# = E8
    assert [th["A2"][0][j] for j in range(3)] == [1, 6, 0]


def test_spaces_span():
    th = thetas(3)
    pairs, rank, dim = tb.spaces_for(3, 12, 6, th)
    assert rank == dim and tb.consistent(pairs)


def test_types_and_counts():
    # the types of p = 3, 5, 7 in rank 32 that pass the genus conditions and the density sieve (17, 6 and 2)
    # 種の条件と密度の篩を通る階数 32 の p = 3, 5, 7 の型（17、6、2 個）
    assert [len(lp.types_for(p)) for p in (3, 5, 7)] == [17, 6, 2]
    for p, d, s in ((3, 12, 6), (5, 16, 4), (7, 14, 3)):
        assert sum(lp.class_counts(p, d, s).values()) == p ** s - 1
    assert lp.m_alpha(sage.QQ(2) / 3) == 1 and lp.m_alpha(sage.QQ(4) / 3) == 3


def test_exclusion_7_3_14_3():
    # the type 7-(3,14)-3 is excluded by the linear program / 型 7-(3,14)-3 は線形計画で除かれる
    p, z, d, s = 7, 3, 14, 3
    th = thetas(p)
    bZ = lp.basis(tb.spaces_for(p, z * (p - 1), s, th)[0])
    bF = lp.basis(tb.spaces_for(p, d, s, th)[0])
    assert lp.solve(p, d, z * (p - 1), s, bF, bZ, lp.class_counts(p, d, s)) is None
