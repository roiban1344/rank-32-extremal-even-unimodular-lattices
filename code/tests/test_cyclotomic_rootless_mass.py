"""Tests of rank32.cyclotomic_rootless_mass: the tables of root systems sum to mass_H, and the value E^2 = 1/720000.

Time: a few seconds.

rank32.cyclotomic_rootless_mass の試験：ルート系の表の和が mass_H になること、E^2 = 1/720000 の値。

時間：数秒。
"""
import sys
from fractions import Fraction as Fr
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")

from rank32 import cyclotomic_rootless_mass as cm  # noqa: E402


@pytest.mark.parametrize("p,z,s", [(5, 4, 4), (5, 5, 5), (5, 6, 4), (7, 4, 4), (17, 2, 0), (5, 4, 0)])
def test_table_sums_to_mass(p, z, s):
    M = cm.table(p, z, s)
    assert sum(M.values()) == cm.mass_H(p, z, s)
    assert M[(0, 0)] == cm.rootless(p, z, s)


def test_e8_squared():
    # p = 5, genus (4, 0): the lattice E8 + E8 with its structure gives m(E^2) = 1/720000 = 1/(600^2 * 2)
    # p = 5、種 (4, 0)：構造つきの E8 + E8 が m(E^2) = 1/720000 = 1/(600^2 * 2) を与える
    assert cm.table(5, 4, 0)[(0, 2)] == Fr(1, 720000)


def test_finite_spaces():
    assert cm.order_O(5, 2, 1) == 8 and cm.sing(5, 2, 1) == 8 and cm.gl(5, 2) == 480
    assert cm.leg(2, 5) == -1 and cm.leg(4, 5) == 1 and cm.leg(10, 5) == 0
    assert cm.mass_H(5, 3, 0) == 0
