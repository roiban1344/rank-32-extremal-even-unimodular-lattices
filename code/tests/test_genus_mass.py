"""Tests of rank32.genus_mass and rank32.rootless_unimodular_mass: known masses of Hermitian and Z-genera.

Time: a few seconds.

rank32.genus_mass と rank32.rootless_unimodular_mass の試験：エルミートの種と Z の種の既知の質量。

時間：数秒。
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import QQ  # noqa: E402

from rank32 import genus_mass as gm  # noqa: E402
from rank32 import rootless_unimodular_mass as rm  # noqa: E402


def test_hkn09():
    # [HKN09, Theorem 1]: the masses of the Eisenstein lattices of rank 4 and 8 / Eisenstein 格子の階数 4 と 8 の質量
    assert gm.mass_H(3, 4, 0) == QQ(1) / (2 ** 7 * 3 ** 5 * 5)
    assert gm.mass_H(3, 8, 0) == QQ(1) / (2 ** 15 * 3 ** 10 * 5 ** 2)


def test_rank_one_and_two():
    for p in (3, 5, 7, 11, 13):
        assert gm.mass_U(p, 1) == QQ(1) / (2 * p)
    assert gm.mass_U(5, 2) == QQ(1) / (8 * 25)
    assert gm.mass_H(3, 3, 0) is None   # z - s odd / z - s が奇数


def test_e8_pairs():
    # E8: 2240 elements of type 3-(1,6)-1, |C| = 600 for 5-(2,0)-0 / E8：型 3-(1,6)-1 の元は 2240 個、5-(2,0)-0 は |C| = 600
    assert gm.M_pair(3, 8, 1, 1) == QQ(2240) / 696729600
    assert gm.M_pair(5, 8, 2, 0) == QQ(1) / 600
    assert gm.mass_Z_genus(3, 8, 0) == QQ(1) / 696729600


def test_genus_conditions():
    assert gm.genus_exists(3, 12, 6) and gm.genus_exists(3, 8, 8) and not gm.genus_exists(3, 6, 6)
    assert gm.order_O(3, 2, 1) == 4 and gm.order_O(3, 1, 0) == 2
    assert gm.isotropic_count(3, 4, 2) == 8   # 2(q + 1) totally isotropic planes of a split space of dimension 4 / 分裂する 4 次元空間の全等方平面は 2(q + 1) 個


def test_rootless_masses():
    m = rm.noroot_masses(3, 7)
    # p = 3: no rootless unimodular lattice in rank 1..5 and 7; in rank 6 only K_12
    # p = 3：階数 1..5 と 7 にルート無しのユニモジュラー格子は無い。階数 6 は K_12 だけ
    assert m[0] == 1 and all(m[z] == 0 for z in (1, 2, 3, 4, 5, 7))
    assert m[6] == QQ(1) / 39191040
    assert all(x == 0 for x in rm.noroot_masses(5, 3)[1:]) and rm.noroot_masses(5, 4)[4] > 0
