"""Tests of rank32.density_sieve: the density bound and the conditions of [Neb16, Remark 3.2]
(the reference is given in rank32.density_sieve).

Time: under a second.

rank32.density_sieve の試験：密度の上界と [Neb16, Remark 3.2] の条件（文献は rank32.density_sieve に書いてある）。

時間：1 秒未満。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

from rank32 import density_sieve as ds  # noqa: E402


def test_unimodular_minimum_4():
    # the density bound allows an even unimodular lattice of minimum 4 in rank 32 and excludes one in rank 16;
    # rank 0 is never excluded
    # 密度の上界は、最小 4 の偶ユニモジュラー格子を階数 32 では許し、階数 16 では除く。階数 0 は除かない
    assert ds.ok(32, 3, 0, 4)
    assert not ds.ok(16, 3, 0, 4)
    assert ds.ok(0, 3, 5, 4)


def test_known_parts():
    # sqrt(3) E8 (rank 8, minimum 6, determinant 3^8) and K_12 (rank 12, minimum 4, determinant 3^6) pass
    # sqrt(3) E8（階数 8、最小 6、行列式 3^8）と K_12（階数 12、最小 4、行列式 3^6）は通る
    assert ds.ok(8, 3, 8, 6)
    assert ds.ok(12, 3, 6, 4)
    assert not ds.ok(8, 3, 0, 4)


def test_remark32():
    assert not ds.remark32(3, 1, 30, 0, 4)       # s = 0 needs z even / s = 0 なら z は偶数
    assert ds.remark32(3, 12, 8, 8, 4)           # s = d = 8: sqrt(3) E8 / s = d = 8：sqrt(3) E8
    assert not ds.remark32(3, 13, 6, 6, 4)       # s = d = 6 is not a multiple of 8 / s = d = 6 は 8 の倍数でない
    # minimum 10: F/sqrt(3) would be even unimodular of rank 8 with minimum >= 4, which does not exist
    # 最小 10：F/sqrt(3) は階数 8・最小 4 以上の偶ユニモジュラー格子になるが、それは無い
    assert not ds.remark32(3, 12, 8, 8, 10)
    assert ds.remark32(3, 8, 24, 24, 10)


def test_gamma_upper():
    assert ds.gamma_upper(8) == 2.0 and ds.gamma_upper(24) == 4.0
    assert len(ds.CE) == 37 and ds.CE[24] == 1.01998
