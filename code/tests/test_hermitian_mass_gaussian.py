"""Tests of rank32.hermitian_mass_gaussian: the dynamic programming against brute force, known values, small ranks.

Checks: the statistics of isotropic submodules by dynamic programming agree with full enumeration on
small 2-adic forms; the mean number of roots (A_r(A_1) / mu_r) for r = 12, 16; A_8(E_8^g) = A_12(E_8[i]) = mu_4;
m(emptyset) = 0 for r = 4, 8; the masses of rank 8 sum to mu_8; the cache survives the JSON round trip.
Needs: Python 3 and numpy. Time: about 10 seconds.

rank32.hermitian_mass_gaussian の試験：動的計画と全列挙の照合、既知の値、小さい階数。

確かめるもの：等方部分加群の統計が、小さい 2 進の形式で動的計画と全列挙とで一致すること。ルートの本数の平均
（A_r(A_1) / mu_r）の r = 12、16 の値。A_8(E_8^g) = A_12(E_8[i]) = mu_4。r = 4、8 で m(emptyset) = 0。
階数 8 の質量の和が mu_8。キャッシュが JSON の往復で変わらないこと。
必要なもの：Python 3 と numpy。時間：約 10 秒。
"""
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

from rank32 import hermitian_mass_gaussian as hm  # noqa: E402

D11, D13, D21 = ("d", 1, 1), ("d", 1, 3), ("d", 2, 1)
H1, P1, P0 = ("H", 1), ("P", 0, 1), ("P", 0, 0)


@pytest.mark.parametrize("b", [[D11], [D11] * 3, [D11, D13], [D21, D11], [H1], [H1, D11], [P1], [P1, P0], [P0, D11]])
def test_dynamic_programming_against_brute_force(b):
    assert hm.stats_cached(tuple(sorted(b))) == hm.brute_stats(b)


def test_known_values():
    assert hm.A_value({("A", 1): 1}, 16) / hm.S_mass(16) == Fr(16320, 3617)
    assert hm.A_value({("A", 1): 1}, 12) / hm.S_mass(12) == Fr(65520, 691)
    assert hm.A_value({("Eg", 8): 1}, 8) == hm.S_mass(4) == hm.A_value({("E", 8): 1}, 12)


def test_small_ranks(tmp_path):
    hm.OUT = tmp_path
    assert hm.solve(4, verbose=False) == 0 and hm.solve(8, verbose=False) == 0
    M = hm.table(8)
    assert sum(M.values()) == hm.S_mass(8)
    nz = {k: v for k, v in M.items() if v}
    assert len(nz) == 3 and all(v > 0 for v in nz.values())
    stats, trans, trans2 = hm._cache_from_json(json.loads((tmp_path / "cache-r8.json").read_text()))
    assert stats == hm._stats_store and trans == hm._trans and trans2 == hm._trans2
