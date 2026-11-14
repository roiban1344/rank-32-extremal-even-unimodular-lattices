"""Tests of rank32.bliss_automorphisms: group orders and isometry tests.

Small cases: D_4 and A_4 scaled to minimum 4 (|Aut| = 1152 and 240), a change of basis of D_4
(isometric) and D_4 against A_4 (not isometric). The slow case (an extremal lattice of rank 32,
tests/data/extremal32_sample.json, against the order recorded when the lattices were classified) runs only
when RANK32_SLOW_TESTS is set.
Needs gcc, bliss and SageMath (with GAP). Time: a few seconds; the slow case about 1 minute.

rank32.bliss_automorphisms の試験：群の位数と同型判定。

小さい場合：最小 4 に拡大した D_4 と A_4（|Aut| = 1152 と 240）、D_4 の基底の取り替え（同型）、
D_4 と A_4（同型でない）。遅い場合（ランク 32 の極値格子、tests/data/extremal32_sample.json。格子を分類したときに記録した位数と比べる）は、
RANK32_SLOW_TESTS が設定されているときだけ走る。
gcc、bliss、SageMath（GAP を含む）が要る。時間：数秒。遅い場合は約 1 分。
"""
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, matrix, identity_matrix  # noqa: E402

from rank32 import bliss_automorphisms as ba  # noqa: E402
from rank32.isometry import is_isometry, all_automorphisms  # noqa: E402

D4 = matrix(ZZ, [[4, -2, 0, 0], [-2, 4, -2, -2], [0, -2, 4, 0], [0, -2, 0, 4]])
A4 = matrix(ZZ, [[4, -2, 0, 0], [-2, 4, -2, 0], [0, -2, 4, -2], [0, 0, -2, 4]])
DATA = Path(__file__).with_name("data") / "extremal32_sample.json"


@pytest.fixture(autouse=True)
def work_dir(tmp_path):
    ba.OUT = tmp_path
    yield
    ba.OUT = None


def test_small_orders():
    for G, order in ((D4, "1152"), (A4, "240")):
        s = ba.Shell(G, "small", 2, test=True)
        r = ba.automorphisms(s)
        assert r["order"] == order
        assert all_automorphisms([matrix(ZZ, g) for g in r["generators"]], G)


def test_small_isometry():
    U = identity_matrix(ZZ, 4)
    U[0, 1] = 1
    s = ba.Shell(D4, "D4", 2, test=True)
    t = ba.Shell(U.T * D4 * U, "D4-rebased", 2, test=True)
    r = ba.isometry(s, t)
    assert r["isometric"] and is_isometry(matrix(ZZ, r["transform"]), s.G, t.G)
    assert not ba.isometry(s, ba.Shell(A4, "A4", 2, test=True))["isometric"]


def test_rejects_non_extremal():
    with pytest.raises(ValueError):
        ba.Shell(D4, "D4", 2)            # not of rank 32 / 階数が 32 でない


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("RANK32_SLOW_TESTS"), reason="slow; set RANK32_SLOW_TESTS=1")
def test_extremal32():
    d = json.loads(DATA.read_text())
    G = matrix(ZZ, d["gram"])
    r = ba.automorphisms(ba.Shell(G, "rank 32", 6))
    assert r["order"] == d["aut_order"]
    assert all_automorphisms([matrix(ZZ, g) for g in r["generators"]], G)
