"""Tests of rank32.two_elementary_min4: the enumeration in ranks 8 and 12.

In rank 8 the only lattice is sqrt(2) E_8 (s = 8, |Aut| = |W(E_8)|); in rank 12 it is sqrt(2) D_12^+
(s = 12). Needs SageMath (PARI). Time: a few seconds.

rank32.two_elementary_min4 の試験：階数 8 と 12 の列挙。

階数 8 では sqrt(2) E_8 だけ（s = 8、|Aut| = |W(E_8)|）、階数 12 では sqrt(2) D_12^+ だけ（s = 12）である。
SageMath（PARI）が要る。時間：数秒。
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

pytest.importorskip("sage.all")

from rank32 import two_elementary_min4 as m4  # noqa: E402


def test_rank8():
    res = m4.enumerate_classes(8)
    assert res["n"] == 1
    c = res["classes"][0]
    assert (c["s"], c["delta"], c["aut"], c["norm4"]) == (8, 0, 696729600, 240)
    assert m4.output_name(8) == "min4-twoel-8.json" and m4.output_name(17, 15) == "min4-twoel-17-s15.json"


def test_rank12():
    res = m4.enumerate_classes(12)
    assert res["n"] == 1 and res["classes"][0]["s"] == 12


def test_extensions_of_e8_squared():
    U = m4.unimodular_min2(16)
    assert sorted(U) == ["D16+", "D8^2+", "E8^2"]
    assert all(G.det() == 1 for G in U.values())
