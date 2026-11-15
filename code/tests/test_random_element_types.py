"""Tests of rank32.random_element_types: types of trivial elements, conjugacy classes, and reading the Gram files.

Needs SageMath and numpy. Time: a few seconds.

rank32.random_element_types の試験：自明な元の型、共役類、Gram のファイルの読み込み。

SageMath と numpy が要る。時間：数秒。
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
import numpy as np  # noqa: E402
from sage.all import ZZ, matrix  # noqa: E402

from rank32 import random_element_types as rt  # noqa: E402

GRAM = json.loads((Path(__file__).with_name("data") / "extremal32_sample.json").read_text())["gram"]
G = matrix(ZZ, GRAM)


def test_fixed_type_and_powers():
    assert rt.fixed_type(G, rt.I, 3) == (32, 0)            # unimodular fixed lattice / 固定格子はユニモジュラー
    assert rt.fixed_type(G, -rt.I, 3) == (0, 0)
    assert np.array_equal(rt.mpow(-rt.I, 5), -rt.I)


def test_class_size():
    x = -rt.I
    gens = [rt.I.copy()]
    assert rt.class_size(x, gens, gens) == 1


def test_grams(tmp_path):
    (tmp_path / "code.json").write_text(json.dumps({"c1": {"gram": GRAM}}))
    (tmp_path / "kmax24.jsonl").write_text(json.dumps({"idx": 7, "gram": GRAM}) + "\n")
    (tmp_path / "reps-x.jsonl").write_text(json.dumps({"id": "r1", "gram": GRAM}) + "\n")
    out = rt.grams(tmp_path / "code.json", tmp_path / "kmax24.jsonl", tmp_path)
    assert set(out) == {"code:c1", "kmax24:7", "r1"}
