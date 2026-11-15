"""Tests of rank32.frame_fingerprint: the derived a_k and fingerprint from a recorded histogram, and the full run.

The histogram hist0 and the classes of minimal vectors of a lattice of the walk give the fingerprint
f816fbb83ea18073 (the recorded value). The full run on the extremal lattice of
tests/data/extremal32_sample.json (about one minute and 1 GB with PL = 4 and 2 threads) runs only
when RANK32_SLOW_TESTS is set; it checks the identities that the design forces.
Needs: numpy; the slow case SageMath and gcc with OpenMP.

rank32.frame_fingerprint の試験：記録のヒストグラムから導く a_k と指紋、全体の実行。

ランダムウォークの格子の hist0 と最短ベクトルの類から、指紋 f816fbb83ea18073（記録された値）が出ること。
tests/data/extremal32_sample.json の極値格子での全体の実行（PL = 4、2 スレッドで約 1 分・1 GB）は、
RANK32_SLOW_TESTS が設定されているときだけ走る。デザインで強制される恒等式を確かめる。
必要なもの：numpy。遅い場合は SageMath と OpenMP つきの gcc。
"""
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

from rank32 import frame_fingerprint as fh  # noqa: E402

HIST0 = [[1, 442302336], [3, 160081920], [6, 62225280], [10, 11460096], [15, 2039040], [21, 829440], [28, 925920], [45, 17280], [66, 432], [120, 270],
         [136, 1152], [190, 432]]
MINCLS = [[620, 73440]]


def test_derive_recorded():
    d = fh.derive(HIST0, sum(c for _, c in MINCLS))
    assert d["fp"] == "f816fbb83ea18073"
    assert sum(k * c for k, c in d["ak"].items()) == fh.N8_LINES
    assert d["a0"] + sum(d["ak"].values()) == fh.A_ADM


@pytest.mark.slow
@pytest.mark.skipif(not os.environ.get("RANK32_SLOW_TESTS"), reason="slow; set RANK32_SLOW_TESTS=1")
def test_run_one(tmp_path, monkeypatch):
    pytest.importorskip("sage.all")
    monkeypatch.setenv("FRAME_PL", "4")
    monkeypatch.setenv("FRAME_SCRATCH", str(tmp_path))
    fh.GRAMS["sample"] = json.loads((Path(__file__).with_name("data") / "extremal32_sample.json").read_text())["gram"]
    r = fh.run_one("sample", 2)
    assert r["n"] == 73440 and r["saturated"] == 0
    assert sum(c for _, c in r["mincls"]) == 73440       # one entry per class of a minimal line / 最短の直線の類ごとに 1 つ
    d = fh.derive(r["hist0"], 73440)
    assert sum(k * c for k, c in d["ak"].items()) == fh.N8_LINES and d["ak"][1] > 0
