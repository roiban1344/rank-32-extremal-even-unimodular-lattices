"""Tests of rank32.even_neighbour_walk: the replay of a step from an extremal lattice of rank 32.

From the extremal lattice of tests/data/extremal32_sample.json, replay(G, w, s) is deterministic and gives
an even unimodular lattice; an admissible class has norm divisible by 4 and its mod-8 lift norm divisible by 8.
Needs SageMath (with fpylll). Time: about 10 seconds.

rank32.even_neighbour_walk の試験：階数 32 の極値格子からの 1 歩の再現。

tests/data/extremal32_sample.json の極値格子から、replay(G, w, s) は決定的で、偶ユニモジュラー格子を与える。
許容類のノルムは 4 で割れ、その mod 8 の持ち上げのノルムは 8 で割れる。
SageMath（fpylll を含む）が要る。時間：約 10 秒。
"""
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, matrix  # noqa: E402

from rank32 import even_neighbour_walk as ew  # noqa: E402

G = matrix(ZZ, json.loads((Path(__file__).with_name("data") / "extremal32_sample.json").read_text())["gram"])


def test_admissible_class():
    rng = random.Random(3)
    v = ew.random_admissible_class(G, rng)
    assert (v * G * v) % 4 == 0 and v != 0
    w = ew.mod8_rep(G, v)
    assert (w * G * w) % 8 == 0 and all((a - b) % 2 == 0 for a, b in zip(v, w))


def test_replay_deterministic():
    a, b = ew.replay(G, 0, 1), ew.replay(G, 0, 1)
    assert a == b
    assert a.det() == 1 and all(a[i, i] % 2 == 0 for i in range(32))
