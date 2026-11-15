"""Tests of rank32.oddp_e_counts: the lemma for e(h), the tables mod p, and the GAP function EofH on small groups.

C = C_15 = <g>, h = g^3 of order 5: only h itself counts (elements of order 15 have smallest odd prime 3),
so e(h) = 1. C = C_10 = <g>, h = g^2: h and g^3 count, e(h) = 2 = lemma_e(10, 5). The GAP tests need the
executable `gap`. Time: a few seconds.

rank32.oddp_e_counts の試験：e(h) の補題、法 p の表、小さい群での GAP の関数 EofH。

C = C_15 = <g>、h = g^3（位数 5）：数えるのは h だけ（位数 15 の元の最小の奇素数は 3）なので e(h) = 1。
C = C_10 = <g>、h = g^2：h と g^3 を数え、e(h) = 2 = lemma_e(10, 5)。GAP の試験には実行ファイル `gap` が要る。
時間：数秒。
"""
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

pytest.importorskip("sage.all")

from rank32 import oddp_e_counts as k  # noqa: E402


def test_lemma_e():
    assert k.lemma_e(62, 31) == 2
    assert k.lemma_e(10, 5) == 2
    assert k.lemma_e(25, 5) is None          # p^2 divides |C| / p^2 が |C| を割る
    assert k.lemma_e(15, 5) is None          # |H| = 3 has an odd prime smaller than 5 / |H| = 3 は 5 より小さい奇素数を持つ


def test_tables_mod_p():
    lines = k.line_table(3, 5)
    assert len(lines) == (5 ** 3 - 1) // 4 and len(set(lines)) == len(lines)
    assert all(k.norm_line([(c * x) % 5 for x in v], 5) == v for v in lines[:10] for c in range(1, 5))
    assert len(k.vec_table(2, 7)) == 49


@pytest.mark.skipif(shutil.which("gap") is None, reason="needs the GAP executable")
def test_eofh():
    scr = f'Read("{k.GLIB}");\ng := (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);;\nPrint("RESULT ", EofH(Group(g), g^3, 5), "\\n");\n'
    scr += 'g := (1,2,3,4,5,6,7,8,9,10);;\nPrint("RESULT ", EofH(Group(g), g^2, 5), "\\n");\nQUIT;\n'
    out = k.run_gap(scr, timeout=120)
    assert [l.split()[1] for l in out.splitlines() if l.startswith("RESULT")] == ["1", "2"]
