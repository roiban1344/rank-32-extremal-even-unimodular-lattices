"""Tests of rank32.forms: known group orders, consistency of the returned matrices, recovery from a stalling basis.

Time: about 30 seconds.

rank32.forms の試験：既知の群の位数、返す行列の整合、無応答になる基底からの回復。

時間：約 30 秒。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

sage = pytest.importorskip("sage.all")
from sage.all import ZZ, CartanMatrix, block_diagonal_matrix, matrix, pari  # noqa: E402

from rank32 import forms  # noqa: E402

A2 = matrix(ZZ, [[2, -1], [-1, 2]])
W2 = matrix(ZZ, [[1, 1], [-1, 0]])
# a gluing in the genus (17, 3, 1); in this basis qfauto does not finish in 150 seconds, in other bases it takes about 1 second
# 種 (17, 3, 1) の貼り合わせ。この基底では qfauto が 150 秒で終わらず、他の基底では約 1 秒
STALL = matrix(ZZ, pari("[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,2,-1,1,1,-1,-1,0,0,0,0,0,0,-1,0,0,-1;0,-1,2,0,0,1,1,0,0,0,0,0,0,1,0,0,0;0,1,0,2,0,0,0,0,0,0,0,0,0,-1,0,0,0;0,1,0,0,2,0,0,0,0,0,0,0,0,0,0,0,-1;0,-1,1,0,0,2,1,0,0,0,0,0,0,1,0,0,0;0,-1,1,0,0,1,2,0,0,0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,2,-1,1,-1,-1,1,1,-1,-1,1;0,0,0,0,0,0,0,-1,2,0,0,1,0,-1,1,1,-1;0,0,0,0,0,0,0,1,0,2,0,0,1,1,0,0,1;0,0,0,0,0,0,0,-1,0,0,2,1,0,0,1,1,0;0,0,0,0,0,0,0,-1,1,0,1,2,0,-1,1,1,-1;0,0,0,0,0,0,0,1,0,1,0,0,2,1,0,0,0;0,-1,1,-1,0,1,0,1,-1,1,0,-1,1,4,-1,0,0;0,0,0,0,0,0,0,-1,1,0,1,1,0,-1,2,1,0;0,0,0,0,0,0,0,-1,1,0,1,1,0,0,1,2,-1;0,-1,0,0,-1,0,1,1,-1,1,0,-1,0,0,0,-1,4]"))
STALL2 = matrix(ZZ, pari("[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,2,-1,1,1,-1,0,0,0,0,0,0,0,-1,1,0,0;0,-1,2,0,0,1,0,0,0,0,0,0,0,1,-1,0,0;0,1,0,2,0,0,0,0,0,0,0,0,0,-1,0,0,0;0,1,0,0,2,0,0,0,0,0,0,0,0,0,1,0,0;0,-1,1,0,0,2,0,0,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,2,-1,1,1,-1,1,-1,1,1,-1,-1;0,0,0,0,0,0,-1,2,0,0,1,0,1,-1,-1,1,0;0,0,0,0,0,0,1,0,2,0,0,1,0,1,1,0,-1;0,0,0,0,0,0,1,0,0,2,0,1,0,0,0,0,-1;0,0,0,0,0,0,-1,1,0,0,2,0,1,-1,-1,1,0;0,0,0,0,0,0,1,0,1,1,0,2,0,1,1,0,-1;0,0,0,0,0,0,-1,1,0,0,1,0,2,-1,-1,1,0;0,-1,1,-1,0,1,1,-1,1,0,-1,1,-1,4,1,-1,0;0,1,-1,0,1,0,1,-1,1,0,-1,1,-1,1,4,0,-1;0,0,0,0,0,0,-1,1,0,0,1,0,1,-1,0,2,0;0,0,0,0,0,0,-1,0,-1,-1,0,-1,0,0,-1,0,2]"))  # isometric to STALL; qfisom(STALL, STALL2) also takes more than 30 seconds / STALL と同型。qfisom(STALL, STALL2) も 30 秒超


@pytest.fixture(autouse=True)
def keep_pari_random_state():
    """qfauto of PARI uses random numbers internally, which changes the generators it picks; restore the state so that later tests are not affected.

    PARI の qfauto は内部で乱数を消費し、生成元の取り方が変わる。後続の試験に影響を与えないよう、状態を戻す。"""
    state = pari.getrand()
    yield
    pari.setrand(state)


def test_orders():
    """|W(E8)| = 696729600."""
    E8 = matrix(ZZ, CartanMatrix(["E", 8]))
    order, gens = forms.qfauto(E8)
    assert order == 696729600
    assert all(g.transpose() * E8 * g == E8 for g in gens)


def test_two_forms():
    for k, order in ((1, 6), (2, 72)):
        G = block_diagonal_matrix([A2] * k, subdivide=False)
        W = block_diagonal_matrix([W2] * k, subdivide=False)
        assert forms.qfauto([G, G * W.transpose()])[0] == order


def test_isometry_returns_original_basis():
    G = matrix(ZZ, CartanMatrix(["E", 6]))
    T = matrix(ZZ, 6, 6, lambda i, j: 1 if i == j else (1 if j == i + 1 else 0))
    H = T.transpose() * G * T
    S = forms.qfisom(G, H)
    assert S.transpose() * H * S == G and abs(S.det()) == 1
    assert forms.qfisom(G, matrix(ZZ, CartanMatrix(["D", 6]))) is None


def test_random_basis_keeps_lattice():
    import random
    rng = random.Random(1)
    for _ in range(5):
        T = forms._random_basis([STALL], rng)
        H = T.transpose() * STALL * T
        assert abs(T.det()) == 1 and forms._maxdiag(H) <= forms._maxdiag(STALL)


def test_recovers_from_stalling_basis():
    order, gens = forms.qfauto(STALL, limit=2)
    assert order == 85614133248000
    assert all(g.transpose() * STALL * g == STALL for g in gens)
    S = forms.qfisom(STALL, STALL2, limit=2)
    assert S.transpose() * STALL2 * S == STALL


def test_interrupt_does_not_keep_frames():
    """A call that succeeds after a timeout leaves no frame holding a PARI result.

    cysignals keeps the interrupt exception until the next interrupt and releases it inside the
    next sig_on. If its traceback holds the frame of _retry (the local r is a result on the PARI
    stack), the release happens inside sig_on and cypari2 raises
    SystemError('calling remove_from_pari_stack() inside sig_on()') (3 of 2400 rows in one run of the classification).

    時間切れの後に成功した呼び出しが、PARI の結果を抱えたフレームを残さない。

    cysignals は割り込みの例外を次の割り込みまで握り、次の sig_on の中で手放す。そのトレースバックが
    _retry のフレーム（局所変数 r は PARI スタック上の結果）を抱えていると、解放が sig_on の中で起き、
    cypari2 が SystemError('calling remove_from_pari_stack() inside sig_on()') を出す（分類の計算の 1 回の実行で、2400 行のうち 3 行）。"""
    import gc
    forms.qfauto(STALL, limit=2)  # the original basis times out, another basis succeeds / 元の基底は時間切れ、別の基底で成功する
    frames = [o for o in gc.get_objects() if type(o).__name__ == "frame" and o.f_code is forms._retry.__code__]
    assert frames == []


def test_system_error_is_retried(monkeypatch):
    """A SystemError of cypari2 is retried; if it repeats, QfTimeout (no verdict is returned).

    cypari2 の SystemError は再試行し、続いたら QfTimeout（判定を返さない）。"""
    E8 = matrix(ZZ, CartanMatrix(["E", 8]))
    pari_ = forms._pari
    calls = []

    def flaky(forms):
        calls.append(1)
        if len(calls) == 1:
            raise SystemError("calling remove_from_pari_stack() inside sig_on()")
        return pari_(forms)

    monkeypatch.setattr(forms, "_pari", flaky)
    assert forms.qfauto(E8)[0] == 696729600 and len(calls) == 2

    def broken(forms):
        raise SystemError("calling remove_from_pari_stack() inside sig_on()")

    monkeypatch.setattr(forms, "_pari", broken)
    with pytest.raises(forms.QfTimeout):
        forms.qfauto(E8)
    with pytest.raises(forms.QfTimeout):
        forms.qfisom(E8, E8)


def test_timeout_is_reported():
    with pytest.raises(forms.QfTimeout):
        forms.qfauto(STALL, limit=0.5, tries=1)


def test_short_vectors():
    E8 = matrix(ZZ, CartanMatrix(["E", 8]))
    V = forms.short_vectors(E8, 2)
    assert V.shape == (120, 8)
    G = np.array(E8, dtype=np.int64)
    assert (np.einsum("ij,jk,ik->i", V, G, V) == 2).all()
    assert len(np.unique(V, axis=0)) == 120
    assert forms.short_vectors(2 * E8, 2).shape == (0, 8)
    # a single pair +-v (PARI prints Mat([...])) / +-1 本だけのとき（PARI の表示が Mat([...]) になる）
    V1 = forms.short_vectors(matrix(ZZ, [[2, 0], [0, 4]]), 2)
    assert V1.shape == (1, 2) and abs(V1).tolist() == [[1, 0]]
