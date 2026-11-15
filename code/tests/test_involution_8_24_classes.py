"""Tests of rank32.involution_8_24_classes: the C part on the root graph of E_8, the hyperbolic basis, and the gluing with sqrt(2) E_8.

The graph of the 240 roots of E_8 (edges coloured by inner products, vertices by the norm) has
automorphism group W(E_8) of order 696729600. The maximal singular subspaces of F_2^8 of plus type
are 270. Needs gcc, bliss and SageMath (PARI, GAP). Time: a few seconds.

rank32.involution_8_24_classes の試験：E_8 のルートのグラフでの C の部分、双曲基底、sqrt(2) E_8 との貼り合わせ。

E_8 の 240 本のルートのグラフ（辺の色は内積、頂点の色はノルム）の自己同型群は位数 696729600 の W(E_8) である。
+ 型の F_2^8 の極大特異部分空間は 270 個である。gcc、bliss、SageMath（PARI、GAP）が要る。時間：数秒。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

pytest.importorskip("sage.all")

from sage.all import ZZ, CartanMatrix, matrix  # noqa: E402

from rank32 import forms  # noqa: E402
from rank32 import involution_8_24_classes as cl  # noqa: E402

C = matrix(ZZ, CartanMatrix(["E", 8]))


def test_standard_subspaces():
    assert cl.SUBS.shape == (270, 15)


def test_canonical_form_of_the_root_graph():
    V = forms.short_vectors(C, 2)
    V = np.vstack([V, -V])
    Cn = np.array(C.list(), dtype=np.int64).reshape(8, 8)
    lib = cl.load_canon()
    h, lab, gens = cl.certificate(lib, V @ Cn @ V.T)
    assert cl.group_order(gens, len(V)) == 696729600
    perm = np.random.default_rng(1).permutation(len(V))         # the hash does not depend on the order of the vertices / ハッシュは頂点の順に依らない
    W = V[perm]
    assert cl.certificate(lib, W @ Cn @ W.T)[0] == h


def test_hyperbolic_basis():
    De = np.array(C.inverse().change_ring(ZZ).list(), dtype=np.int64).reshape(8, 8)
    b = cl.hyperbolic_basis(De)
    assert len(b) == 8
    q = lambda a: int(a @ De @ a) // 2 & 1
    bb = lambda a, c: int(a @ De @ c) & 1
    for i in range(4):
        assert q(b[i]) == 0 and q(b[4 + i]) == 0
        for j in range(4):
            assert bb(b[i], b[4 + j]) == (i == j) and bb(b[i], b[j]) == 0 and bb(b[4 + i], b[4 + j]) == 0
