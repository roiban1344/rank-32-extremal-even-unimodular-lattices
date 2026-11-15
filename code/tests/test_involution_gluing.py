"""Tests of rank32.involution_gluing, rank32.involution_gluing_checks and their entry points on the pair F = K = sqrt(2) E_8.

The image of Aut(sqrt(2) E_8) = W(E_8) in O(D) is O^+(8, 2) of order 348364800 (the kernel is +-1),
so there is a single double coset, with |Stab| = 348364800 in the image and mass 1/348364800. The
glued lattice is even unimodular of rank 16 (minimum 2). The test writes the list of components and
groups.json in a temporary directory with produce/involutions/discriminant_action.py, closes the
pair with the search, and rechecks it with crosscheck/involutions/involution_gluing_verify.py.
Needs gcc and SageMath (PARI, GAP). Time: under a minute.

rank32.involution_gluing、rank32.involution_gluing_checks と入口の試験（組 F = K = sqrt(2) E_8 で）。

Aut(sqrt(2) E_8) = W(E_8) の O(D) への像は位数 348364800 の O^+(8, 2)（カーネルは +-1）なので、両側剰余類は 1 つで、
像の中の |Stab| = 348364800、質量 1/348364800 である。貼り合わせた格子は階数 16 の偶ユニモジュラー格子
（最小 2）になる。試験は一時ディレクトリに成分の一覧と groups.json を discriminant_action.py で書き、探索で組を
閉じ、involution_gluing_verify.py で recheck する。gcc と SageMath（PARI、GAP）が要る。時間：1 分以内。
"""
import importlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))   # code/ as the import root / code/ を import の起点にする

pytest.importorskip("sage.all")

from sage.all import ZZ, CartanMatrix, matrix  # noqa: E402

from rank32 import discriminant_action as da  # noqa: E402
from rank32 import involution_gluing as gl  # noqa: E402
from rank32 import involution_gluing_checks as vv  # noqa: E402

ORDER = 348364800


def entry(rel, name):
    """Import an entry script as a module of the namespace packages produce/ and crosscheck/.

    入口のスクリプトを、名前空間パッケージ produce/、crosscheck/ のモジュールとして import する。"""
    return importlib.import_module(rel[:-3].replace("/", "."))


@pytest.fixture(scope="module")
def tables(tmp_path_factory):
    d = tmp_path_factory.mktemp("gluing")
    G = 2 * matrix(ZZ, CartanMatrix(["E", 8]))
    rec = dict(gram=[[int(x) for x in r] for r in G.rows()], s=8, delta=0, origin="sqrt(2) E8", aut=696729600, norm4=240)
    (d / "parts.json").write_text(json.dumps(dict(d=8, smin=0, n=1, classes=[rec])))
    old = dict(da.PARTS)
    da.PARTS.clear(); da.PARTS["e8"] = ("parts.json", 0)
    try:
        entry("produce/involutions/discriminant_action.py", "discriminant_action_entry").main(d, d / "groups.json")
    finally:
        da.PARTS.clear(); da.PARTS.update(old)
    gr = json.loads((d / "groups.json").read_text())
    assert gr["e8"]["image"] == str(ORDER) and gr["e8"]["kernel"] == 2 and gr["e8"]["orbits"] == [1, 120, 135]
    gr["e8"]["O_q"] = str(ORDER)
    (d / "groups.json").write_text(json.dumps(gr))
    gl.load_tables(d / "groups.json", d / "none.json", d / "none.json", d)
    gl.PAIRS["e8-e8"] = ("e8", "e8")
    return d


def test_search_closes(tables):
    gl.search("e8-e8", 60, tables)
    res = json.loads((tables / "e8-e8.json").read_text())
    assert res["closed"] and res["n_classes"] == 1
    assert Fraction(res["target_image"]) == Fraction(1, ORDER) and Fraction(res["M_up"]) == Fraction(1, 4 * ORDER)
    assert res["classes"][0]["stab_image"] == ORDER and not res["classes"][0]["extremal"]
    assert (tables / "e8-e8.ckpt.json").exists() and (tables / "e8-e8.runs.jsonl").exists()


def test_build_gram_and_glue_of(tables):
    res = json.loads((tables / "e8-e8.json").read_text())
    GF_, F, gF = gl.part("e8")
    G, s = vv.build_gram(GF_, GF_, res["classes"][0]["phi"])
    assert G.nrows() == 16 and G.det() == 1 and vv.min_info(G)[0] == 2
    assert vv.element_order(s) == 2
    phi = vv.glue_of(G, s, GF_, GF_)
    g = gl.glue.Glue(F, F, gF["gens_disc"], gF["gens_disc"], int(gF["image"]), int(gF["image"]))
    assert g.is_antiisometry(phi)
    assert g.canon(phi)[0] == g.canon(np.array(res["classes"][0]["phi"], dtype=np.int64))[0]


def test_recheck(tables):
    cv = entry("crosscheck/involutions/involution_gluing_verify.py", "involution_gluing_verify_entry")
    vv.LISTS = vv.OUT = tables
    cv.mode_recheck("e8-e8")
    r = json.loads((tables / "verify-recheck-e8-e8.json").read_text())
    assert r["classes"] == 1 and Fraction(r["M_up"]) == Fraction(1, 4 * ORDER)
