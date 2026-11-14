"""Tests of rank32.io: the exact bytes of each format, and no temporary file left behind.

rank32.io の試験：書式ごとの正確なバイト列と、一時ファイルが残らないこと。
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # code/ as the import root / code/ を import の起点にする

from rank32 import io as rio  # noqa: E402

OBJ = {"b": [1, 2, {"x": "1/3"}], "a": "é", "n": None}


def test_formats(tmp_path):
    cases = {
        "write_json_atomic": json.dumps(OBJ, ensure_ascii=False, indent=1) + "\n",
        "write_json_atomic_indent1": json.dumps(OBJ, ensure_ascii=False, indent=1),
        "write_json_atomic_compact": json.dumps(OBJ),
        "write_json_atomic_indent0": json.dumps(OBJ, ensure_ascii=False, indent=0),
    }
    for name, text in cases.items():
        p = tmp_path / f"{name}.json"
        getattr(rio, name)(p, OBJ)
        assert p.read_text() == text
        assert json.loads(p.read_text()) == OBJ
    p = tmp_path / "plain.txt"
    rio.write_text_atomic(p, "line\n")
    assert p.read_text() == "line\n"
    assert not list(tmp_path.glob("*.tmp"))


def test_overwrite(tmp_path):
    p = tmp_path / "a.json"
    rio.write_json_atomic(p, [1])
    rio.write_json_atomic(p, [2])
    assert json.loads(p.read_text()) == [2]
