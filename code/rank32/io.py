"""Atomic writing of text and JSON files, and the place of the release assets.

Each function writes to PATH.tmp and then renames it to PATH, so a reader (or a run that
is stopped and started again) never sees a half-written file. The functions differ only in
the format of the JSON text; each keeps the format of the programs that use it, so that
their outputs stay byte for byte the same:
  write_text_atomic(path, text)            the text as given
  write_json_atomic(path, obj)             indent 1, non-ASCII kept, newline at the end
  write_json_atomic_indent1(path, obj)     indent 1, non-ASCII kept, no newline at the end
  write_json_atomic_compact(path, obj)     json.dumps defaults (one line, ASCII), no newline
  write_json_atomic_indent0(path, obj)     indent 0, non-ASCII kept, no newline at the end
asset(*parts) is the default place of the inputs that come from the release assets.

テキストと JSON のファイルをアトミックに書く。また、Release のアセットの置き場所を与える。

どの関数も PATH.tmp に書いてから PATH に名前を変える。そのため、読む側（や、止めて
やり直した実行）が書きかけのファイルを見ることはない。関数の違いは JSON のテキストの書式だけ
である。使うプログラムの出力がバイト単位で変わらないよう、それぞれの書式を保つ。
  write_text_atomic(path, text)            テキストをそのまま
  write_json_atomic(path, obj)             字下げ 1、ASCII 以外はそのまま、末尾に改行
  write_json_atomic_indent1(path, obj)     字下げ 1、ASCII 以外はそのまま、末尾に改行なし
  write_json_atomic_compact(path, obj)     json.dumps の既定（1 行、ASCII）、改行なし
  write_json_atomic_indent0(path, obj)     字下げ 0、ASCII 以外はそのまま、末尾に改行なし
asset(*parts) は、Release のアセットから来る入力の既定の場所である。
"""
import json
import os
from pathlib import Path


def asset(*parts):
    """A path under the directory where the release assets are unpacked: $RANK32_ASSETS, or assets/ under the current directory when RANK32_ASSETS is not set.

    Release のアセットを展開したディレクトリの下のパス。そのディレクトリは $RANK32_ASSETS で、設定されていなければ
    今のディレクトリの下の assets/ である。"""
    return Path(os.environ.get("RANK32_ASSETS", "assets"), *parts)


def write_text_atomic(path, text):
    """Write text to path atomically.

    テキストを path にアトミックに書く。"""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text)
    os.rename(tmp, path)


def write_json_atomic(path, obj):
    """Write obj as JSON (indent 1, non-ASCII kept, final newline) to path atomically.

    obj を JSON（字下げ 1、ASCII 以外はそのまま、末尾に改行）で path にアトミックに書く。"""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=1) + "\n")
    os.rename(tmp, path)


def write_json_atomic_indent1(path, obj):
    """Write obj as JSON (indent 1, non-ASCII kept, no final newline) to path atomically.

    obj を JSON（字下げ 1、ASCII 以外はそのまま、末尾に改行なし）で path にアトミックに書く。"""
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=1))
    os.replace(tmp, path)


def write_json_atomic_compact(path, obj):
    """Write obj as JSON (json.dumps defaults, no final newline) to path atomically.

    obj を JSON（json.dumps の既定、末尾に改行なし）で path にアトミックに書く。"""
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj))
    os.replace(tmp, path)


def write_json_atomic_indent0(path, obj):
    """Write obj as JSON (indent 0, non-ASCII kept, no final newline) to path atomically.

    obj を JSON（字下げ 0、ASCII 以外はそのまま、末尾に改行なし）で path にアトミックに書く。"""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=0))
    tmp.replace(path)
