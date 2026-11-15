"""Shared settings of the tests of code/.

Tests marked `slow` (a minute or more, more memory) run only when RANK32_SLOW_TESTS is set; otherwise they are
skipped here.

code/ の試験の共通の設定。
マーカー `slow` の付いた試験（1 分以上かかり、メモリも多く使う）は、RANK32_SLOW_TESTS が設定されているときだけ走る。
設定されていなければ、ここで飛ばす。
"""
import os

import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: takes a minute or more; runs only when RANK32_SLOW_TESTS is set")


def pytest_collection_modifyitems(config, items):
    if os.environ.get("RANK32_SLOW_TESTS"):
        return
    skip = pytest.mark.skip(reason="slow; set RANK32_SLOW_TESTS=1")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)
