# Code

English | [日本語](README.ja.md)

The programs behind the computations of the papers in this repository. They are kept as they were
run, with the same functions, constants and random seeds. Comments are in English and Japanese.
For now `code/` holds the shared library and its tests; the programs of the computations follow.

## Layout

| Directory | Contents |
|---|---|
| `rank32/` | The shared library. It imports nothing else from `code/`. |
| `tests/` | Small cases, run with `pytest`. |

## Requirements

The programs were run with:

- SageMath 10.9 on Python 3.12, with PARI/GP 2.17.3 and GAP 4.15.1 (through `libgap`);
- bliss 0.77 and gcc, for the C parts;
- numpy and pytest.

Environment variables:

| Variable | Use | Default |
|---|---|---|
| `RANK32_BUILD` | where the C parts are compiled on first use | `$XDG_CACHE_HOME/rank32`, that is `~/.cache/rank32` |
| `BLISS_PREFIX` | where bliss is installed (`include/` and `lib/`) | `$CONDA_PREFIX`, else the prefix of the running Python |
| `RANK32_TMP` | where `rank32.canonical` writes its temporary graph files | the system's temporary directory |
| `RANK32_SLOW_TESTS` | when set to any non-empty value, the tests marked `slow` run too | not set |

## Tests

From the repository root:

```
sage -python -m pytest -q code/tests
```

Tests marked `slow` take a minute or more and more memory; they are skipped unless `RANK32_SLOW_TESTS` is set.
Tests that need SageMath are skipped when it is not available.

## License

The code in `code/` is licensed under the MIT license; see [`code/LICENSE`](LICENSE).
