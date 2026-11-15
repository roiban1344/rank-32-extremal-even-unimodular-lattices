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

`environment.lock.txt` fixes this environment: it lists, with their MD5 digests, the 433 conda-forge packages (linux-64)
of the environment in which the programs were run, including some that the programs do not use, such as JupyterLab.
There are two ways to set it up; run the commands from the repository root.

1. Create the environment with micromamba (with conda: `conda create -n rank32 --file code/environment.lock.txt`):

   ```
   micromamba create -n rank32 -f code/environment.lock.txt
   micromamba run -n rank32 python -m pytest -q code/tests
   ```

2. Build a container from `Dockerfile`, which installs the same packages on a pinned micromamba image (about 10 GB
   unpacked):

   ```
   docker build -t rank32 code
   docker run --rm -v "$PWD":/work:ro -w /work rank32 python -m pytest -q -p no:cacheprovider code/tests
   ```

`-p no:cacheprovider` is there because the repository is mounted read-only (`:ro`), so pytest cannot write its cache.

The packages are for linux-64, so the first way works only on linux-64 machines. On other machines, add
`--platform linux/amd64` to both `docker build` and `docker run`; the container then runs under emulation, more slowly.

Environment variables:

| Variable | Use | Default |
|---|---|---|
| `RANK32_BUILD` | where the C parts are compiled on first use | `$XDG_CACHE_HOME/rank32`, that is `~/.cache/rank32` |
| `BLISS_PREFIX` | where bliss is installed (`include/` and `lib/`) | `$CONDA_PREFIX`, else the prefix of the running Python |
| `RANK32_TMP` | where `rank32.canonical` writes its temporary graph files | the system's temporary directory |
| `RANK32_SLOW_TESTS` | when set to any non-empty value, the tests marked `slow` run too | not set |

## Tests

In the environment above, from the repository root:

```
python -m pytest -q code/tests
```

Tests marked `slow` take a minute or more and more memory; they are skipped unless `RANK32_SLOW_TESTS` is set.
Tests that need SageMath are skipped when it is not available.

## License

The code in `code/` is licensed under the MIT license; see [`code/LICENSE`](LICENSE).
