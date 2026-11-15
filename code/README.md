# Code

English | [日本語](README.ja.md)

The programs behind the computations of the papers in this repository. They are kept as they were
run, with the same functions, constants and random seeds. The one exception is that labels the outputs wrote in
Japanese (such as the verdict of a linear program) are now English constants with the same meaning.
Comments are in English and Japanese.
Not all programs are here yet; the others follow.

## Layout

| Directory | Contents |
|---|---|
| `rank32/` | The shared library. It imports nothing else from `code/`. |
| `produce/<topic>/` | The programs that produced the lists, certificates and data, by topic of the first paper; so far `element_types` (elements of odd prime order), `involutions` (involution types), `gaussian` (Hermitian lattices over Z[i]) and `known_lattices`. The docstring of each file names its computations (C1, C2, ... of Table 11), its inputs (if any), its outputs and its running time. |
| `crosscheck/<topic>/` | Checks of the results of `produce/`, some of them by independent methods. They do not import `produce/`. |
| `toolcheck/` | Checks of the tools on known cases (computations C37 and C38) and auxiliary computations that are not in Table 11; so far only `frame_fingerprint.py` (auxiliary). |
| `tests/` | Small cases, run with `pytest`. |

## Inputs and outputs

The programs read their default inputs from the release asset `computation-outputs.tar.gz`, unpacked in the
directory `$RANK32_ASSETS` (by default `assets/` under the current directory): the outputs of computation C<n> are
in `computation-outputs/C<n>/`. This asset is not attached to a release yet; it comes with a later version. Until
then, give the inputs with the options of each program. The programs write their outputs, by default, to
`results/C<n>/` under the current directory; an auxiliary computation that is not in Table 11 writes to `results/`.

Run the programs with the Python of the environment below, from the repository root. For example,

```
python code/produce/element_types/prime_order_type_masses.py test
```

checks the mass formulas of the types in about 10 seconds; `--help` lists the modes and options of each program.

## Notation

The programs write the type p-(z, f)-s of an element of odd prime order (§3.3 of the first paper) as p-(z, d)-s:
d is the rank f of the fixed lattice.

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
| `RANK32_ASSETS` | where the release assets are unpacked; the programs read their default inputs from it (`rank32.io.asset`) | `assets/` under the current directory |
| `RANK32_TMP` | where `rank32.canonical` writes its temporary graph files | the system's temporary directory |
| `FRAME_PL` | the C part of `rank32.frame_fingerprint` splits the lines by their top PL bits; its table has 2^(32-PL) entries | `12` |
| `FRAME_SCRATCH` | where `rank32.frame_fingerprint` writes the binary input of its C part | the system's temporary directory |
| `CC` | the C compiler of `rank32.frame_fingerprint` | `gcc` |
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
