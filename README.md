# The extremal even unimodular lattices of rank 32

English | [日本語](README.ja.md)

This repository accompanies two papers, published together in each version:

- *The extremal even unimodular lattices of rank 32: their number, and the classification of those with
  automorphisms other than ±1* (`paper/extremal-lattices-rank-32.pdf`, [P1] in the second paper);
- *Automorphism groups of extremal even unimodular lattices of rank 32* (`paper/automorphism-groups-rank-32.pdf`),
  on the 99 lattices with the largest groups.

**Results.** There are exactly 11,044,268 isometry classes of extremal even unimodular lattices of rank 32
(Theorem A). Exactly 143,329 of them have automorphisms other than ±1; the paper classifies these lattices and
determines the order of each automorphism group (Theorem B). The proof is computer-assisted. Table 7 of the
paper lists the statements on which it rests, and Table 10 lists what it assumes about the software.

The second paper describes the groups of the 99 lattices whose groups have order larger than 2,000 or equal to 512, 1,024 or 1,536 (the groups of the other lattices are identified by `IdGroup` in the Small Groups Library):
generalized Fitting subgroup, chief series, a presentation of proved order, conjugacy classes and orbits on the
minimal vectors; 27 of the 99 are in the literature, and the other 72 were not found in the literature we read.

Japanese translations of the two papers and of the note are in `paper/extremal-lattices-rank-32.ja.pdf`, `paper/automorphism-groups-rank-32.ja.pdf` and `paper/assignment-count.ja.pdf`; the numbering of theorems, formulas, tables and computations is the same in both languages.

## Contents

| Path | Content |
|---|---|
| `paper/` | The two papers (Typst sources and PDFs) with their computation ledgers (`computations.typ`, Table 11 of the first; `computations-p2.typ`, Table 11 of the second), the data of Table A1 of the second (`p2-lattices-99.json`), and an auxiliary note on an earlier organization of the count, which shares the component lists (`assignment-count.pdf`, [Aux] in the first paper) |
| `certificates/` | The certificates read by the checks, one file per statement or computation; see `certificates/README.md` |
| `data/classified-lattices.jsonl` | The 143,329 lattices with automorphisms other than ±1: label, order and (where identified) structure of the automorphism group, and the part of Theorem B to which the lattice belongs; see `data/README.md` |
| `data/largest-groups.jsonl` | The 99 lattices of the second paper: label, Gram matrix, generators of the automorphism group, the presentation and the columns of its Table A1; see `data/README.md` |
| `verify/` | `verify.py`, which checks the certificates and data against the first paper, and `verify_p2.py`, which checks the data of the second in integers; see `verify/README.md` |
| `LICENSE` | MIT for the code in `verify/`; CC BY 4.0 for text and data |
| `MANIFEST.json` | SHA-256 digests of all files of this version |

Release assets (Table 12 of the paper) hold the larger data and the research record:

| Asset | Content |
|---|---|
| `research-record.tar.gz` | The programs that produced the lists and certificates, unchanged, with an index by computation and the verification notes ([Rec] in the paper) |
| `research-runs.tar.gz` | The recorded runs of those programs |
| `gaussian-classes.jsonl.gz` | The 76,742 rootless Hermitian Z[i]-lattices of Proposition 4.6 |
| `gaussian-walk-grams-1.jsonl.gz`, `gaussian-walk-grams-2.jsonl.gz` | Gram matrices met by the two neighbour walks that built the Gaussian list (fields in `data/README.md`) |
| `type-8-24-classes.tar` | The 376,306 orbits of good subspaces W for the Niemeier lattice A_1^24 (Appendix B.3 of the first paper), each with a Gram matrix of K_W, \|Aut K_W\| and the order of the stabilizer of W |
| `general-part-groups.tar` | Generators, presentations, conjugacy classes and group identifications for the 12,536 lattices of S_gen with \|Aut L\| <= 10^7, by label, with a README |
| `general-part-records.tar.gz` | The 17,736 full-group records of S_gen (§6.3 of the first paper): for each, the label, \|Aut L\|, the conjugacy classes and the rows of entries, with a README |
| `database-v2.tar.gz` | The database [D], version 2: for each of the 143,329 lattices of S its Gram matrix and generators of its automorphism group, with a README, an index by label and a manifest; `verify/verify_database.py` checks it against this repository and looks up a lattice by its label |

## Checking

```
python3 verify/verify.py
python3 verify/verify_p2.py
```

Both need only Python 3 with its standard library (tested with Python 3.12 and 3.14) and the files of this repository, and exit with status 0 when every
check passes. `verify.py` checks the finite certificates, re-adds the completeness identities from their terms, replays
the bookkeeping of §6 of the first paper from the saved rows, and checks the agreement of the data with that paper;
it does not repeat the group, list and mass computations (see `verify/README.md` and §7.4 of the first paper).
`python3 verify/verify_database.py <path to database-v2.tar.gz>` checks the release asset of the database: its
manifest and index, every label against `data/classified-lattices.jsonl`, every Gram matrix (even, positive
definite, determinant 1), every generator (g^T G g = G) and, for |Aut L| <= 10^6, that the generators generate a
group of order |Aut L|, in about half an hour; with `--label 32.<|Aut L|>.<i>`
it prints the record of one lattice.
`verify_p2.py` checks, in integers, that the given generators are automorphisms of the 99 lattices, that the relations
of each presentation hold for them and that the order of each presentation is |Aut L|; the orders obtained by
coset enumeration are taken as stated (§2.5 of the second paper).

## Reproducing

There are three levels, from cheap to expensive.

1. **The repository.** `python3 verify/verify.py` and `python3 verify/verify_p2.py`, as above (about half a minute and one minute). `verify/README.md` lists every check and what
   it takes as input.
2. **One computation.** Table 11 of the paper (`paper/computations.typ`) names, for each computation C1, C2, ..., its
   certificate here and its place `C<n>/` in the release asset `research-record.tar.gz`. `INDEX.md` there lists the
   programs of each computation (C1, C2, ... of the first paper and A1, ..., A16 of the second, the latter under
   `A<n>/`) and, where one was written, its rerun note `REPRODUCE.md` with the command, the
   inputs and their digests, the outputs and what a rerun must reproduce. The programs are kept as they were run,
   with their research paths; they need SageMath, PARI/GP, GAP and bliss, and some ran for hours to days (Table 11).
3. **The data.** The release assets hold the lists and records the certificates were exported from; their sizes and
   SHA-256 digests are in `MANIFEST.json` (field `release_assets`) and in Table 12 of the first paper. The certificates name the computations they come from in their field
   `computations`.

## Typesetting

The PDFs in `paper/` were typeset from the sources next to them with Typst 0.14.2, for example
`typst compile paper/extremal-lattices-rank-32.typ`. The PDFs use the fonts New Computer Modern, New Computer Modern
Math, Libertinus Serif and DejaVu Sans Mono, which come with Typst, and DejaVu Sans, which must be installed. With another version of Typst or other fonts the text is the same but line and page breaks may differ.
The Japanese versions (`*.ja.typ`) also use IPAex Mincho and IPAex Gothic, which must be installed.

## Versions

Versions are tagged `vYYYY-MM-DD.N`: the date (Japan Standard Time) of the research commit from which the version was
built, which is also the `source_date_epoch` of `MANIFEST.json`, and the number of the version on that day. A later
version may correct or extend an earlier one; the paper states its version on the title page.
`MANIFEST.json` also records the commit of the research repository from which a version was built; that repository is
not public, and the commit serves only to identify the build.

## Authorship and licence

roiban1344 directed the work and is responsible for publication. AI systems, including Claude (Anthropic) and OpenAI
Codex, contributed to the research, code, checking, and writing.

Code is distributed under the MIT licence; text and data under CC BY 4.0. See [LICENSE](LICENSE).
