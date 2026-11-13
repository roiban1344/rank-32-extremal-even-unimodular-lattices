# verify/

English | [日本語](README.ja.md)

```
python3 verify/verify.py [report.json]
```

Run from the repository root. Only the Python standard library is used. The exit status is 0 when every check
passes; the optional argument writes a JSON report.

### Which statement is checked where

The rows are those of Table 7 of the first paper. "Input" means that the check takes the value from a certificate
and does not recompute it.

| Table 7 row | statement | checked by | not recomputed here |
|---|---|---|---|
| 1 | King's rootless mass m | D (the value printed in the paper equals the certificate) | the mass itself |
| 2 | possible types (Prop. 3.3) | C (the six LP certificates, rows rebuilt from the theta bases) | the bases of modular forms |
| 3 | involution components | I (identity (B.1), Table B3) | the enumeration of the components |
| 4 | components of odd prime order | I (Table B4 against the masses, the steps of (B.2)); H (App. D.4) | masses of genera, rootless masses |
| 5 | gluing classes | I (identity (4.1) for the 83 pairs and the 19 types; \|O(q_F)\| from Lemma 2.1) | the stabilizers \|Stab(φ)\| |
| 6 | Gaussian rootless mass M_2 | E (every term of (C.2) from the root data) | the local 2-adic statistics |
| 7 | Gaussian list (Prop. 4.6) | I (Σ f/\|U\| = M_2; the counts by method) | \|U\| and the distinctness of the classes |
| 8, 9 | keys; finite groups | — (written proofs) | — |
| 10 | the lattices of S_gen (§6.3, Lemma 6.3) | F, G, I | that the groups are the full groups |
| 11 | subtraction (Lemma 6.6) | G | the identifiers and keys of the rows |
| 12 | identity (5.1) | D | — |

### The sections of verify.py

| Section | Check |
|---|---|
| A | every file listed in `MANIFEST.json` has the recorded SHA-256 digest |
| B | every computation cited in the paper appears in its ledger (Table 11) |
| C | the six rational certificates of §B.1.1: multipliers y ≥ 0 with yA = 0 and yb = −1; and the inequalities themselves, rebuilt from the stored bases of theta series and dual theta series (normalizations, non-negative coefficients, gluing equations), with the numbers N(j) of classes with q = 2j/p recomputed from (p, d, s) by Lemma 2.1(i). The bases, that is the spaces of modular forms, are inputs |
| D | identity (5.1): N = 2m + Σ (1 − 2/\|Aut L\|) is an integer for King's mass m and the classified lattices; the totals stated in Theorems A and B; Table 8; Table 9 read from the paper and compared with the data cell by cell |
| E | every embedding average A_16(T) recomputed from the root data and the saved local 2-adic statistics by `check_m2_terms.py` (the statistics, the exceptional coefficients and the closed formulas of §C.3 are inputs); every term of (C.2) (`gaussian-rootless-mass.json`), and the counts of `root-coefficients.json` against it: coefficient as a product of irreducible coefficients, \|U(T)\| from Proposition C.1, term = c·A/\|U\|; their sum is M_2 of Theorem 4.5; the coefficients of A_n, D_n, G_m from the closed formulas of §C.3 |
| F | in integers: the generators of the two lattices with the largest groups preserve their Gram matrices; the 106 isometries of §6.3 satisfy TᵀGT = G′ and det T = ±1 |
| G | the counts of §6 against the classified lattices: 12,546 classes of S_gen; for Lemma 6.3, the saved counts (12,519 classes with exact identifiers) and the 27 saved profiles, pairwise different in each bucket; the 67,284 lattices of S_(C_2^2); for the subtraction of Lemma 6.6, the saved totals of the certificate. Then, from the rows of `general-part-witnesses.json`: the classes after each joining step of §6.3 (13,396 by identifiers, 12,585 with the inverse identifiers, 12,546 with the isometries); the 17,736 full-group records, their identifiers, the inverse edges and the 106 isometries join into exactly 12,546 classes, each with one lattice of the general part and the same \|Aut L\| as its records, and the entries and isometries of `general-part-identification.json` lie in the class of their lattice; for Lemma 6.3, the exact types (types with no row lacking an identifier) are read off the rows, each of 12,519 classes has a witness row of an exact type, and the other 27 classes are exactly those of the 27 saved profiles, pairwise different; for Lemma 6.6, the subtraction as set differences of individual entries: the extremal involution entries with \|Stab\| = 4 (and the entries of type (8, 24) with \|Aut K\| = 2) minus the self-inverse entries and the entries of S_gen are exactly the oriented entries of the 67,284 lattices of S_(C_2^2), and the 76,742 Gaussian keys minus the 11,545 keys of S_gen are exactly the keys of the 63,499 lattices of S_CQ (one key with f = 2 and \|U\| = \|Aut L\| for a cyclic group, three keys with f = 1 and \|U\| = 4, 4, \|Aut L\|/2 for a quaternion group). For the 849 quaternion lattices, in integers from `quaternion-triples.json`: each of the 2,547 Gram matrices G is symmetric, even, positive definite (leading minors by Bareiss elimination) and has det G = 1; the matrix W of multiplication by i satisfies W G Wᵀ = G and W² = −1 (W acts on row vectors); in each triple two isometries T with TᵀG_first T = G_second and det T = ±1 join its three classes (1,698 isometries, every class in exactly one triple); the three keys of each triple are exactly the keys of its label in `cq_lattices` of `general-part-witnesses.json`, with \|U\| and f as in `gaussian-list-weights.json`. Inputs, not checked here: that the groups of the records are the full groups, that their conjugacy-class rows are complete (types, centralizer orders), that the identifiers and Gaussian keys of the rows are correct, the inverse and self-inverse identifiers, and the positive isometries joining the entries (for the 106 of §6.3 and the 1,698 of the quaternion triples only T^T G T = G' and det T = ±1 are checked, in F and G); that the key and \|U\| of a Gaussian class belong to its (G, W), and the distinctness of the Gaussian classes, by which different triples are different lattices (Proposition 6.5) |
| H | the finite checks rerun to their certificates: the 9-torsion reduction (`check_exponent9.py`) and the orbits of the elementary quadratic 2-forms of Appendix D.4.6 (`check_quadratic2_orbits.py`, all 230 forms of dimension 1 to 20) |
| I | the completeness identities, added again from their terms: \|O(q_F)\| of every pair of components from the closed formulas of Lemma 2.1 (from the rank of F and s, δ of its Gram matrix, or from the type p-(z,f)-s); (4.1) for the 83 pairs of involution components (with s, δ and 2-elementarity of the components from their Gram matrices) and for the pairs of the 19 types of odd prime order (term by term for 258 pairs, from the numbers of decompositions for the 4 pairs of 5-(6,8)-6, and from the recorded mass for the 2 pairs of 3-(6,20)-6 without extremal gluings); the component lists of Table B4 against the masses of their genera or the rootless masses (for the lists completed by (B.2) the rootless mass is a check), the steps of (B.2) used for Table B4, term by term from \|U(Z_+)\| and g(Z_+), and that these steps join the lists of consecutive levels back to a list (z, z) of rootless mass; identity (B.1) for type (8, 24) (for A_2^12 from the terms of its 106,624 orbits of good 3-dimensional subspaces), and that the four lattices K whose \|Aut K\| was computed with qfauto have orders that no other K has; Σ f/\|U\| = M_2 over the 76,742 Gaussian classes (Proposition 4.6), and the number of classes whose weight came from each method (canonical form, backtrack search, bliss) as stated in its proof, the canonical-form classes being exactly those with a canonical key, with \|U\| = 4 and f = 2; the counts of Tables 5, 6, B3 and B4; Lemma 6.3(ii): for every list of §4.1–§4.2, the rows of the 12,546 lattices of S_gen (one record each) add up to the extremal mass of the list less its trivial classes assigned to S_(C_2^2), none of which is the identifier of a row, the rows with identifiers give distinct classes with \|C\| = \|Stab\|, and the pairs saved by their mass only are the non-extremal pairs of 3-(6,20)-6; the list entries of S_(C_2^2) and S_gen are extremal classes of their lists, and every identifier of order 3 in the full-group records of S_gen is an extremal double coset of its list with \|Stab(φ)\| equal to the centralizer order of its row |

**What is not checked here.** The completeness identities are added again from their terms (section I), but the
terms themselves are not recomputed: the automorphism groups and stabilizers, the numbers of decompositions and of good lines g(Z_+), the masses of genera and rootless masses,
and the isometry decisions that make the listed classes distinct. They rest on the arguments of the paper (Table 7) and on the assumptions about
the software (Table 10). The programs that computed them are in the research record (release asset
`research-record.tar.gz`).

## verify_p2.py (the second paper)

`python3 verify/verify_p2.py [report.json]` checks `data/largest-groups.jsonl`, the data of the 99 lattices of the
second paper, in exact integer arithmetic (standard Python only; about one minute):

| Section | What is checked |
|---|---|
| A | the 99 labels: distinct, of the form `32.<\|Aut L\|>.<i>`, in the order of Table A1, each with the same \|Aut L\| as its row of `data/classified-lattices.jsonl` |
| B | each Gram matrix: integer, symmetric, even, det 1, positive definite |
| C | each generator g of Aut L: an integer matrix with g^T G g = G and det g = ±1 |
| D | each presentation (§2.3 of the second paper): its matrices preserve G; every relation required by its lemma is present, with right-hand side in the later stages, and holds on the integer matrices; every generator of Aut L is the product of its stored word, so the matrices of the presentation generate Aut L; the bound of the lemma equals \|Aut L\| |
| E | the columns of Table A1 in the data equal `paper/p2-lattices-99.json`; the factorization of \|Aut L\| and the lemma of each presentation agree with the table |
| F | the profiles: for every conjugacy class, size × centralizer order = \|Aut L\|, and its frame shape has degree 32, multiplicity of the eigenvalue 1 equal to the rank of the fixed lattice and a_1 equal to the trace; the class sizes add up to \|Aut L\| with one identity class; the orbits on the minimal vectors add up to 146,880. The classes themselves (computed with GAP) are inputs |

Not recomputed (the assumptions of §2.5 of the second paper): the orders \|Aut L\| themselves, and the orders of the
finitely presented quotients and simple sections, counted by coset enumeration directly or along the chain of
Lemma 2.4; the script prints these orders as inputs and checks the chains of Lemma 2.4 for consistency.

## verify_database.py (the release asset of the database)

`python3 verify/verify_database.py database-v2.tar.gz [report.json]` checks the release asset of the database [D],
version 2, against this repository (standard Python only; about half an hour and 1.2 GB of memory). Without the asset it prints
SKIP and exits with status 0. With `--label 32.<|Aut L|>.<i>` it prints the record of one lattice.

| Section | What is checked |
|---|---|
| A | every file of rows has the size, SHA-256 digest and number of rows of the manifest of the archive; 143,329 rows in all |
| B | the index has every label of `data/classified-lattices.jsonl` exactly once, in label order, and every row is where the index says, in a file of its part of Theorem B, with the order \|Aut L\| of its label |
| C | every Gram matrix is even, positive definite and of determinant 1 |
| D | every generator is an integer matrix g with g^T G g = G |
| E | for the 143,316 lattices with \|Aut L\| <= 10^6: the exact orbit of an integer vector under the group generated by the stored generators has \|Aut L\| elements; by D that group is a subgroup of Aut L, so it is all of Aut L, given the order of the first paper |

Not checked: that the lattices have no roots, that the generators of the 13 lattices with \|Aut L\| > 10^6 generate all
of Aut L, and that different labels are different lattices; these are the statements of the first paper and rest on
its proofs.
