# data/

English | [日本語](README.ja.md)

`classified-lattices.jsonl` has one JSON object per line, one line for each of the 143329 extremal even
unimodular lattices L of rank 32 with Aut L ≠ {±1} (Theorem B).

| Field | Meaning |
|---|---|
| `label` | `32.<\|Aut L\|>.<i>`, the label of the database [D]; §6.8 of the paper fixes the numbering |
| `aut` | the order of Aut L |
| `part` | the part of Theorem B: `general` (S_gen), `c2xc2` (S_(C_2^2)), `cyclic-quaternion` (S_CQ) |
| `group.small_group_id` | `[order, number]` in the Small Groups Library of GAP, when identified; `null` otherwise |
| `group.structure` | a structure description of Aut L where one was recorded; `null` for most lattices of S_gen (12466 of the 12546), whose groups are then identified by `small_group_id` alone (12447), by the chief series (9), or not at all (the ten lattices with \|Aut L\| > 10^7, see below) |
| `group.chief_series_orders` | orders of a chief series, for groups identified only in this way |
| `group.presentation_checked` | a presentation of Aut L was checked on its generators and its order confirmed (`true` for 12536 lattices of S_gen); `false` for S_(C_2^2) and S_CQ, whose groups are determined in §§6.4–6.5 of the first paper, and for the ten lattices of S_gen with \|Aut L\| > 10^7, whose presentations are proved in the second paper and stored in `largest-groups.jsonl` |
| `entry_types` | the kinds of entries of L (§6.1), as strings: `involution <list>` for an involution entry, with the list identifier of its pair of components (for example `involution 16-16d8`; see `certificates/README.md`, section Identifiers); `odd prime 3-<z>-<f>-<s>` for an element of order 3 of type 3-(z, f)-s; `odd prime (type not recorded)` (174 lattices of S_gen) for an element of odd prime order whose type the database does not record; and `order 2^a with g^(2^(a-1)) = -1` (a = 2, 3, 4, 5) for a Hermitian entry of order 2^a. For S_CQ only the entry of order 4 is listed; the entries of order 8, 16 and 32 of the 238 lattices of S_CQ with such elements are not recorded |
| | For the ten lattices with \|Aut L\| > 10^7 all four `group` fields are empty (`null`, `false`): their groups are too large for these identifications and are described in §6.3 of the paper and in the second paper |
| `frame` | for S_(C_2^2), the frame shape of the involution σ: its characteristic polynomial written as ∏ (x^k − 1)^(a_k), for example `2^24 1^-16` for (x² − 1)^24 (x − 1)^−16, and `2^16` for (x² − 1)^16 (a factor with exponent 0 is not written) |

The Gram matrices and generators of the automorphism groups are in the release asset `database-v2.tar.gz`. Its
`index.jsonl.gz` is "sorted by label" in the numerical sense: by \|Aut L\|, then by i, both as numbers (so 32.8.2
comes before 32.8.10, and 32.16.1 after all labels 32.8.i).

## largest-groups.jsonl

`largest-groups.jsonl` has one JSON object per line, one line for each of the 99 lattices of the second paper:
the lattices of S_gen with |Aut L| > 2000 or |Aut L| = 512, 1024 or 1536 (the groups of the other lattices are
identified by `IdGroup` in the Small Groups Library). The lines are in the order of Table A1: by |Aut L|, then by the index i of the label.
`verify/verify_p2.py` checks this file. Its Gram matrices are those of the database; its generators are those to which
the words of the presentations refer, and for six lattices (32.4608.2, 32.5760.1, 32.12288.1, 32.23328.1, 32.196608.2,
32.839808.1) the two lists differ: in two of them the database has one matrix more, in the other four one matrix of each list is
not in the other. Both lists generate Aut L.

| Field | Meaning |
|---|---|
| `label` | `32.<\|Aut L\|>.<i>`, the label of the database [D]; the same label as in `classified-lattices.jsonl` |
| `id` | the construction identifier of L: a list entry `<list>:<index>` as in `certificates/README.md` (section Identifiers), or `known:<family>:<name>`, the name of a construction of a previously known lattice (Appendix A of the first paper) from its published description (a code, or a Hermitian or quaternionic structure) |
| `aut` | the order of Aut L (an integer) |
| `gram` | the Gram matrix G of L (32 × 32 integers) |
| `generators_convention` | `columns: x -> g x, g^T G g = G` |
| `generators` | generators g of Aut L in the column convention g^T G g = G (those of the first paper, §7.1) |
| `table_a1` | the columns of Table A1 other than the label and the order: `aut_factor`, `fstar`, `g_mod_fstar`, `chief`, `presentation`, `origin`, as in `paper/p2-lattices-99.json` |
| `presentation` | a presentation of Aut L whose order is proved by the lemmas of §2.3 (below) |
| `profile` | the profile of the second paper (§2.2) as far as it is data: `structure` (orders of F^*(G), of G/F^*(G), of the centre and of the derived group, the abelian invariants, the orders of a chief series and, for each chief factor from the top, its order, `kind` (`abelian` with `p`, `dim`, the order of the image of the action, absolute irreducibility, splitting and the number of classes of complements; or `non-abelian` with `simple` and `copies`), and `split_note` where splitting was not computed); `action_on_L` (for each conjugacy class: order, size, centralizer order, frame shape, trace and the fixed lattice with rank, determinant, discriminant group and minimum, `null` for rank 0; and the lengths of the orbits on the 146880 minimal vectors); `construction` (for an involution σ in a small class: the size of the class, the order of C(σ), whether σ ~ −σ, the invariants and Gram matrices of F = L^σ and K = L^(−σ), the orders of their groups and the images of C(σ); `null` for the two largest groups, which §4.5 of the second paper leaves out (\|Aut L\| > 2 · 10^10); for `known:zeta17:s1-1`, whose only involution is −1, a `note` says so) |

**The presentation.** Its matrices are in the row convention: a matrix S acts on row vectors, v ↦ v S, and
satisfies S G S^T = G; S is the transpose of an automorphism in the column convention. In this convention the
product of group elements is the matrix product, and a word is evaluated as the product of its factors from left to
right. Generators are numbered from 0, except in the letter representation of GAP (below). Fields of every
presentation:

| Field | Meaning |
|---|---|
| `kind` | `pc` (Lemma 2.1), `radical` (Lemma 2.3 with one block stage for Q = G/R, R the solvable radical, followed by the cyclic stages of R) or `chain` (Lemma 2.3 with stages along a chief series) |
| `lemma` | `2.1` or `2.3` |
| `convention` | the row convention, in words |
| `generator_words` | for each matrix g of `generators` of this file (for six lattices not the same set as in the database, see above), in the same order, a word in the generators of the presentation whose product is g^T |

`kind` = `pc`:

| Field | Meaning |
|---|---|
| `generators` | the matrices x_0, …, x_(m−1), a pc sequence along a chief series |
| `relative_orders` | r_0, …, r_(m−1); their product is \|Aut L\| |
| `relations` | `{"kind": "power", "i", "r", "word"}` for x_i^(r_i) = word and `{"kind": "conjugate", "i", "j", "word"}` for x_i^(−1) x_j x_i = word (i < j). A word is an exponent vector (e_0, …, e_(m−1)) meaning x_0^(e_0) x_1^(e_1) ⋯ x_(m−1)^(e_(m−1)), with e_k = 0 for k ≤ i |
| `generator_words` | exponent vectors as in `relations` |

`kind` = `radical`:

| Field | Meaning |
|---|---|
| `radical_generators`, `radical_relative_orders`, `radical_relations` | a pc presentation of R, with the fields of `kind` = `pc`; the product of the relative orders is \|R\| |
| `lifts` | the matrices y_1, …, y_n: lifts to G of the generators of a presentation of Q |
| `lift_conjugates` | `{"lift": a, "j", "word"}` for y^(−1) x_j y = word, where y is the lift with index a (numbered from 0) and the word is an exponent vector in the x |
| `quotient_relators` | the relators of the presentation of Q in GAP's letter representation: the letter a > 0 is y_a and −a is y_a^(−1) (numbered from 1) |
| `quotient_relator_values` | for each relator, the exponent vector in the x of its value on the lifts |
| `quotient_order` | \|Q\| = \|Aut L\| / \|R\| |
| `quotient_presentation_order` | the order of the group defined by `quotient_relators`, with `method`: coset enumeration on the trivial subgroup or, for `12-20-c63:236`, Lemma 2.4 with X_0 = {y} and R_0 = {y^r} (`lemma_2_4`: the generator y, the relator y^r, and the index of ⟨y⟩ counted by coset enumeration). This order is an input (§2.5) |
| `generator_words` | `{"radical": e, "quotient": w}`: g^T is the product of the exponent vector e in the x and the letter word w in the y, in this order |

`kind` = `chain`:

| Field | Meaning |
|---|---|
| `generators` | all matrices, stage after stage |
| `stages` | `{"kind": "cyclic", "r", "generators": [a]}`, or `{"kind": "block", "k", "T_order", "components"}` with one component for each of the k minimal normal subgroups of a non-abelian chief factor T^k. A component has its `generators` (indices), its `relators` (letter representation, numbered from 1 within the component) and `presentation_order`: the order of the group defined by these relators, with `method`, either coset enumeration on the trivial subgroup or the chain of Lemma 2.4 in `lemma_2_4_chain` (at each step the generator subset, the number of relators in its letters, the next subset and the index counted by coset enumeration; at the last step the order counted by coset enumeration). These orders are inputs (§2.5) |
| `relations` | `{"kind": "power", "stage", "a", "r", "word"}` for x_a^r = word; `{"kind": "relator", "stage", "component", "relator", "word"}` for a relator of a component, evaluated on its generators, = word; `{"kind": "commutator", "stage", "a", "b", "word"}` for x_a^(−1) x_b^(−1) x_a x_b = word (a and b in different components); `{"kind": "conjugate", "stage", "a", "b", "word"}` for x_a^(−1) x_b x_a = word (b in a later stage). A word is a list of pairs `[generator, exponent]`, the product of the powers from left to right; all its generators are in the stages after `stage` |
| `order_bound` | the bound of Lemma 2.3: the product of r over the cyclic stages and of \|T\|^k over the block stages; it equals \|Aut L\| |
| `generator_words` | lists of pairs as in `relations` |

The right-hand sides of the relations and the words for the generators were found with GAP (polycyclic exponents,
homomorphisms onto the finitely presented groups, and sifting through the stages). `verify/verify_p2.py` evaluates
them with integer matrices, so an error in them would be detected; how they were found does not affect the
conclusions.

## Reading the data

With the Python standard library only, from the repository root: read one lattice of `largest-groups.jsonl`
and check that its Gram matrix is even and symmetric and that every generator is an isometry.

```python
import json
rows = {r['label']: r for r in map(json.loads, open('data/largest-groups.jsonl'))}
r = rows['32.512.1']
G, gens = r['gram'], r['generators']            # 32 x 32 integers; generators act on columns: x -> g x
n = len(G)
mul = lambda A, B: [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
tr = lambda A: [list(c) for c in zip(*A)]
assert all(G[i][j] == G[j][i] for i in range(n) for j in range(n)) and all(G[i][i] % 2 == 0 for i in range(n))
assert all(mul(mul(tr(g), G), g) == G for g in gens)   # each generator is an isometry of L
print(r['aut'], len(gens))
```

## Fields of the release assets

The assets `database-v2.tar.gz`, `general-part-records.tar.gz` and `general-part-groups.tar` contain their own README.
For the others:

`type-8-24-classes.tar` has four files `classes-0.jsonl.gz`, ..., `classes-3.jsonl.gz` with the 376306 orbits of
good subspaces W for P = A_1^24 (Appendix B.3 of the first paper); orbit i is in the file `classes-<i mod 4>`.

| Field | Content |
|---|---|
| `i` | the index of the orbit, 0, ..., 376305 |
| `stab` | \|Stab(W)\| in Aut P |
| `gram` | the Gram matrix of K_W (rank 24, even, determinant 256, minimum 4), LLL-reduced |
| `aut` | \|Aut K_W\| |
| `index` | [N : N'], where N = √2 K_W^# and N' is spanned by the 720 vectors of norm 2 of K_W^# |
| `method`, `cert` | how K_W was identified: `graph` (the canonical form of the graph of the 720 vectors; `cert` is its hash), `orbit` (the same with the orbit of N / N'; `cert` adds a hash of its least element), or `qfauto` (`cert` is a hash of the invariants, and isometries were decided with `qfisom`) |
| `orbit` | the index of Aut K_W in the automorphism group of that graph (0 when not computed, for `qfauto`) |
| `census` | for the 270 maximal totally singular subspaces T' of K_W^# / K_W, the number of T' by the number of roots of K_W + T' |
| `np` | for the 135 non-zero singular classes p of K_W^# / K_W, the number of classes by the number of vectors of norm 2 in p |

`gaussian-classes.jsonl.gz` has one line for each of the 76742 rootless Hermitian Z[i]-lattices of Proposition 4.6:

| Field | Content |
|---|---|
| `key` | the key of §7.2 of the first paper: `k` and 24 hexadecimal digits (the canonical Gram matrix of a discrete pair) or `D` and 16 hexadecimal digits (the invariant D) |
| `D` | the invariant D (16 hexadecimal digits) |
| `gram`, `W` | the Gram matrix G of the trace lattice and the matrix W of multiplication by i, in the row convention: W G W^T = G, W^2 = −1 |
| `U` | the order of the unitary group |
| `antilinear` | whether the lattice is isometric to its complex conjugate |
| `f`, `w` | f = 1 if `antilinear`, otherwise 2, and the weight w = f / \|U\| of Proposition 4.6 |

`gaussian-walk-grams-1.jsonl.gz` and `gaussian-walk-grams-2.jsonl.gz` have the fields `key`, `gram` and `W` as above
for the classes met by the two neighbour walks, and `id`: `k11-<n>` for the first walk and `k12-<n>` for the second
(11 and 12 are the seeds of the walks), where n is the number of rootless lattices the walk had visited when it met the
class.
