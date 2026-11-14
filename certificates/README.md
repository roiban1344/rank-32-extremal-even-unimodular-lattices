# certificates/

English | [日本語](README.ja.md)

Each file is a JSON object with the fields `certificate` (title),
`statement` (what it certifies), `computations` (rows of Table 11 of the paper) and `scope` (what it
does and does not show), followed by the data; in the larger files each field, and each element of a list,
is on a line of its own. `verify/verify.py` compares every file with its digest in
`MANIFEST.json` and checks the content of every certificate against the paper, within the scope stated in the
file. For most certificates the check is the recomputation of an identity from its terms, the replay of the
bookkeeping from saved rows (`general-part-witnesses.json` for Lemmas 6.3 and 6.6), the check of integer
matrices (isometries and Gram matrices, as in `quaternion-triples.json`) or the rerun of a finite check.
`verify/README.md` says which, and what each check takes as input.

| File | Certificate | Computations |
|---|---|---|
| `c2xc2-part.json` | The 67284 lattices of S_(C_2^2) | C26, C28 |
| `exceptional-full-groups.json` | Full automorphism groups of the two known lattices with the largest groups | C22 |
| `exponent-9-reduction.json` | The two Z/9 blocks of the 9-torsion reduction | — |
| `gaussian-list-weights.json` | Weights of the 76742 rootless Hermitian Z[i]-lattices of rank 16 | C12 |
| `gaussian-rootless-mass.json` | The rootless mass M_2 of Hermitian Z[i]-lattices of rank 16 | C11, C19 |
| `general-part-distinctness.json` | Distinctness of the 12546 lattices of S_gen | C24, C25 |
| `general-part-identification.json` | Identification of the 12546 lattices of S_gen | C1, C23, C24 |
| `general-part-witnesses.json` | Row-level witnesses for Lemmas 6.3 and 6.6 | C24, C25, C29 |
| `involution-gluing-masses.json` | Gluing classes of the involution types (15, 17), (14, 18), (12, 20) and (16, 16) | C4, C5, C7 |
| `king-rootless-mass.json` | King's rootless mass in rank 32 | C14, C15 |
| `odd-prime-masses.json` | Components and gluing classes of the 19 types of odd prime order | C2, C8, C9, C16, C17 |
| `quadratic-2-form-orbits.json` | Orbits of elementary quadratic 2-forms used in the 2-adic counts | C18 |
| `quaternion-triples.json` | The quaternion lattices of S_CQ | C27, C28 |
| `root-coefficients.json` | Closed formulas for the inclusion-exclusion coefficients over Z[i] | C20 |
| `subtraction-records.json` | Records used by the subtraction Lemma 6.6 | C29 |
| `type-8-24-masses.json` | The 14482 lattices K of the involution type (8, 24) | C6 |
| `type-exclusion-lp.json` | Rational infeasibility certificates for six types of odd prime order | C3 |

## Identifiers

Lattices are named by their database labels `32.<|Aut L|>.<i>` (§6.8 of the paper).
Entries of the lists of §4 are named by list identifiers `<list>:<index>`, where `<list>` names the type
and the pair of components: `15-17`, `14-18<K>`, `12-20-c<K>`, `16-16-s<s>d<delta><F,K>`, `16-16s14` (s = 14,
delta = 1) and `16-16d8` (s = 16, delta = 1) for the involution types of Tables 5 and 6, `8-24` for type (8, 24),
and `<z>-<f>-<s>` for p = 3 and `<p>-<z>-<f>-<s>` for p >= 5 for the types of odd prime order of Table 4.
For an involution list, `<index>` is the position of the class in the list of its pair in
`involution-gluing-masses.json`, and for `8-24` the position of K in `type-8-24-masses.json`.
For p = 3 an entry is `<z>-<f>-<s>:<i>-<j>/<k>`: the k-th double coset of the pair of the i-th lattice F of
minimum at least 4 and the j-th lattice Z of their genera in `odd-prime-masses.json` (`genera_F`, in list
order, counting only the classes with `min4`; `genera_Z`), which is the k-th entry of `stabilizers` of that pair.
Exceptions, used only as names of the constructions behind full-group records: `13-6-5:<j>` is the pair with the
j-th lattice Z of 3-(13,6)-5 (F is unique and there is one double coset), `6-20-6:<h>/<k>` names a pair of
3-(6,20)-6 by the first 12 hexadecimal digits of a hash of F and the coset number, and `<p>-<z>-<f>-<s>:<n>` is
the n-th extremal pair of its type in the research record. `record:<i>` is a full-group record of a lattice
built otherwise (the previously known lattices of Appendix A).

The rows of a full-group record (`general-part-witnesses.json`) have a type: `inv:<list>` for an involution,
`3-<z>-<f>-<s>` for order 3, and `odd:p=<p>:f=<f>` for p >= 5. The rows for p >= 5 do not record s, and none
of them carries an identifier, so the coarser name never makes a type exact (Lemma 6.3). Gaussian classes are named
by their keys (§7.2): `k` followed by 24 hexadecimal digits is a hash of the canonical Gram matrix, `D` followed
by 16 hexadecimal digits the invariant D, used where the selection of lines of §7.2 does not give a canonical
Gram matrix (the selected lines do not generate L, or their colours do not order them).
The release assets use the same identifiers.

## c2xc2-part.json

Trivial involution entries not met by S_gen, with [phi] and [phi^-1] identified for F = K (Lemma 6.6(i)); each gives one lattice with Aut L = {+-1, +-sigma}.

*Scope.* The list entries and their lattices. Completeness of the gluing lists and correctness of the identifiers are inputs (§4, Lemma 6.6).

## exceptional-full-groups.json

|Aut L| = 2 |orbit of a minimal line| |stabilizer of a minimal vector|, computed from the Gram matrix alone by the algorithm of §7.1.

*Scope.* Gram matrices, the computed orders and generators (column convention: g^T G g = G). The generators are checked in integers; that they generate the full group is the output of the algorithm of §7.1.

## exponent-9-reduction.json

All 48 maps of GL_2(F_3) are checked; 8 graphs are admissible (§D.4.6).

*Scope.* Finite check by integer modular arithmetic; rerun by verify/check_exponent9.py.

## gaussian-list-weights.json

Proposition 4.6: sum over the list of f/|U| equals M_2 of Theorem 4.5, where f = 2 if the lattice is not isometric to its complex conjugate and f = 1 otherwise; hence the list is complete (Lemma 7.1).

*Scope.* For each class: its key (canonical-form key or invariant key of §7.2), |U| and f. The orders |U|, the flags f and the distinctness of the classes are inputs here (Proposition 4.6 and Table 10); the Gram matrices are in the release asset of Gaussian classes under the same keys. weight_methods[i] names the methods that gave (|U|, f) of classes[i] with this value (canon: the canonical form; search: the backtrack search; bliss: Aut L from bliss and the centralizer from GAP; joined by +); the counts by method are those of the proof of Proposition 4.6.

## gaussian-rootless-mass.json

M_2 = sum over root types T of c(T) A_16(T) / |U(T)| (C.2), with the coefficients c of Proposition C.3 and the embedding averages A_16(T) of Theorem C.4.

*Scope.* All 12934 terms with their coefficients, group orders and embedding averages, and the local 2-adic statistics they use. The local statistics are inputs here; their computation is the finite transition of Appendix C.6.

## general-part-distinctness.json

Lemma 6.3: an identification class with an exact identifier on a row of an exact type is distinct from every other class; the remaining 27 classes have pairwise different minimum-shell profiles within each full-group signature bucket.

*Scope.* Counts of exact identifiers by type in the saved full-group records, the 27 profiles and the signature buckets. Conditional on full groups, complete conjugacy-class rows and correct identifier lookup (Table 10).

## general-part-identification.json

Every constructed lattice of S_gen belongs to exactly one class below; classes are joined by shared exact identifiers, inverse gluing identifiers and the explicit integer isometries listed here (§6.3).

*Scope.* Positive identifications. Distinctness of the classes is the content of general-part-distinctness.json and Lemma 6.3.

## general-part-witnesses.json

Lemma 6.3: the saved full-group records of S_gen, their exact identifiers, the inverse identifiers and the integer isometries join into 12546 classes, one per lattice of the general part; 12519 classes carry a row of an exact type with an identifier (the witness below) and the other 27 are the classes with saved profiles. Lemma 6.6: the extremal involution entries with |Stab| = 4 (and the entries of type (8, 24) with |Aut K| = 2), minus the self-inverse entries and the entries met by S_gen, are exactly the oriented entries of S_(C_2^2); the Gaussian keys of Proposition 4.6 minus those met by S_gen are exactly the keys of the lattices of S_CQ.

*Scope.* For each successful full-group record of S_gen (record:<i> names the i-th record where the lattice is not a list entry): its list entry, |Aut L|, the rows of its entries (§6.1) as [row, kind, type, identifier or null, centralizer order, self-conjugate]: one row for each involution entry {[sigma], [-sigma]} with sigma != +-1 (self-conjugate: sigma is conjugate to -sigma, so that the entry is one class) and one for each conjugacy class of elements of odd prime order; and its Hermitian rows of order 4 as [row, key, |U|, f]. Further the inverse edges [phi] -- [phi^-1] (the 8716 inverse identifiers of general-part-identification.json and the reverses of 6 of them, 5005 unordered pairs), the isometry edges of general-part-identification.json, the database source of each lattice of the general part, the witness row of each class, the rows completed by high-shell keys, the self-inverse entries and, for each lattice of S_CQ, its Gaussian keys. verify/verify.py replays the bookkeeping from these rows and the other certificates, including the identity of Lemma 6.3(ii): for each list, the rows of the lattices of S_gen, one record each, weighted 1/|C| (2/|C| for an involution entry of a list with F = K that is not self-conjugate), add up to the extremal mass of the list less its trivial classes assigned to S_(C_2^2). Conditional on full groups, complete conjugacy-class rows, correct identifiers and keys, correct self-inverse and inverse identifiers and the positive isometries (Table 10); none of these is recomputed here.

## involution-gluing-masses.json

For each of the 83 pairs of components (F, K), the listed double cosets [phi] satisfy identity (4.1): sum 1/|Stab(phi)| = |O(q_F)| / (|Aut F| |Aut K|), so the list is complete (Proposition 4.1, §4.1.3).

*Scope.* For each pair: the components (Gram matrices, |Aut|, s, delta), |O(q_F)|, the kernels of Aut F -> O(q_F) and Aut K -> O(q_K), and for each class, in list order, |Stab(phi)| and whether L_phi is extremal. The list entry <list>:<i> is the i-th class of the pair. The orders |Aut F|, |Aut K| and |Stab(phi)| are inputs (Table 10); |O(q_F)| is the value of the recursion of §4.1.3, an upper bound, which equality in (4.1) shows to be exact. The identity itself is checked here.

## king-rootless-mass.json

m = sum over extremal even unimodular lattices L of rank 32 of 1/|Aut L|.

*Scope.* The exact value of [Kin03, Table 1, PDF p. 6], row of the empty root system (Weyl group of order 1). Adopted from the source. This certificate holds the value only: the inclusion-exclusion of [Kin03] was recomputed with our implementation (Table 10 of the paper), and verify does not repeat it.

## odd-prime-masses.json

For each type p-(z,f)-s of Table 4 and each pair of components (F, Z), the double cosets [phi] satisfy identity (4.1): sum 1/|Stab(phi)| = |O(q_F)| / (|Aut F| |U(Z)|) (Proposition 4.1, Section 4.2). The component lists of Table B4 are complete: for F, sum 1/|Aut F| over the classes equals the mass of the genus, or for (20,6), p = 3 and (16,4), p = 5 over the rootless classes the rootless mass; for Z, sum 1/|U(Z)| equals the rootless Hermitian mass of the genus, and the steps of identity (B.2) hold (Appendix B.4, Appendix D).

*Scope.* genera_F and genera_Z: for each genus its classes in list order (|Aut F| and minimum >= 4; |U(Z)|), the target mass and how completeness is obtained; genera_Z ends with the upper levels of the (B.2) chains that are not in Table B4 (table_B4 false), so that each chain can be followed to its root (z, z). types: for each pair the indices into the class lists (F, Z), |Aut F|, |U(Z)| and, for each double coset, |Stab(phi)| and extremality; cosets_saved says what is recorded when not all double cosets are listed. b2_steps: the two sides of (B.2), term by term (|U(Z_+)|, g(Z_+)). The group orders, the numbers g(Z_+), the masses of the genera and the rootless masses are inputs (Table 10); the identities are checked here.

## quadratic-2-form-orbits.json

For every diagonal elementary quadratic form over F_2 of dimension 1 to 20, the element orbits used by the compressed counting are orbits of actual isometries (Appendix D.4.6).

*Scope.* Finite check of all 230 forms.

## quaternion-triples.json

The 2547 Gaussian classes with g ~ -g not met by S_gen fall into 849 triples, each one lattice (Proposition 6.5(ii), Lemma 6.6(iii)).

*Scope.* For each of the 2547 Gaussian classes: its key, |U|, f, a Gram matrix G of the underlying lattice of rank 32 and the matrix W of multiplication by i on row vectors (W G W^T = G, W^2 = -1); for each of the 849 lattices: its database label, its three keys and two integer isometries T (T^T G_first T = G_second, det T = +-1) joining them. These are the Gram matrices against which the isometries were found; for 64 keys they are other bases than those of the release asset of Gaussian classes. That the keys and |U| belong to (G, W) and that the Gaussian list is complete and distinct are inputs (Table 10); given these, different triples are different lattices by Proposition 6.5.

## root-coefficients.json

c(A_n) = (-1)^n n!, c(D_n) = (-1)^n (n-1)! (2^(n-1) - 1), c(G_m) = (-1)^m 2^(m-1) (m-1)! agree with all 12934 terms of gaussian-rootless-mass.json.

*Scope.* Comparison of the closed formulas with the term table; the classes with g^2 = -1 in W(F_4) and W(E_8) are checked separately with character tables.

## subtraction-records.json

In the full-group records of S_gen every trivial involution entry has an exact pair identifier, at most one per lattice and only for |Aut L| = 4 mod 8; every Hermitian entry of order 4 has a key, and all keys belong to the list of Proposition 4.6 with equal |U| and f.

*Scope.* Read-only consistency of saved records; full groups and complete class rows are assumptions (Table 10).

## type-8-24-masses.json

Identity (B.1): for P = Leech, A_1^24 and A_2^12, sum over Aut P-orbits of good subspaces W of 1/|Stab(W)| = sum over the 14482 lattices K of a(K, P)/|Aut K| (Appendix B.3, Theorem B.3).

*Scope.* For each K, in list order (list entry 8-24:<i>), |Aut K|, how K was separated from the others and |Aut K| obtained (aut_from: graph, the 720 vectors of norm 2 of K^# generate it and the canonical form of their graph is a complete invariant; orbit, the orbit of N/N' completes that form; qfauto, four lattices whose |Aut K|, computed with qfauto, is the order of no other K, which separates them), and the number of the 270 maximal totally singular subspaces T' for which K + T' has a given number of roots (the root numbers 0, 48 and 72 determine the Niemeier lattice: Leech, A_1^24, A_2^12; the others are shared by several). For the Leech lattice and A_1^24 the stabilizer orders of all orbits of good subspaces; for A_2^12, for each of the 106624 orbits of good 3-dimensional subspaces W_3, the pair [|Stab(W_3)|, g(W_3)] with g(W_3) the number of good W containing W_3, so that the left side is (1/15) sum g/|Stab| (Appendix B.3), and that value. The group orders and the numbers g are inputs (Table 10).

## type-exclusion-lp.json

For each type, nonnegative multipliers y of the inequalities a.x <= b with y.A = 0 and y.b = -1, so the system has no rational solution (Proposition 3.3, reason 4; Appendix B.1.1).

*Scope.* The finite systems and their certificates. That the systems encode the constraints on theta series is the argument of Appendix B.1.1.
