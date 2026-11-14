
#let calcs2 = (
  (id: "A1", files: ("research record: A1/",),
   en: [Selection of the 99 lattices from $S_"gen"$ of [P1, Thm B] (`list` mode): orders outside the range of `IdGroup`; their invariants (order, entries of [P1, §6.1], number of conjugacy classes) and pairwise distinction],
   time: [—]),
  (id: "A2", files: ("research record: A2/",),
   en: [The profile of each lattice: $F^*(G)$, chief series with the action on each factor and splitting, conjugacy classes with centralizers, frame shapes and fixed lattices, orbits on the 146880 minimal vectors],
   time: [about 1 h]),
  (id: "A3", files: ("research record: A3/",),
   en: [Presentations of verified order by Lemmas 2.1–2.4 (`pc_presentation_check`, `radical_presentation_check`, `chain_presentation_check`); all relations and the generation of $G$ checked with integer matrices; orders of quotients and simple sections by coset enumeration on presentations rebuilt from the stored relators (`fpsize1`, `fpsizes`, `fpsimp`); the right-hand sides of the relations and the words for the generators recomputed from the stored matrices for the data file of Table 12, which `verify/verify_p2.py` evaluates in integers],
   time: [—]),
  (id: "A4", files: ("research record: A4/",),
   en: [The sections $Omega_8^+(2)$ and $Omega_10^+(2)$ by the chain of sub-presentations (Lemma 2.4): export of the presentation on a small-degree permutation representation, the chain of indices, and the assembly of the presentation of $G$ (`export`, `count`, `chain`); conjugacy classes of the two largest groups on the permutation representations on minimal vectors (`classes`)],
   time: [about 15 min]),
  (id: "A5", files: ("research record: A5/",),
   en: [For an involution $sigma$ of smallest class (and for chosen elements of order 3): the fixed and anti-fixed lattices, their automorphism groups (bliss), and the images of $C(sigma)$ (`glue`, `glues`); the structure of $"Aut" F$ for the rank-16 component of §5 (`autF`)],
   time: [—]),
  (id: "A6", files: ("research record: A6/",),
   en: [Group names confirmed by explicit isomorphisms (`names` mode; bijective homomorphisms returned by `IsomorphismGroups` and checked)],
   time: [—]),
  (id: "A7", files: ("research record: A7/",),
   en: [The two lattices with an element of type $3 - (9, 14) - 7$ are 3-neighbours: reflections of the discriminant group and the stabilizer orders that identify the double cosets],
   time: [—]),
  (id: "A8", files: ("research record: A8/",),
   en: [The invariant quadratic form on $Z(O_3(G))$ for $|"Aut" L| = 839808$ and the type of each of its 26 non-zero elements (`zform` mode)],
   time: [—]),
  (id: "A9", files: ("research record: A9/",),
   en: [The lattice with layer $"SL"(2,9) compose "SL"(2,9)$: decomposition of the character (`tensor`), the commutant and the two rank-16 sublattices (`split16`), and their identification in the genus $"II"_16 (3^(+8))$],
   time: [—]),
  (id: "A10", files: ("research record: A10/",),
   en: [Fixed-point-free elements of order 5, the matrices $sqrt(5)$ and $phi$ of Lemma 7.1 checked in integers, the action of the generators on $sqrt(5)$, and the icosian part of one lattice (`parts`)],
   time: [—]),
  (id: "A11", files: ("research record: A11/",),
   en: [Minimal norms of the non-zero discriminant classes of the components (`parts`, `allclasses`), and the comparison with the 83 pairs of components of involution type [P1, C4] (second implementation)],
   time: [—]),
  (id: "A12", files: ("research record: A12/",),
   en: [The 79 pages of rank 32 of [NS]: evenness, determinant and roots; for the extremal pages $|"Aut" L|$ (bliss) and isometries to the lattices among the 99 of the same order, with $T^top G T = G'$ and $|det T| = 1$ checked in integers],
   time: [65 min]),
  (id: "A13", files: ("research record: A13/",),
   en: [Tables 2–4 of §3: the numbers of lattices of $S$ by $|"Aut" L|$, subset of [P1, Thm B] and `IdGroup`, counted from the identification lists of [P1, C30] (no expected value is written in the script), with the class of each `IdGroup` (solvable, 2-group, centre, $O_3$, `StructureDescription`) from the Small Groups Library; checked order by order against [P1, Table 9] and against the 99 of §4],
   time: [40 s]),
  (id: "A14", files: ("research record: A14/",),
   en: [The orders of the presentations of the identification run for the nine lattices of $S_"gen"$ with $|"Aut" L| > 10^5$ whose order had not been counted there: generators checked to preserve the Gram matrix and to generate a group of order $|"Aut" L|$ on the minimal vectors, relators checked, then a chain of point stabilizers in a permutation representation of small degree (coset enumeration and Reidemeister–Schreier)],
   time: [about 25 min]),
  (id: "A15", files: ("research record: A15/",),
   en: [The commutant $"End"_H (QQ^32)$ for $H = "Aut" L$, its derived group and the centralizer of the element of order 4 of $F(G)$, for the lattice of [Neb98q, §21]; the element $sqrt(-5)$ of the commutant checked in integers],
   time: [2 min]),
  (id: "A16", files: ("research record: A16/",),
   en: [The maximal subgroups of $"Sp"(6,2)$ and their indices from the GAP character table library (`Maxes`), compared with [ATLAS, $S_6(2)$]],
   time: [—]),
)

#let _calc2-ids = calcs2.map(c => c.id)

#let acalc(..ids) = {
  let xs = ids.pos()
  for x in xs { assert(x in _calc2-ids, message: "unknown computation " + x) }
  [\[]
  xs.map(x => link(label("calc2-" + x), x)).join(", ")
  [\]]
}

#let calc2-table(caption: none, lang: "en") = {
  block(above: 0.8em, below: 0.4em, text(size: 0.92em, caption))
  table(columns: (auto, 1fr, auto), align: (left, left, right),
    table.header([ID], [computation and code], [CPU time]),
    ..calcs2.map(c => (
      [#c.id #label("calc2-" + c.id)],
      {
        c.en
        linebreak()
        text(size: 0.8em, c.files.map(f => raw(f)).join(linebreak()))
      },
      c.time,
    )).flatten()
  )
}
