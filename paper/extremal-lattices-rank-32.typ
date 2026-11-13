#import "style.typ": *
#import "computations.typ": *
#show: paper.with(
  lang: "en",
  title: [The extremal even unimodular lattices of rank 32],
  subtitle: [Their number, and the classification of those with automorphisms other than $plus.minus 1$],
  authors: [roiban1344],
  date: [Version 2026-11-14.1],
)


#show figure.where(kind: table): set block(breakable: true)

#policy-box[
*Version 2026-11-14.1.* A versioned preprint; later versions may correct or extend it. The statements on which the proofs rest and their computational assumptions are listed in Tables 7 and 10.]

#abstract-box[
There are exactly $"11,044,268"$ isometry classes of extremal even unimodular lattices of rank 32. Exactly $"143,329"$ have automorphisms other than $plus.minus 1$; we classify these lattices and determine the order of each automorphism group. The proof first classifies the lattices with additional automorphisms, using complete lists of pairs with an automorphism of odd prime order, an involution other than $plus.minus 1$, or a Gaussian structure. The number of all classes is then obtained from King's rootless mass $m$ by $N=2m+sum_(L in S)(1-2 slash |"Aut" L|)$. Representatives of the lattices with automorphism group ${plus.minus 1}$ are not constructed. The proof assumes the correctness of the finite computations specified in §7.
]

= 1. Introduction <s1>

#policy-box[
*What this section does.* §1.1 recalls what was known about extremal lattices of rank 32. §1.2 states the two main theorems and a corollary. §1.3 explains the method in one page. §1.4 says what the proof assumes. §1.5 describes the organisation.
]

== 1.1 Extremal even unimodular lattices of rank 32

The minimum of an even unimodular lattice of rank $n$ satisfies $min L <= 2 floor(n slash 24) + 2$ [SPLAG, Ch. 7, Cor. 21], and $L$ is _extremal_ if equality holds; extremal lattices do not exist for large $n$ [MOS75, Cor. to Thm 1]. In rank 32 an extremal lattice is one without vectors of norm 2. Rank 32 is the first rank in which the genus of even unimodular lattices is far too large to be listed: a mass-based estimate suggests about $1.2 times 10^9$ classes [ChT26, §1]. It is also the first rank in which extremal lattices are numerous.

Several families were known. Koch and Venkov showed that there are exactly five extremal lattices of neighbour defect 0, the code lattices of the five doubly even self-dual codes of length 32 and minimum weight 8 [KV89, Satz 10] (the classification of the codes is quoted there from Conway and Pless), and continued the study of rank 32 in [KV91]. Koch and Nebe found that there are exactly ten extremal lattices of neighbour defect 8 [KN93, Thm 1]; the twelve codes of their Table 1 give these ten and two of the code lattices of neighbour defect 0 (the non-isometry is stated there with reference to an appendix by Nebe to [KV91]). Lattices with additional structure were classified by mass formulas: the eleven extremal lattices with an automorphism of order 3 without fixed vectors [Hen13, Thm 1, Table 4, Rem. 3(a)], using the mass of Eisenstein lattices [HKN09, Thm 1]; the structures over the Hurwitz order [BN97, Thm 5.1]; and the four extremal lattices from a genus of Hermitian lattices of rank 2 over $ZZ[zeta_17]$ [Jür15, Satz 3.6.9]. The mass of the even unimodular lattices of rank 32 with an automorphism $g$ with $g^2 = -1$ is given in [HKN12, Thm 2]. Chenevier and Taïbi indicate that a classification of the whole genus of rank 32 may now be within reach and announce a forthcoming study [ChT26, §1]. Appendix A locates the 27 previously known lattices with $|"Aut" L| > 2$ in our lists.

Since every lattice has the automorphisms $plus.minus 1$, King's mass $m$ [Kin03] gives the lower bound $N >= 2 m approx 1.0969 times 10^7$ (§2.2). No exact count or classification had been published as of November 2026.

== 1.2 Main results

Let $cal(L)$ be the set of isometry classes of extremal even unimodular lattices of rank 32 (§2.1). Every $L in cal(L)$ has the automorphisms $plus.minus 1$, so $|"Aut" L| >= 2$.

#thm[
*Theorem A.* $|cal(L)| = "11,044,268"$.
]

#thm[
*Theorem B.* Let $S = {L in cal(L) : |"Aut" L| > 2}$. Then $|S| = "143,329"$, and $S = S_"gen" union.sq S_(C_2^2) union.sq S_"CQ"$, where
+ $S_(C_2^2)$ consists of the $"67,284"$ lattices with $"Aut" L = {plus.minus 1, plus.minus sigma} tilde.equiv C_2 times C_2$;
+ $S_"CQ"$ consists of the $"63,499"$ lattices for which $"Aut" L$ is a 2-group with $-1$ as its only involution, that is, cyclic of order $4, 8, 16$ or generalized quaternion of order $8, 16, 32, 64$, with the numbers in Table 8;
+ $S_"gen"$ consists of the remaining $"12,546"$ lattices; each has an automorphism of odd prime order or an involution $sigma != plus.minus 1$ with $|C(sigma)| > 4$ or $sigma$ conjugate to $-sigma$.
The lattices of $S$, with $"Aut" L$, are listed in the database [D], version 2, under the numbering of §6.8.
]

Here $C(sigma)$ is the centralizer of $sigma$ in $"Aut" L$.

#thm[
*Corollary C.* Exactly $"10,900,939"$ extremal even unimodular lattices of rank 32 have automorphism group ${plus.minus 1}$. The distribution of $|"Aut" L|$ over $cal(L)$ is given in Table 9.
]

#ptable(
  columns: 4, align: (left, right, left, left),
  caption: [Table 1. The parts of $S$ in Theorem B.],
  table.header([part], [number], [$"Aut" L$], [proved in]),
  [$S_"gen"$], [$"12,546"$], [computed for each $L$ #calc("C1")], [§6.3],
  [$S_(C_2^2)$], [$"67,284"$], [$C_2 times C_2$], [§6.4, Proposition 6.4\; the number by Lemma 6.6],
  [$S_"CQ"$], [$"63,499"$], [$C_(2^n)$ ($n = 2, 3, 4$) or $Q_(2^n)$ ($n = 3, dots, 6$)], [§6.5, Table 8\; the numbers by Lemma 6.6],
  [$S$], [$"143,329"$], [], [],
  [$cal(L) without S$], [$"10,900,939"$], [${plus.minus 1}$], [Corollary C],
)

#let _x-hand = rgb("#e3f1e6")
#let _x-finite = rgb("#e2edf9")
#let _x-both = gradient.linear((_x-hand, 0%), (_x-hand, 50%), (_x-finite, 50%), (_x-finite, 100%))
#let _x-fill = (hand: _x-hand, finite: _x-finite, cited: luma(236), both: _x-both, result: white, plain: white, note: none)
#let _x-tag = (hand: [W], finite: [F], both: [W+F], cited: [C])
#let _x-anchor(n, a, h) = {
  let (x, y, w, k, b) = n
  let (s, f) = if type(a) == array { a } else { (a, 0.5) }
  if s == "t" { (x + w * f, y) } else if s == "b" { (x + w * f, y + h) } else if s == "l" { (x, y + h * f) } else { (x + w, y + h * f) }
}
#let _x-arrow(p, q) = {
  let (x1, y1) = p
  let (x2, y2) = q
  let l = std.calc.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
  let ux = (x2 - x1) / l
  let uy = (y2 - y1) / l
  place(line(start: (x1 * 1pt, y1 * 1pt), end: (x2 * 1pt, y2 * 1pt), stroke: 0.6pt + luma(80)))
  place(polygon(fill: luma(80), stroke: none, (x2 * 1pt, y2 * 1pt),
    ((x2 - 5 * ux + 2.2 * uy) * 1pt, (y2 - 5 * uy - 2.2 * ux) * 1pt),
    ((x2 - 5 * ux - 2.2 * uy) * 1pt, (y2 - 5 * uy + 2.2 * ux) * 1pt)))
}
#let _x-diagram(width, height, nodes, edges, h: 38) = block(width: width * 1pt, height: height * 1pt, breakable: false, {
  for (a, sa, b, sb) in edges { _x-arrow(_x-anchor(nodes.at(a), sa, h), _x-anchor(nodes.at(b), sb, h)) }
  for (key, n) in nodes {
    let (x, y, w, k, body) = n
    place(dx: x * 1pt, dy: y * 1pt, box(width: w * 1pt, height: h * 1pt, radius: 2pt, inset: 2pt,
      fill: _x-fill.at(k), stroke: if k == "result" { 1pt + luma(40) } else if k == "note" { none } else { 0.5pt + luma(120) },
      align(center + horizon, text(size: 7.5pt, {
        body
        let tag = _x-tag.at(k, default: none)
        if tag != none { std.h(3pt); text(size: 6pt, fill: luma(70), [(#tag)]) }
      }))))
  }
})
#let _x-key(k, body) = box(baseline: 20%, box(width: 16pt, height: 8pt, radius: 1pt, fill: _x-fill.at(k), stroke: 0.5pt + luma(120))) + h(3pt) + text(size: 8pt, [#_x-tag.at(k): #body])

#figure(kind: image, numbering: none, gap: 0.8em,
  _x-diagram(440, 168, h: 46, (
    all: (130, 0, 180, "result", [*all extremal lattices* \ $N = |cal(L)| = "11,044,268"$ \ Theorem A]),
    pm: (0, 64, 200, "plain", [$"Aut" L = {plus.minus 1}$: $"10,900,939"$ \ Corollary C: the difference $N - |S|$; \ these lattices are not constructed]),
    s: (240, 64, 200, "plain", [$S$ ($|"Aut" L| > 2$): $"143,329"$ \ Theorem B: each lattice is constructed \ and listed in [D] with $"Aut" L$]),
    gen: (168, 122, 88, "plain", [$S_"gen"$: $"12,546"$ \ §6.3 \ $"Aut" L$ computed]),
    c22: (260, 122, 88, "plain", [$S_(C_2^2)$: $"67,284"$ \ Proposition 6.4, Lemma 6.6(i) \ $"Aut" L tilde.equiv C_2^2$]),
    cq: (352, 122, 88, "plain", [$S_"CQ"$: $"63,499"$ \ §6.5, Lemma 6.6(ii)(iii) \ cyclic or quaternion]),
  ), (
    ("all", ("b", 0.3), "pm", "t"), ("all", ("b", 0.7), "s", "t"),
    ("s", ("b", 0.2), "gen", "t"), ("s", "b", "c22", "t"), ("s", ("b", 0.8), "cq", "t"),
  )),
  caption: [Figure 1. Where the number $N$ comes from. The lattices of $S$ are classified (Theorem B); the lattices with $"Aut" L = {plus.minus 1}$ are only counted. Identity (5.1) splits $N$ in a different way: each lattice $L$ contributes $2 slash |"Aut" L|$ to $2 m approx "10,968,923.3"$, and each lattice of $S$ contributes the rest, $1 - 2 slash |"Aut" L|$. These remainders add up to about $"9,687.0"$ over $S_"gen"$, to exactly $"33,642"$ over $S_(C_2^2)$ (each lattice gives $1 slash 2$) and to $"512,251" slash 16 approx "32,015.7"$ over $S_"CQ"$ (Table 8), in all about $"75,344.7"$ (§6.7).],
)

== 1.3 The method

Let $G$ be a finite group containing the central involution $-1$. If $G!={plus.minus 1}$, choose $x$ outside that subgroup. If its order has an odd prime divisor, a power of $x$ has odd prime order. Otherwise its order is $2^a$. For $a=1$, $x$ is an involution other than $plus.minus 1$. For $a>=2$, $x^(2^(a-1))$ is an involution; if it is $-1$, a further power has square $-1$. Thus these three kinds of structures detect every lattice of $S$.

We classify the pairs carrying these structures. For odd primes and involutions this means complete lists of components and of gluing classes (§4, Appendix B); Gaussian pairs are complete by their exact rootless mass (§4.3, Appendix C). Completeness requires valid candidates, pairwise non-isometry, correct group orders in the weights, and an independently justified mass. An equality of masses alone does not establish all four conditions.

We then identify the underlying integer lattices, compute their full automorphism groups where needed, and use finite-group arguments for the remaining cyclic, quaternionic and $C_2^2$ cases (§6). This establishes Theorem B without using the total $N$. Finally, summing $1=2 slash |G|+(1-2 slash |G|)$ over all lattices gives the formula $N = 2 m + sum_(L in S) (1 - 2 slash |"Aut" L|)$ of (5.1) and Theorem A (§5 and §6.7).

An earlier organization of the count by assigned group elements is kept in the auxiliary note [Aux]; it is not used in the proof.

== 1.4 What the proof assumes

The proof is computer-assisted. Its numerical inputs are King's rootless mass, complete lists of components and gluing pairs, the Gaussian rootless mass and list, exact identifications of the underlying integer lattices, and their full automorphism-group orders. It relies on complete short-vector enumeration, group and orbit computations, canonical forms and isometry decisions, and our implementations of the stated mass formulas. These dependencies are listed in §7. Returned automorphisms and positive isometries are checked in integers; completeness of groups and negative isometry decisions remain software assumptions.

The two exceptional full-group computations from Gram matrices alone, for $"BW"_32$ and for the code lattice with $|"Aut" L| = 2^31 3^5 5^2 7$, are separated from the structural calculation of their conjugacy classes (§6.3). Group presentations or the order of a subgroup generated by construction matrices do not by themselves provide an upper bound for the full lattice automorphism group. The aggregation (5.1) reads only King's mass and the classified lattices. This arithmetic check is conditional on the mathematical completeness and the group orders of the classification; it is not an independent proof of those inputs.

== 1.5 Organisation

§2 fixes conventions. §3 detects additional automorphisms and gives their possible types. §4 constructs complete lists of structured pairs. §5 reduces the total to the classification, proved in §6. §7 records algorithms and reproducibility. Appendices B, C and D supply the component classifications and mass formulas. The auxiliary note [Aux] keeps an earlier organization of the count. Statements and displayed equations are numbered within each section, so that (4.1) is an equation and Proposition 4.1 a statement; tables are numbered consecutively, with the letter of the appendix in the appendices. The figures are numbered in the same way.

*How to read this paper.* The paper can be read in several ways.
- _For the numbers only._ §1.2 states them, and Figure 1 shows how $N$ splits. §6.8 and Table 9 give the distribution of $|"Aut" L|$, and the database [D] lists the lattices of $S$.
- _To check the proof._ Figure 2 (§5.2) and Table 7 show which statement rests on which. Read §3 and §4 for the complete lists, then §6 for the classification and §5.1 for the count. The appendices prove the component lists (B) and the masses (C, D). Several long proofs begin with a paragraph headed "Strategy of the proof".
- _To use the data._ §6.8 explains the labels and §7.4 the files; the README files of the repository describe them.
- _For the groups themselves._ The second paper of this repository describes the 99 lattices with the largest groups.
The remarks headed "Points that are easy to misread" (§2 and §5) collect conventions on which misreadings are likely.

= 2. Preliminaries <s2>

#policy-box[
*What this section does.* We fix conventions (§2.1), recall the mass and King's rootless mass (§2.2), discriminant forms, 2-elementary lattices and gluing (§2.3), and Hermitian structures given by automorphisms whose minimal polynomial is cyclotomic (§2.4). Everything here is standard; we give references and the precise forms we use.
]

#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 2. Notation.],
  table.header([symbol], [meaning], [§]),
  [$cal(L)$, $N$], [isometry classes of extremal even unimodular lattices of rank 32; $N = |cal(L)|$], [2.1],
  [$m$], [King's mass $sum_(L in cal(L)) 1 slash |"Aut" L|$], [2.2],
  [$S$], [${L in cal(L) : |"Aut" L| > 2}$], [1.2],
  [$C(h)$], [centralizer of $h$ in $"Aut" L$], [2.1],
  [$F, K$, $Z$], [fixed lattice $L^h$, and its orthogonal complement ($K = L^(-sigma)$ for an involution, $Z$ for odd prime order)], [4.1.1],
  [$h$, $H$, $tilde(H)$], [$h$ is an element of $"Aut" L$, as in a pair $(L, h)$; $H$ is the Hermitian form with $(x, y) = "Tr" H(x, y)$ (§2.4), and $tilde(H)$ its rescaling in Lemma 4.4. Appendices C and D, where no element $h$ occurs, write $h$ for their Hermitian form (Table 3)], [2.4],
  [$KK$, $nu$], [the cyclotomic field $QQ(zeta)$ (§2.4) or $QQ(zeta_p)$ (§4.2, Appendix D), not the lattice $K$; $nu = pi.alt overline(pi.alt)$ (Lemma 4.4, Appendix D)], [2.4],
  [$sigma$, $g$], [$sigma != plus.minus 1$ an involution of $"Aut" L$; $g$ an automorphism with $g^(2^(a-1)) = -1$, a Gaussian structure if $g^2 = -1$], [2.4, 3.3],
  [$omega$], [the alternating form $H'(x, y) mod pi.alt$ on $L slash pi.alt L$ (§2.4, and $h(x, y) mod pi.alt$ in Appendix C); in §7.2 a structure of a Hermitian lattice; $zeta_3$ in Appendix D.4], [2.4],
  [$s$, $delta$], [rank of the discriminant group $FF_p^s$; for $p = 2$ the parity invariant of a 2-elementary lattice], [2.3],
  [$q$, $O(q)$], [discriminant quadratic form and its orthogonal group], [2.3],
  [$[phi]$, $"Stab"(phi)$], [double coset of a gluing $phi$ and its stabilizer], [4.1.1],
  [$p - (z, f) - s$], [type of an element $h$ of odd prime order $p$: $"rank"_(ZZ[zeta_p]) (L^h)^perp = z$, $"rank" L^h = f$, $L^h$ has determinant $p^s$], [3.3],
  [$e(h)$], [the number of elements of $"Aut" L$ assigned to $h$ in the auxiliary count; computed inside a full group as a check, not a counting input], [[Aux, §3.2]],
  [$f$], [2 or 1 according as a Hermitian $ZZ[i]$-lattice is not or is isometric to its complex conjugate (§§4.3, 6.1); in a type $p - (z, f) - s$ the rank of $L^h$], [6.1],
  [$"U"$, $U$], [the unitary group $"U"(L, H) = C(g)$ of a Hermitian lattice $(L, H)$, also $"U"(T)$, $"U"(Z)$ for root lattices and components, and $|"U"|$ its order (also the finite unitary groups $"U"_b (FF_4)$ of Appendix D); the italic $U$ is always a unimodular sublattice (Proposition 3.2, Lemma 4.4, Appendices B and D)], [2.4],
  [$D$], [the isometry invariant of §7.2 of a pair $(L, plus.minus g)$; $D(M)$ is a discriminant group; [D] is the database], [7.2],
  [$pi.alt$, $pi$], [$pi.alt$ is a prime element: $1 - zeta$ in §2.4, $1 - zeta_p$ in §4.2 and Appendix D, $1 + i$ in §7.2 and Appendix C (where also $pi.alt_i = 1 - zeta_(2^i)$); $pi$ is the circle constant], [2.4],
  [$M_2$], [rootless mass of Hermitian $ZZ[i]$-lattices of rank 16], [4.3],
  [$(z, s)$, $"mass"^emptyset$], [type of a Hermitian $ZZ[zeta_p]$-lattice and the rootless mass of its genus], [D],
)

== 2.1 Lattices and automorphisms
A _lattice_ is a free $ZZ$-module $L$ of finite rank $n$ with a positive definite symmetric bilinear form $(x, y)$, given by a Gram matrix in a basis; $(x, x)$ is the _norm_ of $x$, and $min L$ the smallest norm of a non-zero vector. The dual is $L^\# = {y in L times.o QQ : (y, L) subset ZZ}$. $L$ is _integral_ if $L subset L^\#$, _even_ if all norms are even, and _unimodular_ if $L = L^\#$. Even unimodular lattices exist exactly in ranks divisible by 8, and their minimum satisfies $min L <= 2 floor(n slash 24) + 2$ [SPLAG, Ch. 7, Cor. 21]; $L$ is _extremal_ if equality holds. In rank 32 this means $min L = 4$, that is, $L$ has no _roots_ (vectors of norm 2).

Two lattices are isometric if a $ZZ$-linear bijection preserves the forms. $cal(L)$ is the set of isometry classes of extremal even unimodular lattices of rank 32, and $N = |cal(L)|$. $"Aut" L$ is the finite group of isometries of $L$ onto itself; it contains $plus.minus 1$, and $-1$ is central. For $h in "Aut" L$, $C(h)$ is its centralizer, $[h]$ its conjugacy class, and $L^h = {x in L : x h = x}$ its fixed lattice (we let $"Aut" L$ act on the right).

== 2.2 Masses
The _mass_ of a set of lattices is $sum 1 slash |"Aut" L|$ over its isometry classes. The mass of all even unimodular lattices of rank $n$ is given by the Smith–Minkowski–Siegel mass formula [SPLAG, Ch. 16]. King computed the mass of those without roots by inclusion–exclusion over the root systems they may contain [Kin03]; for $n = 32$ it is the mass $m$ of $cal(L)$,

#tagged("2.1", $ m = "1,310,037,331,282,023,326,658,917" / "238,863,431,761,920,000" approx "5,484,461.65" $)

== 2.3 Discriminant forms, 2-elementary lattices, gluing
For an even lattice $M$, the _discriminant group_ $D(M) = M^\# slash M$ carries the _discriminant form_ $q_M (x + M) = (x, x) mod 2 ZZ$ and the bilinear form $b_M (x + M, y + M) = (x, y) mod ZZ$, written $b$ when $M$ is clear [Nik79]; $O(q_M)$ is its orthogonal group, and $"Aut" M -> O(q_M)$ has image $overline("Aut" M)$. $M$ is _$p$-elementary_ if $p M^\# subset M$, so that $D(M) tilde.equiv FF_p^s$. For $p = 2$ we put $delta = 0$ if $q_M$ takes values in $ZZ slash 2 ZZ$ and $delta = 1$ otherwise; the even 2-elementary lattices of given rank, $s$ and $delta$ form one genus when they exist [Nik79, Thm 3.6.2]. For odd $p$ the genus of even $p$-elementary lattices is determined by the rank and $s$ ([SPLAG, Ch. 15, Thm 13]; see reason 1 in the proof of Proposition 3.3).

#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 3. Normalizations of forms. $q$ always denotes the discriminant form of an even $ZZ$-lattice with values in $QQ slash 2 ZZ$.],
  table.header([object], [form], [where]),
  [even $ZZ$-lattice $M$], [$q_M (x + M) = (x, x) mod 2 ZZ$ on $D(M)$], [§2.3],
  [even $ZZ$-lattice $M$], [$b_M (x + M, y + M) = (x, y) mod ZZ$ on $D(M)$], [§2.3],
  [$p$-elementary $D(M)$, $p$ odd], [$Q(a) = (p slash 2) q_M (a) in FF_p$, i.e. $Q(a) = j$ when $q_M (a) = 2 j slash p$], [§4.2, §B.1.1, $beta_3$ in §D.4],
  [the same, convention $3 q$], [$3 q_M (a) mod 3 = 2 Q(a)$; determinant characters differ by $(2 slash 3)^(dim) = (-1)^(dim)$], [§D.4.6],
  [2-elementary $D(M)$], [$(rho, delta, sigma)$: rank, $delta$ as above, Brown invariant], [§2.3, §D.4.6],
  [$ZZ[zeta_k]$-structure (§2.4)], [$(x, y) = "Tr"_(KK slash QQ) H(x, y)$], [§2.4],
  [Gaussian lattices (Appendix C)], [the form written $h$ in Appendix C is $2 H$ for the $H$ of §2.4, so $(x, y) = "Re" h(x, y)$; roots $h(x, x) = 2$; $q(x) = h(x, x) mod 2$; $omega(x, y) = h(x, y) mod pi.alt$ on $L slash pi.alt L$, $pi.alt = 1 + i$], [§2.4 ($H' = d H$, $d = 2$), §4.3, App. C],
  [$ZZ[zeta_p]$, $p$ odd (§4.2, App. D)], [$tilde(H) = nu^((p - 3) slash 2) H$ with $(x, y) = "Tr"_(KK slash QQ) H(x, y)$ and $nu = pi.alt overline(pi.alt)$; $a(x, y) = pi.alt tilde(H)(x, y) mod pi.alt$ on $Z slash Z^*$; Appendix D writes $h$ for $tilde(H)$], [Lemma 4.4, §D.1],
  [Eisenstein lattices, $p = 3$ (§D.4)], [$h = tilde(H) = H$, roots $h(x, x) = 1$; $beta_3 (x, y) = 3 h(x, y) mod pi.alt$ on $D_3$; $beta_2 (x, y) = 2 h(x, y) mod 2 cal(O)_2$ on $D_2$], [§D.4.1],
)


If $F$ and $Z$ are primitive sublattices of an even unimodular lattice $L$ with $Z = F^perp$, the projections identify $D(F)$ and $D(Z)$ by an isomorphism $phi$ with $q_Z compose phi = -q_F$, and $L = {x + y : x in F^\#, y in Z^\#, phi(x + F) = y + Z}$; conversely every such $phi$ defines an even unimodular lattice [Nik79, §1.5, Prop. 1.6.1]. Even overlattices of an even lattice $M$ correspond to the subgroups of $D(M)$ on which $q_M$ vanishes [Nik79, Prop. 1.4.1]. §4.1.1 turns this into a correspondence for pairs $(L, h)$ (Proposition 4.1).

#thm[
*Lemma 2.1 (orders of $O(q)$).* Let $F$ be an even positive definite lattice of rank $d$ with discriminant form $q$ on $V = D(F)$.
#set enum(numbering: "(i)")
+ If $F$ is $p$-elementary for an odd prime $p$, with $V tilde.equiv FF_p^s$, then $O(q)$ is the orthogonal group of the non-degenerate quadratic form $Q$ of Table 3, so $|O(q)| = |"GO"_s (p)|$ for odd $s$ and $|O(q)| = |"GO"^epsilon_s (p)|$ for even $s$, where $epsilon = chi((-1)^(s slash 2) "disc" Q)$ and $chi("disc" Q) = e^(2 pi i d slash 8) epsilon_p^(-s)$; here $chi$ is the Legendre symbol modulo $p$, and $epsilon_p = 1$ or $i$ according as $p equiv 1$ or $3$ modulo 4.
+ If $F$ is 2-elementary with invariants $s$ and $delta$, then $|O(q)|$ is $|"GO"^epsilon_s (2)|$ if $delta = 0$, with $epsilon = +$ if and only if $d equiv 0$ (mod 8); $|"GO"^epsilon_(s-1) (2)|$ if $delta = 1$ and $s$ is odd, with $epsilon = +$ if and only if $d equiv plus.minus 1$ (mod 8); $2 |"Sp"_(s-2) (2)|$ if $delta = 1$, $s$ is even and $d equiv 2$ (mod 4); and $2^(s-2) |"GO"^epsilon_(s-2) (2)|$ if $delta = 1$, $s$ is even and $d equiv 0$ (mod 4), with $epsilon = +$ if and only if $d equiv 0$ (mod 8).
Here $"GO"$ is the full orthogonal group of a non-degenerate quadratic form and $"Sp"$ the symplectic group: $|"GO"^epsilon_(2m) (r)| = 2 r^(m(m-1)) (r^m - epsilon) product_(i=1)^(m-1) (r^(2i) - 1)$, $|"GO"_(2m+1) (r)| = 2 r^(m^2) product_(i=1)^m (r^(2i) - 1)$ for odd $r$, and $|"Sp"_(2m) (r)| = r^(m^2) product_(i=1)^m (r^(2i) - 1)$; for $m = 0$, $|"GO"^epsilon_0 (r)| = |"Sp"_0 (r)| = 1$ (the zero space), and the formula for $"GO"_(2m+1) (r)$ gives $|"GO"_1 (r)| = |{plus.minus 1}| = 2$.
]
#proofbox[
Write $gamma(q) = sum_(x in V) e^(pi i q(x))$; by Milgram's formula, $gamma(q) = |V|^(1 slash 2) e^(2 pi i d slash 8)$. For a non-degenerate quadratic form $P$ of dimension $2m$ over $FF_2$ of type $epsilon$, $sum_x (-1)^(P(x)) = epsilon 2^m$; and for a vector $w$, $P(x) + beta(w, x) = P(x + w) + P(w)$, where $beta$ is the polar form, so $P + beta(w, dot)$ has the type of $P$ times $(-1)^(P(w))$.

_(i)._ As $p$ is odd, $q$ and $Q$ determine each other, so $O(q) = O(Q)$. Diagonalizing $Q = sum_j c_j x_j^2$ and using quadratic Gauss sums, $gamma(q) = chi("disc" Q) epsilon_p^s p^(s slash 2)$; compare with Milgram's formula. A form of even dimension $s$ has type $+$ if and only if $(-1)^(s slash 2) "disc" Q$ is a square.

_(ii)._ Here $q$ takes values in $(1 slash 2) ZZ slash 2 ZZ$, and $beta = 2 b$ is a non-degenerate symmetric bilinear form over $FF_2$. Since $x |-> beta(x, x)$ is linear, there is a unique $c$ with $beta(x, x) = beta(x, c)$ for all $x$; it is fixed by $O(q)$, and $delta = 0$ if and only if $c = 0$. Put $W = c^perp = {x : q(x) in ZZ slash 2 ZZ}$; on $W$, $Q = q|_W$ is an $FF_2$-valued quadratic form with polar form $beta$, and the radical of $beta|_W$ is $W inter chevron.l c chevron.r$.

If $delta = 0$, then $W = V$, $O(q) = O(Q)$ and $gamma(q) = epsilon 2^(s slash 2)$.

If $delta = 1$ and $s$ is odd, then $c in.not W$, since otherwise $W slash chevron.l c chevron.r$ would carry a non-degenerate alternating form in the odd dimension $s - 2$. So $V = W perp chevron.l c chevron.r$ with $q(c) = plus.minus 1 slash 2$, $Q$ is non-degenerate on $W$, and $O(q) = O(Q)$: isometries fix $c$, and an isometry of $Q$ extended by $c |-> c$ preserves $q$. Moreover $gamma(q) = epsilon 2^((s-1) slash 2) (1 + e^(plus.minus pi i slash 2))$, which gives $d equiv plus.minus 1$ (mod 8) for $epsilon = +$ and $d equiv 4 plus.minus 1$ (mod 8) for $epsilon = -$.

If $delta = 1$ and $s$ is even, then $c in W$, since $beta|_W$ is alternating on a space of odd dimension; so $q(c) in {0, 1}$, and every $x in.not W$ has $q(x) in {1 slash 2, 3 slash 2}$. Restriction gives $rho : O(q) -> O(Q)$. For $x in.not W$ put $P_x (w) = q(x + w) - q(x) = Q(w) + beta(x, w)$ on $W$. Fix $e in.not W$ and $g in O(Q)$. The vectors $v$ with $beta(v, g w) = beta(e, w)$ for all $w in W$ are $v_0$ and $v_0 + c$, both outside $W$; the map equal to $g$ on $W$ with $e |-> v$ lies in $O(q)$ if and only if $q(v) = q(e)$, and $q(v_0 + c) = q(v_0) + q(c) + 1$.

_Case $q(c) = 0$._ Exactly one of $v_0$, $v_0 + c$ has the value $q(e)$, so $rho$ is bijective. $Q$ induces a non-degenerate form $overline(Q)$ on $W slash chevron.l c chevron.r$; $O(Q) -> O(overline(Q))$ is onto, with kernel ${w |-> w + lambda(w) c}$, $lambda in (W slash chevron.l c chevron.r)^*$, of order $2^(s-2)$. Hence $|O(q)| = 2^(s-2) |"GO"^epsilon_(s-2) (2)|$ with $epsilon$ the type of $overline(Q)$. As $P_e (c) = 1$, the sum over $e + W$ in $gamma(q)$ vanishes, so $gamma(q) = 2 sum_(W slash chevron.l c chevron.r) (-1)^(overline(Q)) = epsilon 2^(s slash 2)$: $d equiv 0$ (mod 4), and $epsilon = +$ if and only if $d equiv 0$ (mod 8).

_Case $q(c) = 1$._ $v_0$ and $v_0 + c$ have the same value, so $ker rho = {1, t}$ with $t(x) = x + beta(x, c) c$. For $x in.not W$, $P_x (c) = 0$, so $P_x$ is a non-degenerate form on $W slash chevron.l c chevron.r$, and $P_(x + w) = P_x + beta(w, dot)$. By the remark above, the type of $P_(e + w)$ is that of $P_e$ times $(-1)^(q(e + w) - q(e))$; hence $q(x)$ is determined by the type of $P_x$, for all $x in.not W$. Now $P_(v_0) compose g = P_e$, so $P_(v_0)$ and $P_e$ have the same type, $q(v_0) = q(e)$, and $rho$ is onto. Since $Q(c) = 1$, restriction to $W slash chevron.l c chevron.r$ is an isomorphism $O(Q) -> "Sp"(W slash chevron.l c chevron.r) tilde.equiv "Sp"_(s-2) (2)$: a symplectic $h$ lifts uniquely, by sending $w$ to the lift of $h(overline(w))$ with value $Q(w)$. Hence $|O(q)| = 2 |"Sp"_(s-2) (2)|$. Here the sum over $W$ in $gamma(q)$ vanishes (pair $w$ with $w + c$), and the sum over $e + W$ is $plus.minus 2^(s slash 2) e^(pi i q(e))$, so $d equiv 2$ (mod 4). #qed
]

#rem[
*Example (Lemma 2.1(ii) in small cases).* This example is not used in the proofs. For $F = A_1^2$ ($d = 2$), $D(F) tilde.equiv FF_2^2$, and the three non-zero classes have $q = 1 slash 2, 1 slash 2, 1$. So $delta = 1$, $s = 2$ and $d equiv 2$ (mod 4), and the lemma gives $|O(q)| = 2 |"Sp"_0 (2)| = 2$: the only non-trivial isometry exchanges the two classes with $q = 1 slash 2$. For $F = D_4$ ($d = 4$), the three non-zero classes all have $q = 1$. So $delta = 0$, $s = 2$, and since $d equiv.not 0$ (mod 8) the type is $epsilon = -$; the lemma gives $|"GO"^-_2 (2)| = 6$, all permutations of the three classes. For $F = sqrt(2) E_8$, the component of type $(8, 24)$ (§4.1.3), $q(x) = (x, x) slash 2 mod 2$ on $E_8 slash 2 E_8$, so $s = 8$, $delta = 0$, $d = 8$, and $|O(q)| = |"GO"^+_8 (2)| = "348,364,800" = |W(E_8)| slash 2$, in accordance with §4.1.3. (The first two values can be checked by listing the six bijections of the three non-zero classes.)
]

== 2.4 Hermitian structures
If $g in "Aut" L$ has minimal polynomial the cyclotomic polynomial $Phi_k$, then $L$ is a torsion-free module over $ZZ[zeta_k]$ by $zeta_k x = x g$, of rank $32 slash phi(k)$, and $(x, y) = "Tr" H(x, y)$ for a unique Hermitian form $H$. An isometry $L -> L'$ conjugating $g$ to $g'$ is the same as an isometry of Hermitian lattices, so the pairs $(L, g)$ up to conjugacy correspond to the isometry classes of Hermitian lattices, and $C(g)$ is the unitary group $"U"(L, H)$. The pair $(L, g^(-1))$ corresponds to the complex conjugate lattice $(L, overline(H))$; it is isometric to $(L, H)$ exactly when $g$ is conjugate to $g^(-1)$ in $"Aut" L$. For $k = 2^a$, $g^(-1) = -g^(2^(a-1) - 1)$, which is $-g$ for $a = 2$. For $k=4$, the rootless mass $M_2$ is $sum 1 slash |"U"(L,H)|$ over the Gaussian classes whose trace form is extremal (§4.3).

Let $a >= 2$, $KK = QQ(zeta)$ with $zeta$ a primitive $2^a$-th root of unity, $cal(O) = ZZ[zeta]$ its ring of integers, $KK^+ = QQ(zeta + zeta^(-1))$ and $cal(O)^+ = cal(O) inter KK^+$. The element $pi.alt = 1 - zeta$ generates the unique prime of $cal(O)$ above 2, and $cal(O) slash pi.alt cal(O) = FF_2$. Let $g in "Aut" L$ satisfy $g^(2^(a-1)) = -1$. Then $L$ is an $cal(O)$-module via $zeta x = x g$, torsion free of rank $r = 32 slash 2^(a-1)$, and since $g$ is an isometry there is a unique form $H : L times L -> KK$, $cal(O)$-linear in the first variable and Hermitian, with $(x, y) = "Tr"_(KK slash QQ) H(x, y)$. We call $r$ the rank of $L$ over $cal(O)$.

#thm[
*Lemma 2.2.* If $L$ is even unimodular, then $r$ is even.
]
#proofbox[
Let $frak(D)$, $frak(D)^+$ and $frak(D)_(KK slash KK^+)$ be the differents of $KK slash QQ$, $KK^+ slash QQ$ and $KK slash KK^+$, so that $frak(D) = frak(D)_(KK slash KK^+) frak(D)^+$. Since $cal(O) = cal(O)^+[zeta]$ and $zeta$ is a root of $X^2 - (zeta + zeta^(-1)) X + 1$, we have $frak(D)_(KK slash KK^+) = (zeta - zeta^(-1)) = (1 - zeta^2) = (1 - zeta)(1 + zeta) cal(O) = pi.alt^2 cal(O)$, because $-zeta$ is also a primitive $2^a$-th root of unity. Put $frak(p)^+ = pi.alt^2 cal(O) inter cal(O)^+$, the prime of $cal(O)^+$ above 2.

Since $H(x, L)$ is an $cal(O)$-module, $(x, L) subset ZZ$ if and only if $H(x, L) subset frak(D)^(-1)$. Choose a generator $d$ of $frak(D)$ and put $H' = d H$. Unimodularity gives $H'(L, L) subset cal(O)$, and if $H'(x, L) subset pi.alt cal(O)$ then $H(pi.alt^(-1) x, L) subset frak(D)^(-1)$, so $pi.alt^(-1) x in L$. Hence $omega(x, y) = H'(x, y) mod pi.alt$ is a non-degenerate form on $L slash pi.alt L tilde.equiv FF_2^r$. It is bilinear and symmetric, since complex conjugation and every unit act trivially on $cal(O) slash pi.alt = FF_2$.

It remains to show that $omega$ is alternating. Fix $x in L$ and put $c = H(x, x) in KK^+$. For $alpha in cal(O)$, evenness of $L$ applied to $alpha x$ and to $(1 + alpha) x$ gives $2 "Tr"_(KK^+ slash QQ)(alpha overline(alpha) c) in 2 ZZ$ and $2 "Tr"_(KK^+ slash QQ)((alpha + overline(alpha)) c) in 2 ZZ$. The elements $alpha overline(alpha)$ and $alpha + overline(alpha)$ ($alpha in cal(O)$) span $cal(O)^+$ over $ZZ$: the traces $alpha + overline(alpha)$ form an ideal containing $2$ and $zeta + zeta^(-1)$, which is $frak(p)^+$ (for $a = 2$ it is $2 ZZ$, and for $a >= 3$ the element $zeta + zeta^(-1) equiv -pi.alt overline(pi.alt) med (mod 2)$ has valuation 1 at $frak(p)^+$), and $1 = 1 dot overline(1)$ generates $cal(O)^+ slash frak(p)^+ = FF_2$. Hence $"Tr"_(KK^+ slash QQ)(t c) in ZZ$ for all $t in cal(O)^+$, that is $c in (frak(D)^+)^(-1)$, and $H'(x, x) = d c in frak(D) (frak(D)^+)^(-1) = pi.alt^2 cal(O)$. So $omega(x, x) = 0$. A non-degenerate alternating form exists only in even dimension. #qed
]

For $a = 6$ the rank is $r = 1$, so by Lemma 2.2 no even unimodular lattice of rank 32 has an automorphism $g$ with $g^32 = -1$. For $a >= 7$ the degree $2^(a-1)$ of the minimal polynomial exceeds 32.

#rem[
*Three normalizations of a Hermitian form.* Table 3 lists them; here is the reason for each, with the smallest examples. This remark is not used in the proofs.
- §2.4 uses the form $H$ with $(x, y) = "Tr" H(x, y)$, the trace to $QQ$; it compares directly with the integer lattice.
- For $ZZ[i]$ (§4.3, Appendix C) the form $h = 2 H$ is used, so that $(x, y) = "Re" h(x, y)$ and a root has $h(x, x) = 2$, as in the integer lattice. For example, in $A_1 [i] = A_1 times.o ZZ[i]$ (Proposition C.1, an $A_1^2$ over $ZZ$) a root $x$ has $(x, x) = 2$, $H(x, x) = 1$ and $h(x, x) = 2$.
- For an element of odd prime order $p$ (Lemma 4.4, Appendix D) the rescaled form $tilde(H) = nu^((p - 3) slash 2) H$ is used, because the lattices $U$ of Lemma 4.4 are unimodular for it, and the theory of unimodular Hermitian lattices then applies. For $p = 3$, $tilde(H) = H$; for example $A_2 = ZZ[omega]$ with $(x, y) = "Tr"(x overline(y))$ has $H(1, 1) = 1$.
]

#rem[
*Points that are easy to misread (conventions).* These remarks add nothing to the proofs; they collect conventions of this paper on which misreadings are likely.
+ _Extremal means minimum 4._ In rank 32 an extremal lattice is exactly an even unimodular lattice without roots (§2.1). "Rootless" and "extremal" mean the same here.
+ _$|"Aut" L|$ counts $plus.minus 1$._ Every group order in this paper is the order of the full group $"Aut" L$, which contains $-1$. So $|"Aut" L| >= 2$, and the "smallest possible group" is ${plus.minus 1}$, of order 2.
+ _A mass is not a number of lattices._ The mass $m$ counts each lattice $L$ with weight $1 slash |"Aut" L|$. As $|"Aut" L| >= 2$, we have $2 m <= N$, with equality only if every group is ${plus.minus 1}$. Here $2 m approx 1.0969 times 10^7$, while $N approx 1.1044 times 10^7$.
+ _Pairs are not lattices._ The lists of §4 list pairs $(L, h)$ up to isometry, with weight $1 slash |C(h)|$. One lattice occurs once for each conjugacy class of its structures (§3.1), so a list of pairs may contain the same lattice several times. §6 groups the pairs lattice by lattice.
+ _$sigma$ and $-sigma$; $g$ and $-g$._ An involution $sigma$ and $-sigma$ give the same lattice with $F$ and $K$ exchanged (§4.1.2). A Gaussian structure $g$ and $-g$ give complex conjugate Hermitian lattices (§2.4). This is why the weight $f slash |"U"|$ of §4.3 has $f = 2$ or $1$.
+ _Letters._ An element of $"Aut" L$ is $h$ and a Hermitian form is $H$; the lattice $L^(-sigma)$ is $K$ and a cyclotomic field is $KK$; the unitary group is $"U"$ (upright) and a unimodular sublattice is $U$ (italic): in this paper the upright $"U"$ always denotes a unitary group and the italic $U$ always a lattice; a prime element is $pi.alt$ and $pi$ is the circle constant. Two letters keep two meanings in the main text, because the certificates and the data use them: $f$ is the rank of $L^h$ in a type $p - (z, f) - s$ and also the factor 1 or 2 of a Gaussian class, and $D$ is a discriminant group and also the invariant of §7.2 ([D] is the database). Table 2 lists the conventions of the main text, including the three uses of $omega$. Some letters have a further meaning inside one proof or one appendix, stated where it is used: for instance $h$ is the Hermitian form of Appendices C and D, $N$ a rank in Appendix D, and $S$, $T$ the Gram matrices of Siegel's formula in Appendix C.
]

= 3. Detecting additional automorphisms <s3>

== 3.1 The three kinds of structure

A _pair_ $(L,h)$ consists of $L in cal(L)$ and $h in "Aut" L$. An isometry of pairs $(L, h) -> (L', h')$ is an isometry $L -> L'$ that conjugates $h$ to $h'$. For fixed $L$, these classes are exactly the conjugacy classes in $"Aut" L$.

#thm[
*Proposition 3.1.* If $"Aut" L!={plus.minus 1}$, then $L$ has an automorphism of odd prime order, an involution other than $plus.minus 1$, or an automorphism $g$ with $g^2=-1$.
]
#proofbox[
Choose $x in "Aut" L without {plus.minus 1}$. If its order has an odd prime divisor $p$, a power of $x$ has order $p$. Otherwise its order is $2^a$. For $a=1$, $x$ is the required involution. For $a>=2$, $x^(2^(a-1))$ is an involution. If it is not $-1$, it is the required one; if it is $-1$, put $g=x^(2^(a-2))$. #qed
]

== 3.2 From structured pairs to the classification

Take complete lists of the pairs in Proposition 3.1, identify their underlying integer lattices, and retain one representative per isometry class. Every lattice of $S$ occurs, and every valid representative belongs to $S$. Thus proving that these representatives are distinct and determining their full group orders proves Theorem B. §6 gives the group-theoretic reductions that make this identification feasible. §5 then gives Theorem A from this classification and King's mass.

An earlier organization of the count by assigned group elements is given in the auxiliary note [Aux]; its numbers and higher cyclotomic mass formulas are not premises of this argument.

== 3.3 The types
#thm[
*Proposition 3.2.* Let $sigma != plus.minus 1$ be an involution of $L in cal(L)$, $F = L^sigma = {x in L : x sigma = x}$, $K = L^(-sigma)$ and $d = "rank" F$. Then $F$ and $K$ are even 2-elementary lattices of minimum at least 4 with isomorphic discriminant groups $(ZZ slash 2)^s$, and
$ d in {8, 12, 14, 15, 16, 17, 18, 20, 24}. $
]
#proofbox[
$F$ and $K$ are primitive sublattices of $L$ and orthogonal complements of each other, so the projections of the unimodular lattice $L$ to $F times.o QQ$ and $K times.o QQ$ are $F^\#$ and $K^\#$, and $L slash (F perp K) tilde.equiv F^\# slash F tilde.equiv K^\# slash K$ [Nik79, §1.5, Prop. 1.6.1]. Since $2 x = (x + x sigma) + (x - x sigma) in F perp K$ for $x in L$, these groups are 2-elementary. As sublattices of $L$, $F$ and $K$ are even of minimum at least 4.

The form $psi(x, y) = (x, y) mod 2$ on $F slash 2 F^\#$ is well defined and alternating since $F$ is even, and non-degenerate since ${x in F : (x, F) subset 2 ZZ} = F inter 2 F^\# = 2 F^\#$. Let $T$ with $2 F^\# subset T subset F$ be the preimage of a maximal isotropic subspace; then $[T : 2 F^\#] = 2^((d - s) slash 2)$ and $det T = 4^d 2^(-s) 2^(-(d - s)) = 2^d$. Hence $U = 2^(-1 slash 2) T$ is an integral unimodular lattice of rank $d$ with $min U >= min F slash 2 >= 2$, that is, without vectors of norm 1 (this is [Neb13b, Lemma 4.9]). The same holds for $K$ with rank $32 - d$. Unimodular lattices without vectors of norm 1 exist in rank $n$ only for $n = 0, 8, 12$ and $n >= 14$ [SPLAG, Table 16.7]. Since $sigma != plus.minus 1$, $0 < d < 32$, and both $d$ and $32 - d$ must lie in that set. #qed
]

For an element $h$ of odd prime order $p$, let $F = L^h$ of rank $f$ and $Z = F^perp inter L$. Then $Z$ is a module over $ZZ[zeta_p]$ via $h$, of rank $z$ with $f + (p - 1) z = 32$, and $F^\# slash F tilde.equiv Z^\# slash Z tilde.equiv (ZZ slash p)^s$ [Neb13b, Lemma 4.2]. We call $p - (z, f) - s$ the _type_ of $h$ [Neb13b, Def. 4.1].

#thm[
*Proposition 3.3.* The type of an element of odd prime order of a lattice in $cal(L)$ is one of the 19 types in Table 4.
]
#proofbox[
Let $h$ have type $p - (z, f) - s$. Then $p <= 31$ and $1 <= z <= 32 slash (p - 1)$, $f = 32 - (p - 1) z$, $0 <= s <= min(z, f)$ ($s <= f$ since $F^\# slash F tilde.equiv FF_p^s$ is a quotient of $F^\# slash p F^\# tilde.equiv FF_p^f$, and $s <= z$ by Proposition D.1), and $z - s$ is even by Lemma 4.4. There are 93 such quadruples with $p in {3, 5, 7, dots, 31}$. Appendix B.1 (Table B1) excludes 74 of them, each for one of the following reasons; the remaining 19 are the types of Table 4.
+ _No genus_ (7 quadruples): there is no genus of even $p$-elementary lattices of rank $f$ (or $(p - 1) z$) and determinant $p^s$; for odd $p$ such a genus exists if and only if the rank $n$ is even, $n + s (p - 1) equiv 0 med (mod 4)$, and $8 | n$ when $s in {0, n}$. This is [SPLAG, Ch. 15, Thm 13(a) and Table 15.4] for the genus $"II"_(n, 0)(p^(plus.minus s))$ in the notation there, where the signature is $(n, 0)$ and $p^(plus.minus s)$ is the $p$-part of the discriminant group.
+ _Density_ (59): $F$ or $Z$ has minimum at least 4 and determinant $p^s$, so its centre density is at least $p^(-s slash 2)$, which exceeds the Cohn–Elkies bound on the centre density in its dimension [CE03, Table 3]; the bounds of that table are proved by exact rational arithmetic [CE03, §7]. In each of the 59 cases $p^(-s slash 2)$ exceeds the printed bound by at least $3.7 times 10^(-3)$ (a factor of at least $1.034$), far more than the rounding of the table to five decimals #calc("C2").
+ _No rootless unimodular sublattice_ (2): $Z$ contains a Hermitian unimodular lattice $U$ of the same rank (Lemma 4.4), which has no roots since $L$ is extremal; but the mass of rootless Hermitian unimodular lattices of that rank is 0 #calc("C2").
+ _Linear programming_ (6): the theta series of $F$ and $Z$ are modular forms for $Gamma_0(p)$ with character, exchanged with those of $sqrt(p) F^\#$ and $sqrt(p) Z^\#$ by the Fricke involution; with the constraints from $min F, min Z >= 4$, Lemma 4.2 and Lemma 4.3, the exact linear program (over $QQ$) for their coefficients is infeasible; §B.1.1 gives a rational certificate of infeasibility for each of the six #calc("C2", "C3").
]

#rem[
*Remark.* The proof uses only the linear programs for the six quadruples of reason 4. They are also excluded by complete lists of components (proved complete by their masses): for five of them the genus of $F$ contains no lattice of minimum at least 4, and for $3 - (8, 16) - 6$, whose certificate in §B.1.1 uses only the conditions on $F$ (so that this genus of $F$ has no lattice of minimum at least 4 either), no unimodular sublattice $U$ of rank 8 admits a suitable overlattice $Z$ #calc("C2"). This is a cross-check and is not used.
]


#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 4. The possible detecting structures. For an involution, $d = "rank" L^sigma$; $sigma$ and $-sigma$ have types $(d, 32 - d)$ and $(32 - d, d)$. For an element of odd prime order, the type is $p - (z, f) - s$ (Table 2).],
  table.header([kind], [types], [§]),
  [involution $!= plus.minus 1$], [$(8, 24)$, $(12, 20)$, $(14, 18)$, $(15, 17)$, $(16, 16)$], [4.1],
  [order 3], [$3 - (6, 20) - 6$, $3 - (8, 16) - 8$, $3 - (9, 14) - 7$, $3 - (9, 14) - 9$, $3 - (10, 12) - 6$, $3 - (10, 12) - 8$, $3 - (10, 12) - 10$, $3 - (12, 8) - 8$, $3 - (13, 6) - 5$, $3 - (16, 0) - 0$], [4.2],
  [order 5], [$5 - (4, 16) - 4$, $5 - (5, 12) - 5$, $5 - (6, 8) - 4$, $5 - (6, 8) - 6$, $5 - (7, 4) - 3$, $5 - (8, 0) - 0$], [4.2],
  [order 7, 17, 31], [$7 - (4, 8) - 4$, $17 - (2, 0) - 0$, $31 - (1, 2) - 1$], [4.2],
  [$ZZ[i]$], [rank 16; $g^2=-1$], [4.3],
)

= 4. Complete lists of structured pairs <s4>

#policy-box[
*What this section does.* We list components and gluings for involutions and odd prime order (§§4.1–4.2), and Gaussian pairs (§4.3). Completeness follows from exhaustive reductions or exact positive mass identities, provided each candidate is valid, listed classes are distinct and their weights are correct. §6 then identifies the underlying integer lattices.
]
== 4.1 Pairs of prime order and involutions
=== 4.1.1 Gluings of prime order

Let $h$ be an automorphism of prime order $p$ of $L in cal(L)$, $h != -1$ if $p = 2$. Put $F = L^h$ and $Z = F^perp inter L$; for $p = 2$ we write $K = Z = L^(-h)$. Let $G_F = "Aut" F$, let $G_Z$ be the centralizer of $h|_Z$ in $"Aut" Z$ (for $p = 2$, $h|_Z = -1$ and $G_Z = "Aut" Z$), and let $overline(G)_F$, $overline(G)_Z$ be their images in the orthogonal groups of the discriminant forms $q_F$, $q_Z$. The lattice $L$ is the gluing
$ L = L_phi = {x + y : x in F^\#, y in Z^\#, phi(x + F) = y + Z} $
along the isomorphism $phi : F^\# slash F -> Z^\# slash Z$ induced by the projections, which satisfies $q_Z compose phi = -q_F$.

#thm[
*Proposition 4.1.* Fix $F$ and $(Z, h|_Z)$ up to isometry. Let $"Iso"(q_F, -q_Z)$ be the set of isomorphisms $phi$ with $q_Z compose phi = -q_F$.
+ For every $phi in "Iso"(q_F, -q_Z)$, $L_phi$ is an even unimodular lattice and $"id"_F perp h|_Z$ extends to an automorphism $h_phi$ of $L_phi$ of order $p$ with fixed lattice $F$.
+ $phi |-> (L_phi, h_phi)$ induces a bijection between the double cosets $overline(G)_Z backslash "Iso"(q_F, -q_Z) slash overline(G)_F$ and the isometry classes of pairs $(L, h)$ with these $F$ and $(Z, h|_Z)$, where $L$ runs over all even unimodular lattices of rank 32.
+ The stabilizer $"Stab"(phi) = {(u, v) in G_F times G_Z : overline(v) phi overline(u)^(-1) = phi}$ is isomorphic to $C(h_phi)$.
+ Consequently
#tagged("4.1", $ sum_([phi]) 1 / (|"Stab"(phi)|) = (|"Iso"(q_F, -q_Z)|) / (|G_F| |G_Z|). $)
]
#proofbox[
(i) $L_phi$ is even unimodular since $phi$ is an anti-isometry of discriminant forms [Nik79, Prop. 1.4.1, 1.6.1]. The element $h|_Z$ acts trivially on $Z^\# slash Z$: for $p = 2$ because $Z^\# slash Z$ is 2-elementary, and for odd $p$ because $(1 - h) Z^\# subset Z$ [Neb13b, proof of Lemma 4.2]. Hence $"id"_F perp h|_Z$ preserves $L_phi$. (ii) and (iii): an isometry $L_phi -> L_(phi')$ commuting with $h_phi$, $h_(phi')$ maps $F$ to $F$ and $Z$ to $Z$, so it is a pair $(u, v) in G_F times G_Z$ with $overline(v) phi overline(u)^(-1) = phi'$, and conversely. Taking $phi = phi'$ gives (iii). (iv) is the orbit–stabilizer formula for the action of $G_F times G_Z$ on $"Iso"(q_F, -q_Z)$. #qed
]

When $"Iso"(q_F, -q_Z)$ is non-empty, $|"Iso"(q_F, -q_Z)| = |O(q_F)|$. Identity (4.1) is what makes the lists complete: the right side is computed from $|O(q_F)|$, $|G_F|$ and $|G_Z|$ alone, every term of the left side is positive, so a list of double cosets whose left side equals the right side contains every double coset. Every pair $(L, h)$ with $L in cal(L)$ arises from components $F$ and $Z$ of minimum at least 4, so we list all components of minimum at least 4 in the relevant genera, list all double cosets for each pair of components, prove the lists complete by (4.1), and then keep the double cosets for which $L_phi$ has no vector of norm 2. Each component list is in turn proved complete in Appendix B: for involutions by an enumeration that provably reaches every lattice in question (for type $(8, 24)$ together with identity (B.1), Theorem B.3), for elements of odd prime order by an exact equality with the mass of the genus or of its rootless part.

#rem[
*Example (a gluing in rank 16).* This example is not used in the proofs; it shows Proposition 4.1 and identity (4.1) on the smallest familiar case. Let $L = E_8 perp E_8$ and let $sigma$ exchange the two summands, $(x, y) sigma = (y, x)$. Then $F = L^sigma = {(x, x)}$ and $K = L^(-sigma) = {(x, -x)}$ ($x in E_8$) are both isometric to $sqrt(2) E_8$. Their duals are $F^\# = {(x, x) slash 2}$ and $K^\# = {(x, -x) slash 2}$, so $D(F) tilde.equiv D(K) tilde.equiv E_8 slash 2 E_8$. The gluing $phi$ sends $(x, x) slash 2 + F$ to $(x, -x) slash 2 + K$, and $(x, x) slash 2 + (x, -x) slash 2 = (x, 0)$, so $L_phi$ is $E_8 perp E_8$ again.

No other gluing of these two components exists. The argument of Proposition 4.1 works in any rank, and the right side of (4.1) is $|O(q_F)| slash |"Aut" F|^2$, where $|"Aut" F| = |W(E_8)| = "696,729,600"$ and $|O(q_F)| = |"GO"^+_8 (2)| = |W(E_8)| slash 2$ (the example after Lemma 2.1). So the right side is $1 slash (2 |W(E_8)|)$. On the left, $C(sigma)$ consists of the pairs $(w, w)$ with $w in W(E_8)$ and their products with $sigma$, so $|C(sigma)| = 2 |W(E_8)|$. This single term already equals the right side, and every term is positive. Hence an even unimodular lattice of rank 16 with an involution whose fixed lattice and its complement are both $sqrt(2) E_8$ is $E_8 perp E_8$; for instance $D_16^+$ has no such involution. The same mechanism, that $W(E_8) -> O(q_F)$ is onto, gives exactly one double coset for each lattice $K$ of type $(8, 24)$ (§4.1.3).
]

=== 4.1.2 Orientation

The involutions $sigma$ and $-sigma$ give the same integer lattice. If $F tilde.not K$, the inverse gluing belongs to the list with the components interchanged, so one of the two orientations suffices. If $F tilde.equiv K$, inversion sends $[phi]$ to $[phi^(-1)]$ within the same list. These classes coincide exactly when $sigma$ is conjugate to $-sigma$. §6.1 uses this relation to define involution entries. No element-assignment count is needed for this identification.

=== 4.1.3 The five types

Let $sigma != plus.minus 1$ be an involution of $L in cal(L)$ of type $(d, 32 - d)$, $d <= 16$ (§4.1.2). The components $F = L^sigma$ and $K = L^(-sigma)$ are even, 2-elementary ($2 F^\# subset F$), of minimum at least 4, and their discriminant forms are anti-isometric. Hence they have the same invariants $s$ ($F^\# slash F tilde.equiv FF_2^s$) and $delta$ ($delta = 0$ if $q_F (x) = (x, x) mod 2 ZZ$ takes integral values, $delta = 1$ otherwise). Appendix B.2 lists all such lattices in the ranks that occur (Table B2); for $d <= 15$ there is exactly one $F$, of the form $sqrt(2) U$ with $U$ unimodular. Table 5 shows the components of each type.

*Extremality.* By the proof of Proposition 3.2, $F$ contains $sqrt(2) U$ with $U$ unimodular and $min U >= 2$, so $F^\# subset U slash sqrt(2)$ and every non-zero vector of $F^\#$ has norm at least 1; the same holds for $K$. A vector of $L_phi$ outside $F perp K$ is $x + y$ with $x in F^\# without F$, $y in K^\# without K$, and if it has norm 2 then $(x, x) = (y, y) = 1$. Vectors of $F perp K$ have norm 0 or at least 4. Hence _$L_phi$ is extremal if and only if there is no class $a in F^\# slash F$ such that both $a$ and $phi(a)$ contain a vector of norm 1._ This depends only on the double coset of $phi$ and is decided from the lists of the classes of norm 1 in $F^\# slash F$ and $K^\# slash K$.

*Gluings.* For each pair $(F, K)$ of components, the double cosets of Proposition 4.1(ii) were listed by computing a canonical form of $phi$ under $overline(G)_K times overline(G)_F$ (§7.2) and generating new classes by composing $phi$ with reflections $x |-> x + 2 b(x, v) v$ of $F^\# slash F$ in elements with $q_F (v) equiv 1$, and by random walks, among them walks that keep a given element of a known stabilizer #calc("C4").

How the classes were found plays no role: for each of the 83 pairs of components of types $(15, 17)$, $(14, 18)$, $(12, 20)$ and $(16, 16)$, the sum of $1 slash |"Stab"(phi)|$ over the listed classes equals, as a fraction, the value of the right side of (4.1) computed below, where $|"Iso"(q_F, -q_K)| = |O(q_F)|$ since the list is not empty. The value of $|O(q_F)|$ was computed as $R(q_F)$, where $R(q) = n(q) R(q|_(H^perp))$, $n(q)$ is the number of ordered pairs $(x, y)$ with $q(x) = q(y) = 0$ and $b(x, y) = 1 slash 2$, and $H = chevron.l x, y chevron.r$ for one such pair; for forms on a group of rank at most 3 and for forms without such pairs, $R(q) = |O(q)|$ is counted exhaustively. The stabilizer of $(x, y)$ in $O(q)$ is isomorphic, by restriction to $H^perp$, to $O(q|_(H^perp))$, and its orbit has at most $n(q)$ elements, so $|O(q)| <= n(q) |O(q|_(H^perp))|$ and, by induction on the rank of the group, $R(q) >= |O(q)|$, whether or not $O(q)$ is transitive on these pairs. It was computed once for each genus, since the lattices of a genus have isometric discriminant forms [Nik79, Cor. 1.9.4]. The listed classes are distinct, so their sum is at most the true right side of (4.1), which is at most the one computed with $R(q_F)$; since the two ends are equal, the list is complete and $R(q_F) = |O(q_F)|$. Independently of the program for $R$, every recorded value equals the closed formula of Lemma 2.1(ii), which the checks of §7.4 recompute from the rank of $F$ and the invariants $s$, $delta$ of its Gram matrix; so the right side of (4.1) is the true one. _In short_, $R(q)$ counts pairs $(x, y)$ with $q(x) = q(y) = 0$ and $b(x, y) = 1 slash 2$ and recurses on their orthogonal complement; this gives an upper bound for $|O(q)|$ and hence for the right side of (4.1), while the sum over the listed classes is a lower bound for the true right side; when the two agree, both the list and the value are certified. The closed formula of Lemma 2.1(ii) confirms the value independently. The numbers of classes and of extremal classes are in Tables 5 and 6.



*Type $(15, 17)$.* $F = sqrt(2) A_15^+$, where $A_15^+$ is the unimodular lattice of rank 15 without vectors of norm 1 [SPLAG, Table 16.7], and $K$ is the unique lattice of rank 17 with $s = 15$, $delta = 1$ and minimum at least 4. There are $"82,086"$ classes, $"20,502"$ of them extremal. The type $(17, 15)$ needs no list of its own: by §4.1.2 its pairs $(L, -sigma)$ are the inverse gluings $[phi^(-1)]$ of the classes of this list, and Table 5 counts $M$ for both orientations.

*Type $(14, 18)$.* $F = sqrt(2) (E_7^2)^+$, and $K$ is one of five lattices of rank 18 with $s = 14$, $delta = 1$. There are $"54,649"$ classes, $"16,590"$ of them extremal.

*Type $(12, 20)$.* $F = sqrt(2) D_12^+$, and $K$ is one of the 66 lattices of rank 20 with $s = 12$, $delta = 1$ and minimum at least 4 #calc("C5"). There are $"27,434"$ classes, $"13,798"$ of them extremal.

*Type $(16, 16)$.* There are ten components of rank 16, one for each $(s, delta)$ with $s in {8, 10, 12, 14}$ except $(8, 1)$, and three with $s = 16$: $sqrt(2) E_8^2$ and $sqrt(2) D_16^+$ ($delta = 0$) and $sqrt(2) (D_8^2)^+$ ($delta = 1$). Components with equal $(s, delta)$ give the eleven pairs of Table 6. In ten of them $F tilde.equiv K$, and $(L, sigma)$, $(L, -sigma)$ are both in the list (§4.1.2); for the pair $(sqrt(2) E_8^2, sqrt(2) D_16^+)$ the mass with both orientations counts this list twice. The pair $s = 16$, $delta = 1$ alone has $"222,208"$ classes before retaining only the extremal pairs.

*Type $(8, 24)$.* Here the gluing is trivial and the work is in the list of $K$. The only possible $F$ is $sqrt(2) E_8$ (Appendix B.2), with $"Aut" F = W(E_8)$, and $q_F$ is the quadratic form $x |-> (x, x) slash 2 mod 2$ on $E_8 slash 2 E_8$. It is well known that $W(E_8) -> O(q_F) = "GO"^+_8 (2)$ is onto with kernel ${plus.minus 1}$. By Proposition 4.1, $overline(G)_F = O(q_F)$ acts simply transitively on $"Iso"(q_F, -q_K)$, so for each $K$ there is exactly one double coset; and for every $v in "Aut" K$ there are exactly two $u in W(E_8)$ with $(u, v) in "Stab"(phi)$, so $|C(sigma)| = 2 |"Aut" K|$. Thus the pairs $(L, sigma)$ of type $(8, 24)$ correspond to the lattices $K$, and their contribution to $M$, counted for $sigma$ and $-sigma$, is $sum_K 1 slash |"Aut" K|$. By Proposition B.2, $L$ is extremal if and only if $min K >= 4$ and $min K^\# >= 2$. There are exactly $"14,482"$ such $K$ (Theorem B.3 and §B.3), with
$ sum_K 1 / (|"Aut" K|) = "4,770,461,726,538,527" / "792,723,456,000" approx "6,017.813" $
#calc("C6").

#ptable(
  columns: 6, align: (left, left, left, right, right, right),
  caption: [Table 5. Involution pairs: components, gluing classes, extremal classes and mass $M=sum 1 slash |"Stab"(phi)|$ over the extremal classes, with both orientations. For type $(8, 24)$ only the lattices $K$ that give extremal lattices are listed, so its classes are all extremal. Components of type $(16,16)$ are detailed in Table 6 #calc("C4", "C5", "C6", "C7").],
  table.header([type], [$F$], [$K$], [classes], [extremal], [$M$ (approx.)]),
  [$(15, 17)$], [$sqrt(2) A_15^+$; 1, $s = 15$, $delta = 1$], [1 class, $s = 15$, $delta = 1$], [$"82,086"$], [$"20,502"$], [$"9,664.170"$],
  [$(14, 18)$], [$sqrt(2) (E_7^2)^+$; 1, $s = 14$, $delta = 1$], [5 classes, $s = 14$, $delta = 1$], [$"54,649"$], [$"16,590"$], [$"6,876.012"$],
  [$(12, 20)$], [$sqrt(2) D_12^+$; 1, $s = 12$, $delta = 1$], [66 classes, $s = 12$, $delta = 1$], [$"27,434"$], [$"13,798"$], [$"5,507.962"$],
  [$(8, 24)$], [$sqrt(2) E_8$; 1, $s = 8$, $delta = 0$], [$"14,482"$ lattices, $s = 8$], [$"14,482"$], [$"14,482"$], [$"6,017.813"$],
  [$(16, 16)$], [10 classes, $8 <= s <= 16$], [same], [$"233,544"$], [$"46,805"$], [$"9,883.635"$],
)


#ptable(
  columns: (auto, 1fr, auto, auto, auto, auto, auto), align: (left, left, right, right, right, right, right),
  caption: [Table 6. The eleven pairs of components of type $(16, 16)$. In each pair $F tilde.equiv K$ except the one marked $F tilde.not K$, whose $M$ is counted for both orientations. "$C = 4$" is the number of extremal classes with $|"Stab"(phi)| = 4$. A dash in the column of components: the genus of $F$ contains exactly one lattice of minimum at least 4, and $F tilde.equiv K$ is that lattice.],
  table.header([$s$, $delta$], [components], [$|"Aut" F|$], [classes], [extremal], [$C = 4$], [$M$ (exact, extremal classes)]),
  [8, 0], [—], [$"89,181,388,800"$], [2], [2], [0], [$1 slash "22,830,435,532,800"$],
  [10, 0], [—], [$"21,139,292,160"$], [9], [4], [0], [$255 slash "3,758,096,384"$],
  [10, 1], [—], [$"743,178,240"$], [4], [4], [0], [$1 slash "6,193,152"$],
  [12, 0], [—], [$"195,689,447,424"$], [122], [33], [0], [$"7,139,633" slash "7,247,757,312"$],
  [12, 1], [—], [$"754,974,720"$], [72], [33], [0], [$"16,750,411" slash "251,658,240"$],
  [14, 0], [—], [$"53,271,016,243,200"$], [1,461], [135], [0], [$"151,662,361" slash "469,762,048"$],
  [14, 1], [—], [$"32,614,907,904"$], [6,108], [1,702], [667], [$"11,480,973,667" slash "50,331,648"$],
  [16, 0], [$sqrt(2) E_8^2$, $sqrt(2) E_8^2$], [$"970,864,271,032,320,000"$], [1,337], [81], [0], [$"5,800,534,177" slash "13,212,057,600"$],
  [16, 0], [$sqrt(2) E_8^2$, $sqrt(2) D_16^+$ ($F tilde.not K$)], [as in the adjacent rows], [1,080], [51], [0], [$2 dot "3,198,263" slash "5,160,960"$],
  [16, 0], [$sqrt(2) D_16^+$, $sqrt(2) D_16^+$], [$"685,597,979,049,984,000"$], [1,141], [74], [0], [$"925,264,055" slash "1,056,964,608"$],
  [16, 1], [$sqrt(2) (D_8^2)^+$, $sqrt(2) (D_8^2)^+$], [$"53,271,016,243,200"$], [$"222,208"$], [$"44,686"$], [$"34,940"$], [$"8,502,033,728,743" slash "880,803,840"$],
)

== 4.2 Elements of odd prime order
We use the notation of §3.3 and §4.1.1 with $h$ of odd prime order $p$ and type $p - (z, f) - s$.

#thm[
*Lemma 4.2.* If $min Z >= 4$, then $min Z^\# >= 2 (p - 1) slash p$.
]
#proofbox[
Let $0 != y in Z^\#$ and $beta = (y, y)$. Since $(1 - h^i) y in Z$ is non-zero for $1 <= i < p$ [Neb13b, proof of Lemma 4.2], its norm is at least 4. As $sum_(i = 0)^(p - 1) y h^i = 0$, we have $sum_(i = 1)^(p - 1) (y, y h^i) = -beta$, and $4 (p - 1) <= sum_(i = 1)^(p - 1) (y - y h^i, y - y h^i) = 2 (p - 1) beta + 2 beta = 2 p beta$. #qed
]

#thm[
*Lemma 4.3.* Let $min F >= 4$ and $min Z >= 4$. Then $L_phi$ has a vector of norm 2 if and only if there are $a in F^\# slash F$ containing a vector of norm $2 slash p$ and $b in Z^\# slash Z$ containing a vector of norm $2 - 2 slash p$ with $phi(a) = b$.
]
#proofbox[
A vector of norm 2 is $x + y$ with $x in F^\#$, $y in Z^\#$ and $(x, x) + (y, y) = 2$; neither is 0 by the assumption on the minima. Since $F^\# slash F$ has exponent $p$ and $F$ is even, $(x, x) in (2 slash p) ZZ$, and $(y, y) >= 2 - 2 slash p$ by Lemma 4.2. Hence $(x, x) = 2 slash p$ and $(y, y) = 2 - 2 slash p$. Conversely, if $x in a$ has norm $2 slash p$ and $y in b$ has norm $2 - 2 slash p$, then $x + y in L_phi$ has norm 2. #qed
]

In particular, if $F^\#$ has no vector of norm $2 slash p$ or $Z^\#$ has none of norm $2 - 2 slash p$, every gluing is extremal.

#thm[
*Lemma 4.4.* Let $pi.alt = 1 - zeta_p$. There is a Hermitian form $tilde(H)$ on $Z$ over $ZZ[zeta_p]$, a rescaling of the one with trace $(dot, dot)$, such that $Z^* := {x : tilde(H)(x, Z) subset ZZ[zeta_p]}$ equals $pi.alt Z^\#$, $Z^* subset Z subset pi.alt^(-1) Z^*$, and $Z slash Z^*$ carries a non-degenerate alternating form over $FF_p$. Hence $z - s$ is even, and the lattices $U$ with $Z^* subset U subset Z$ corresponding to maximal isotropic subspaces are unimodular for $tilde(H)$, of rank $z$ and index $p^((z - s) slash 2)$ in $Z$.
]
#proofbox[
Let $KK = QQ(zeta_p)$, $cal(O) = ZZ[zeta_p]$, on which $h$ acts as $zeta_p$, and $KK^+$ its real subfield. The different of $KK slash QQ$ is $frak(D) = pi.alt^(p - 2) cal(O)$, and, since $zeta_p equiv 1 med (mod pi.alt)$, complex conjugation acts trivially on $cal(O) slash pi.alt = FF_p$. Let $H$ be the Hermitian form with $(x, y) = "Tr"_(KK slash QQ) H (x, y)$, so that $Z^\# = {x : H (x, Z) subset frak(D)^(-1)}$. Put $nu = pi.alt overline(pi.alt)$, a totally positive element of $KK^+$ with $nu cal(O) = pi.alt^2 cal(O)$, let $beta = nu^((p - 3) slash 2)$, so that $beta cal(O) = pi.alt^(p - 3) cal(O)$, and $tilde(H) = beta H$, again Hermitian since $beta$ is real. Then $tilde(H)(x, Z) subset cal(O)$ if and only if $H (x, Z) subset pi.alt frak(D)^(-1)$, so $Z^* = pi.alt Z^\#$. Since $pi.alt$ acts as $1 - h$ and $(1 - h) Z^\# subset Z$ [Neb13b, proof of Lemma 4.2], $Z^* subset Z$; and $H (Z, Z) subset frak(D)^(-1)$ gives $tilde(H)(Z, Z) subset pi.alt^(-1) cal(O)$, so $pi.alt Z subset Z^*$. The index is $[Z : Z^*] = [Z^\# : pi.alt Z^\#] slash [Z^\# : Z] = p^(z - s)$.

On $Z slash Z^*$ put $a(x, y) = pi.alt tilde(H)(x, y) mod pi.alt in FF_p$. It is well defined since $tilde(H)(Z^*, Z) subset cal(O)$, and non-degenerate since $pi.alt tilde(H)(x, Z) subset pi.alt cal(O)$ means $x in Z^*$. As $overline(pi.alt) = -zeta_p^(-1) pi.alt equiv -pi.alt med (mod pi.alt^2)$, $a(y, x) = -a(x, y)$. Finally $tilde(H)(x, x) in KK^+ inter pi.alt^(-1) cal(O)$, and elements of $KK^+$ have even $pi.alt$-adic valuation because $pi.alt$ is ramified in $KK slash KK^+$; so $tilde(H)(x, x) in cal(O)$ and $a(x, x) = 0$. Thus $a$ is a non-degenerate alternating form, and $z - s$ is even.

The dual of $Z^*$ with respect to $tilde(H)$ is $Z$ (because $tilde(H)(x, pi.alt Z^\#) subset cal(O)$ if and only if $H (x, Z^\#) subset frak(D)^(-1)$). Hence for $Z^* subset U subset Z$ the $tilde(H)$-dual $U^*$ satisfies $U^* subset Z$, and $U^* slash Z^*$ is the orthogonal of $U slash Z^*$ for $a$. If $U slash Z^*$ is maximal isotropic, it equals its orthogonal, so $U^* = U$. #qed
]

*Roots of unimodular sublattices.* Let $U$ be as in Lemma 4.4, with the form $tilde(H)$ of Lemma 4.4. If $x in U$ is a root, then $cal(O) x$ is an orthogonal summand of $U$ isometric to the root lattice $A_(p-1)$ with $h$ acting as a Coxeter element. Indeed, the roots $x h^i$ span $cal(O) x$, a root lattice of rank $p - 1$ with an automorphism $h$ of order $p$ without fixed vectors; $h$ cannot permute $p$ components, so it preserves the component containing $x$, which therefore is all of $cal(O) x$; and the only irreducible root lattice of rank $p - 1$ with an automorphism of order $p$ is $A_(p-1)$, on which all such automorphisms are conjugate (they are $p$-cycles in $W(A_(p-1)) times {plus.minus 1}$). For $A_(p-1) = pi.alt cal(O)$ with $(a, b) = "Tr"_(KK slash QQ)(a overline(b)) slash p$, the form of Lemma 4.4 is $tilde(H)(a, b) = (pi.alt overline(pi.alt))^((p-3) slash 2) a overline(b) slash p$, so $tilde(H)(pi.alt, pi.alt) = (pi.alt overline(pi.alt))^((p-1) slash 2) slash p$ is a unit. So $tilde(H)(x, x)$ is a unit, and since $U$ is unimodular, $cal(O) x$ splits off. A root $y$ of $U$ is $y = a + b$ with $a in cal(O) x$, $b perp cal(O) x$, and $(a, a) + (b, b) = 2$ with both terms in $2 ZZ$, so $y in cal(O) x$ or $y perp x$. Hence
#tagged("4.2", $ U = A_(p-1)^k perp U', quad |"U"(U)| = (2 p)^k k! |"U"(U')|, $)
with $U'$ rootless and determined by $U$. Summing (4.2) over a genus expresses its mass through the rootless masses of smaller rank; this is how the rootless masses used in Proposition 3.3 (reason 3) and below are obtained (Appendix D).

*Components.* For $h$ of type $p - (z, f) - s$, $F$ is an even $p$-elementary lattice of rank $f$ and determinant $p^s$, and $Z$ is a Hermitian $cal(O)$-lattice of rank $z$ on which $h$ acts as $zeta_p$; both have minimum at least 4. For odd $p$ the genus of $F$ is determined by $(f, s)$: a genus is determined by its rank, signature and discriminant form [Nik79, Cor. 1.9.4], and the discriminant form on $FF_p^s$ is fixed by the existence condition of [SPLAG, Ch. 15, Thm 13(a)]; and the genus of $Z$ is determined by $(z, s)$ (Proposition D.1). Appendix B.4 lists the components and proves the lists complete (Table B4): the lattices $F$ by a classification of their genus by Kneser neighbours, completed when $sum 1 slash |"Aut" F|$ reaches the mass of the genus, or, for two genera, by the rootless mass of the genus; the lattices $Z$ through their unimodular sublattices $U$, whose lists are complete by the rootless Hermitian mass, and then $Z supset U$ by Lemma 4.4, step by step in $s$, with a mass identity at each step.

*Gluings.* For each pair $(F, Z)$, the double cosets of Proposition 4.1(ii) were computed with GAP as double cosets in $O(q_F)$, a group whose order $|O(q_F)|$ is given by Lemma 2.1(i) (the recorded orders agree with it, §7.4), realized as a permutation group on $FF_p^s$, or by a canonical form when GAP did not finish #calc("C8", "C9")\; in every case identity (4.1) holds as an equality of fractions. Extremality is decided by Lemma 4.3. If $F^\#$ has no vector of norm $2 slash p$ or $Z^\#$ has none of norm $2 - 2 slash p$, every gluing is extremal; this holds for 8 of the 16 types with $F != 0$.



*Order 3.* The ten types have $"6,617"$ extremal pairs; $"5,383"$ of them have $|C(h)| = 6$. For $3 - (12, 8) - 8$, the largest type, $F = sqrt(3) E_8$ with $overline(G)_F = W(E_8)$, $Z$ runs over 51 classes, and the double cosets are the orbits of $"U"(Z)$ on the $"56,862"$ cosets of $W(E_8)$ in $"GO"^+_8 (3)$; there are $"15,423"$ pairs, of which $"5,872"$ are extremal. For $3 - (6, 20) - 6$, $Z$ is the Coxeter–Todd lattice $K_12$ with its Eisenstein structure and $overline(G)_Z$ has index 2 in $O(q_Z) tilde.equiv "GO"^-_6 (3)$, so each $F$ gives one or two double cosets; GAP shows that two occur for exactly one of the 22 lattices $F$ with $F^\#$ free of vectors of norm $2 slash 3$, which gives 23 extremal pairs. For the other two lattices $F$, whose duals contain vectors of norm $2 slash 3$, no gluing is extremal: a class $a in F^\# slash F$ containing such a vector has $q_F (a) = 2 slash 3$, so $q_Z (phi(a)) = 4 slash 3$ for every $phi$, and every class of $K_12^\# slash K_12$ with $q = 4 slash 3$ contains a vector of norm $4 slash 3$, so Lemma 4.3 applies. Indeed, two distinct vectors of norm $4 slash 3$ in one class differ by a non-zero vector of $K_12$, of norm at least 4, so their inner product is at most $-2 slash 3$; the sum of $k$ of them has norm at most $4 k slash 3 - 2 k (k - 1) slash 3$, so a class contains at most three of them. Hence the $756$ vectors of norm $4 slash 3$ of $K_12^\#$ meet at least $252$ classes, all with $q = 4 slash 3$, while by Milgram's formula ($|K_12^\# slash K_12| = 3^6$, signature 12) exactly $252$ classes have $q = 4 slash 3$ #calc("C10"). For $3 - (16, 0) - 0$, $F = 0$ and $L = Z$ is a rootless Hermitian lattice of rank 16 over $ZZ[zeta_3]$; there are exactly eleven such lattices [Hen13, Thm 1, Table 4, Rem. 3(a)]. We constructed eleven such lattices without the Gram matrices of [Hen13], by neighbour methods and from codes, computed $C(h)="U"(Z)$ ourselves, and checked that the orders $|"U"(Z)|$ agree with [Hen13, Table 4] and that $sum 1 slash |"U"(Z)|$ equals the rootless mass of the genus (Appendix D).

*Orders 5, 7, 17, 31.* The nine types have $218$ extremal pairs. For $5 - (7, 4) - 3$, $F = sqrt(5) A_4^\#$ and the classes of $F^\# slash F$ containing a vector of norm $2 slash 5$ are exactly those with $q_F = 2 slash 5$, so $L_phi$ is extremal if and only if $Z^\#$ has no vector of norm $8 slash 5$, independently of $phi$; 13 of the 23 lattices $Z$ qualify, and since $"Aut" F -> O(q_F)$ is onto, each gives one pair with $C(h) tilde.equiv "U"(Z)$. For $5 - (6, 8) - 6$, $F = F_1 perp F_1$ with $F_1 = sqrt(5) A_4^\#$, the classes of $phi$ correspond to unordered pairs ${V_1, V_2}$ of subspaces with $D(Z) = V_1 perp V_2$, and counting the $"U"(Z)$-orbits of the good pairs over the four lattices $Z$ gives 155 pairs. For $5 - (8, 0) - 0$, $F = 0$ and $L = Z$; the nine classes were found by $pi.alt$-neighbours, which preserve the $ZZ[zeta_5]$-structure, and the list is complete because the sum of $1 slash |"U"(Z)|$ equals the rootless mass of the genus (Appendix D). For $17 - (2, 0) - 0$, $L = Z$ is a Hermitian $ZZ[zeta_17]$-lattice of rank 2; its genus consists of seven Hermitian lattices, six of which have extremal trace form, and these six give exactly four lattices up to isometry [Jür15, Satz 3.6.9 and its proof]; independently, the sum over these four lattices of the number of elements of order 17 divided by $|"Aut" L|$ equals the rootless mass $1 slash 12$ of the genus #calc("C9")\; there are six pairs. For $31 - (1, 2) - 1$, $F$ is the binary form $mat(4, 1; 1, 8)$ with $"Aut" F = {plus.minus 1} = O(q_F)$, $Z$ runs over the rank-1 lattices of its genus other than $A_30$; every lattice of rank 1 has unitary group $chevron.l -zeta_31 chevron.r$ of order 62 and the genus has mass $9 slash 62$ (§D.1), so it has 9 classes (the relative class number $h^-(QQ(zeta_31)) = 9$ [Sch98, p. 1226]); so the eight rootless ones have total weight $4 slash 31$ (Appendix D). The gluing is unique; this gives eight pairs with $|C(h)| = 62$.

#rem[
*Example (identity (4.1) for the type $3 - (9, 14) - 7$).* This is a case in which identity (4.1) can be followed by hand; the example is not used in the proofs. By Table B4, the genus $(14, 7)$ of $F$ and the genus $(9, 7)$ of $Z$ each contain exactly one lattice of minimum at least 4, so there is one pair of components. For it, $|"Aut" F| = "8,491,392" = 2^7 3^6 dot 7 dot 13$ and $|"U"(Z)| = "22,674,816" = 2^7 3^(11)$. The form $Q$ of Table 3 on $F^\# slash F tilde.equiv FF_3^7$ is non-degenerate of odd dimension 7. Hence, by Lemma 2.1(i), $|O(q_F)| = |"GO"_7 (3)| = 2 dot 3^9 (3^2 - 1)(3^4 - 1)(3^6 - 1) = "18,341,406,720"$. The right side of (4.1) is therefore
$ ("18,341,406,720") / ("8,491,392" dot "22,674,816") = 5 / ("52,488") . $
The computed list has two double cosets, with $|"Stab"(phi)| = "11,664"$ and $"104,976"$, and
$ 1 / ("11,664") + 1 / ("104,976") = 10 / ("104,976") = 5 / ("52,488") . $
So the left side already equals the right side. As every term is positive, there is no third double coset, whatever method found the two. Both gluings are extremal #calc("C8"). The two lattices are `32.23328.1` and `32.839808.1` of [D], with $|"Aut" L| = 2^5 3^6$ and $2^7 3^8$; the second paper of this repository describes them (its Theorem 6.1).
]



== 4.3 Gaussian pairs and their rootless mass
Let $cal(O)=ZZ[i]$. By §2.4, a pair $(L,g)$ with $g^2=-1$ is the same as an $cal(O)$-lattice $L$ of rank $r=16$ with a Hermitian form $H$ whose trace form $"Tr" H$ is even unimodular; isometries of pairs are the $cal(O)$-linear isometries of Hermitian lattices, so $C(g) = "U"(L, H)$, the unitary group. $L$ is extremal exactly when $(L, H)$ has no _root_, a vector with $"Tr" H(x, x) = 2$. Hence
$ M_2 = sum_((L, H)) 1 / (|"U"(L, H)|), $
the sum over isometry classes of rootless Hermitian $cal(O)$-lattices of rank $r$ with even unimodular trace form: the _rootless mass_.

#thm[
*Theorem 4.5.* $M_2 = "4,128,748,922,231,417" slash "120,324,096,000"$.
]
The proof is in Appendix C. It follows King's method for $ZZ$-lattices [Kin03]: the rootless mass is a finite alternating combination (C.2) of the masses $A_r (T)$ of embeddings of root lattices $T$, with coefficients from the classification of root lattices over $ZZ[i]$ and their root sublattices; each $A_r (T)$ is given by a Siegel-type formula (C.3) as a product of local densities; and the local densities at the ramified prime 2 reduce to counting on $L slash 2 L$ (Proposition C.7, Theorem C.8). The computation is in exact rational arithmetic #calc("C11")\; the implementation is among the assumptions of Table 10.

#thm[
*Proposition 4.6.* The list #calc("C12") consists of $"76,742"$ rootless Hermitian $ZZ[i]$-lattices of rank 16, no two of which are isometric or complex conjugate to each other. With $f = 2$ for a lattice not isometric to its complex conjugate and $f = 1$ otherwise, the sum of $f slash |"U"|$ over the list equals the value of Theorem 4.5. Hence the list contains every rootless Hermitian $ZZ[i]$-lattice of rank 16 up to isometry and complex conjugation.
]
A pair $(L, g)$ with $g^2 = -1$ and its pair $(L, -g)$ are the Hermitian lattice $(L, H)$ and its complex conjugate, so (§2.4) $M_2 = sum f slash |"U"|$ over the lattices up to isometry and complex conjugation; these are the Hermitian entries of order 4 of §6.1.

_Strategy of the proof._ The search that produced the list does not prove it complete. Completeness follows from Lemma 7.1, which needs three facts: the classes of the list are pairwise distinct, each weight $f slash |"U"|$ is correct, and the weights add up to $M_2$, which Theorem 4.5 computes independently of the list. A missing class would add a positive weight, so none is missing. The proof checks the three facts in this order.

*Proof.* The lattices were collected from the Hermitian entries of order 4 of the lattices met in the lists of §4.1 and §4.2 and in $S_"gen"$, and from random walks by neighbours in the genus #calc("C12"). Two lattices were kept apart only if they differ in a canonical form or in the invariant $D$ of the pair $(L, plus.minus g)$, both invariant under isometry and complex conjugation (§7.2). This proves the first statement. The weights $f slash |"U"|$ were obtained by three methods: from the automorphism group of a canonical form, by a backtrack search for the isometries commuting with $g$ up to sign (it returns their group $N_plus.minus$, and $f slash |"U"| = 2 slash |N_plus.minus|$), and from $"Aut" L$ computed with bliss. Every lattice has a weight from at least one method. For $"11,421"$ lattices two methods were applied (bliss and the search for $"11,144"$, bliss and the canonical form for $277$), and they agree. Of the others, $"62,449"$ have their weight from the canonical form alone: there a canonically selected set of Gaussian lines with distinct colours generates $L$. Each line is preserved by $"U"$; connectedness and the phase-normalization argument of §7.2 force $"U" = {plus.minus 1, plus.minus i}$. For the discrete classes in the final list, the two normalized Gram hashes differ, which excludes an isometry to the complex conjugate; $"2,547"$ have their weight from the backtrack search alone and $325$ from bliss alone #calc("C12"). The sum of the weights equals the value of Theorem 4.5 as a fraction #calc("C12"). Since Theorem 4.5 is proved in Appendix C, Lemma 7.1 shows that the list is complete. #qed

The list supplies a lower bound for the Gaussian rootless mass independently of the mass formula. Equality establishes completeness under the validity, distinctness and group-weight conditions just proved.



= 5. Proof of Theorem A <s5>

#policy-box[
*What this section does.* We reduce the total to the classification (§5.1) and list the statements on which the proof rests, with their proofs, assumptions and computations (§5.2, Table 7).
]

#rem[
*Points that are easy to misread (the count).* These remarks add nothing to the proofs.
+ _(5.1) is not an independent check of Theorem B._ It is an exact identity, true for any finite set of lattices whose automorphism groups contain $plus.minus 1$. It turns the classification of $S$ into the number $N$, and it is only as reliable as the list of $S$ and the orders $|"Aut" L|$ it reads (§1.4).
+ _The lattices with $"Aut" L = {plus.minus 1}$ are counted, not constructed._ Their number $"10,900,939"$ (Corollary C) is $N - |S|$. No representative of them is computed in this paper.
+ _The auxiliary count is a check._ The thirteen-term computation of [Aux] reaches the same $N$ by another route, but it shares the component lists with the main proof and is not a premise (§5.2).
+ _A "trivial" entry is not a trivial element._ An involution entry is trivial when $C(sigma) = {plus.minus 1, plus.minus sigma}$ and $sigma tilde.not -sigma$; a Hermitian entry of order 4 is trivial when $N_G (chevron.l g chevron.r) = chevron.l g chevron.r$ (§6.1, Lemma 6.1). Such an entry still comes from an automorphism other than $plus.minus 1$, and its lattice still belongs to $S$.
+ _Identifiers and keys._ An identifier names the class of an entry in one of the lists of §4, and some rows of the records have none (§6.3). A key is an isometry invariant used to separate and sort (§7.2, §6.8). Equal keys prove isomorphism only in the situation of Lemma 7.2.
+ _Labels._ In a label `32.`$|"Aut" L|$`.`$i$ (§6.8), the index $i$ follows a fixed rule for sorting and has no further meaning.
]

== 5.1 Reduction to the classification

For every $L in cal(L)$, write $G_L="Aut" L$. Since $G_L$ contains $plus.minus 1$,
#tagged("5.1", $ N=2m+sum_(L in S)(1-2 slash |G_L|). $)
Indeed, for each class split $1$ as $2 slash |G_L|+(1-2 slash |G_L|)$ and sum. The second summand vanishes exactly for the classes outside $S$.

Theorem B, proved in §6 from the complete structured-pair lists, supplies a non-redundant list of $S$ and the orders $|G_L|$. Substitution in (5.1), in exact rational arithmetic, gives $N="11,044,268"$. The computation reads only King's mass and the saved classified-lattice table #calc("C13"). This proves Theorem A conditional on precisely the classification and mass inputs described above. The proof of Theorem B does not use this numerical value.

#rem[
*Example (the method in rank 24).* This example is not used in the proofs. In rank 24 the bound of §2.1 is again $min L <= 4$, and the only extremal lattice is the Leech lattice $Lambda$, with $"Aut" Lambda = "Co"_0$ [SPLAG, Table 16.1]. The auxiliary note runs the tools of §4 and of Appendices C and D in rank 24 and compares every value with the conjugacy classes of $"Co"_0$ [Aux, §5.6, Table 7c]. King's rootless mass of rank 24 is $1 slash |"Co"_0|$, with $|"Co"_0| = "8,315,553,613,086,720,000"$. For involutions of type $(16, 8)$ the list of §4.1 has one double coset; it is extremal, and it is $Lambda$, with $|"Stab"(phi)| = |C(sigma)| = "178,362,777,600"$. For type $(12, 12)$ the list has six double cosets, and only one is extremal: $Lambda$ again, with $sigma tilde -sigma$. Since $Lambda$ is the only extremal lattice, $S = {Lambda}$ in rank 24, and (5.1) gives
$ 2 slash |"Co"_0| + (1 - 2 slash |"Co"_0|) = 1 . $
In rank 24 the mass term is negligible, and the count comes from $S$. In rank 32 it is the other way round: $2 m approx "10,968,923.3"$ is 99.3% of $N$, and the $"143,329"$ lattices of $S$ add only about $"75,344.7"$. The comparison in rank 24 checks the programs; it does not check the values of rank 32.
]

== 5.2 Dependencies of the main proof

The classification-based proof is assembled from the statements of Table 7. Each row states its output, the equivalence relation and weight under which a list is complete, how it is proved, which software assumptions of Table 10 it uses, and the computations of Table 11 that support it. Long proofs and finite calculations are in the appendices; the main text cites the statements.

Figure 2 draws the same dependencies as a diagram.

#figure(kind: image, numbering: none, gap: 0.8em,
  {
  _x-diagram(440, 390, (
    det: (0, 0, 100, "hand", [*detection* \ Proposition 3.1 \ row 9]),
    types: (113, 0, 100, "both", [*possible types* \ Propositions 3.2, 3.3 \ row 2]),
    keys: (339, 58, 100, "hand", [*keys* \ §7.2, Lemmas 7.1, 7.2 \ row 8]),
    inv: (0, 58, 100, "both", [*involution components* \ App. B.2–B.3, (B.1) \ row 3]),
    odd: (113, 58, 100, "finite", [*odd-order components* \ App. B.4, D, (B.2) \ row 4]),
    mass2: (226, 58, 100, "both", [*Gaussian mass* $M_2$ \ Theorem 4.5, App. C \ row 6]),
    glue: (56, 116, 100, "both", [*gluing classes* \ Proposition 4.1, (4.1) \ row 5]),
    glist: (226, 116, 100, "both", [*Gaussian list* \ Proposition 4.6 \ row 7]),
    groups: (339, 116, 100, "hand", [*finite groups* \ Lemmas 6.1, 6.2, \ Props. 6.4, 6.5; row 9]),
    gen: (56, 174, 100, "both", [*lattices of* $S_"gen"$ \ §6.3, Lemma 6.3 \ row 10]),
    sub: (282, 174, 100, "hand", [*subtraction* \ Lemma 6.6: $S_(C_2^2)$, $S_"CQ"$ \ row 11]),
    theorem-b: (170, 232, 100, "result", [*Theorem B* \ the classification of $S$]),
    king: (339, 232, 100, "cited", [*King's mass* $m$ \ §2.2, (2.1) \ row 1]),
    agg: (226, 290, 140, "both", [*identity (5.1)* \ $N = 2 m + sum_(L in S) (1 - 2 slash |"Aut" L|)$ \ row 12]),
    theorem-a: (246, 348, 100, "result", [*Theorem A* \ and Corollary C]),
  ), (
    ("det", "r", "types", "l"), ("types", ("b", 0.3), "inv", ("t", 0.6)), ("types", "b", "odd", "t"),
    ("inv", ("b", 0.7), "glue", ("t", 0.3)), ("odd", ("b", 0.3), "glue", ("t", 0.7)),
    ("mass2", "b", "glist", "t"), ("keys", ("b", 0.15), "glist", ("r", 0.35)),
    ("glue", "b", "gen", "t"), ("glue", ("r", 0.8), "sub", ("t", 0.15)),
    ("glist", ("b", 0.7), "sub", ("t", 0.4)), ("groups", ("b", 0.4), "sub", ("t", 0.8)),
    ("gen", "r", "sub", "l"),
    ("gen", ("b", 0.7), "theorem-b", ("t", 0.3)), ("sub", ("b", 0.3), "theorem-b", ("t", 0.7)),
    ("theorem-b", ("b", 0.7), "agg", ("t", 0.3)), ("king", ("b", 0.5), "agg", ("t", 0.8)),
    ("agg", "b", "theorem-a", "t"),
  ))
  v(4pt)
  [#_x-key("hand")[written proof] #h(10pt) #_x-key("finite")[exact finite computation] #h(10pt) #_x-key("both")[both] #h(10pt) #_x-key("cited")[cited result]]
  },
  caption: [Figure 2. The dependencies of the proofs of Theorems A and B, drawn from Table 7 ("row $k$" is its row). An arrow means "is used by". Colours give the kind of proof, as in the column "proof" of Table 7; the computations behind each box are listed there. Proposition 3.1 also gives the three parts of $S$ (§6.1), and the keys of §7.2 also separate the lattices of $S_"gen"$; these two uses are not drawn. Theorem B does not use King's mass or the value of $N$.],
)

#[
#set text(size: 0.86em)
#ptable(
  columns: (1.15fr, 1.1fr, 1.45fr, 1.1fr, 0.62fr), align: (left, left, left, left, left),
  caption: [Table 7. The statements used by the proofs of Theorems A and B. "Hand" means a written proof in this paper; "cited" a published result; "finite" an exact finite computation. Theorem B uses rows 2–11, Theorem A in addition rows 1 and 12.],
  table.header([statement], [classes, weights], [proof], [assumptions (Table 10)], [computations]),
  [1. King's rootless mass $m$ (§2.2)], [$ZZ$-isometry; $1 slash |"Aut" L|$], [cited [Kin03, Table 1]; recomputation of the same inclusion–exclusion], [King's table; our code], [C14, C15],
  [2. Possible types (Propositions 3.2, 3.3)], [—], [hand; cited [Nik79], [SPLAG], [CE03], [Neb13b]; six finite LP certificates; rootless Hermitian masses of Appendix D (reason 3)], [dimensions of spaces of modular forms; theta coefficients; our implementation of the mass formulas of Appendix D], [C2, C3],
  [3. Involution components (Appendix B.2–B.3)], [$ZZ$-isometry; $1 slash |"Aut"|$], [hand (Lemma B.1, Proposition B.2, Theorem B.3): complete by an enumeration that reaches every lattice, except the lattices $K$ of type $(8, 24)$ for $A_2^(12)$, which are complete by identity (B.1)], [`qfminim`, `qfauto`, `qfisom`, bliss, GAP (the orbits of §B.3)], [C5, C6, C7],
  [4. Components of odd prime order (Appendix B.4, D)], [$F$: $ZZ$-isometry; $(Z, h|_Z)$: isometries commuting with $h$, not identified under Galois; $1 slash |"Aut" F|$, $1 slash |"U"(Z)|$], [finite: equality with genus and rootless masses; the 2-adic and 9-torsion counts of §D.4; cited [Hen13], [Jür15]], [SageMath masses; our code for Appendix D; `qfauto`], [C2, C8, C9, C16, C17, C18],
  [5. Gluing classes (Proposition 4.1, §§4.1.3, 4.2)], [pairs $(L, h)$ up to isometry; orientation as in §4.1.2; $1 slash |"Stab"(phi)|$], [hand (Proposition 4.1, Lemmas 4.2–4.4); finite: identity (4.1) as an equality of fractions for each pair of components], [GAP double cosets; canonical form of gluings; `qfminim` for extremality], [C4, C6, C8, C9],
  [6. Gaussian rootless mass $M_2$ (Theorem 4.5, Appendix C)], [$ZZ[i]$-linear isometry; $1 slash |"U"|$; conjugates not identified], [hand (C.1–C.10); finite: exact sum over $"12,934"$ root types], [our code for Appendix C; GAP for the Weyl group classes], [C11, C19, C20],
  [7. Gaussian list (Proposition 4.6)], [isometry and complex conjugation; $f slash |"U"|$], [Lemma 7.1 with row 6; keys of §7.2], [bliss; backtrack search; canonical forms], [C12],
  [8. Keys (Lemmas 7.1, 7.2, §7.2)], [one relation per key], [hand: invariance of each key], [bliss canonical labelling], [—],
  [9. Finite groups (Proposition 3.1, Lemmas 6.1, 6.2, Propositions 6.4, 6.5)], [—], [hand; cited [Gor80]], [none], [C21 (check)],
  [10. The lattices of $S_"gen"$ (§6.3, Lemma 6.3)], [$ZZ$-isometry; full $|"Aut" L|$], [shared exact identifiers and integer isometries for equality; Lemma 6.3 (exact types) and minimum-shell profiles for distinctness], [bliss and GAP for full groups and classes; `qfminim`; exactness of identifiers (Table 10, our code)], [C1, C22, C23, C24, C25],
  [11. Subtraction (Lemma 6.6)], [entries of each lattice], [hand, from rows 3–10], [the records of §6.3 list all involution and Hermitian entries of the full group], [C25, C26, C27, C28, C29],
  [12. Identity (5.1)], [—], [hand; finite exact sum], [—], [C13],
)
]

The integer obtained in (5.1), the auxiliary thirteen-term sum and the distribution of group orders do not replace any row of Table 7.

The auxiliary thirteen-term computation of [Aux] uses rows 1–6 but neither row 10 nor row 11: it never computes a full automorphism group of rank 32. Its agreement with (5.1) is therefore a global check by a different arithmetic route of the full-group orders and the identifications of rows 10–11. It does not test the rows shared by both routes, and it cannot detect errors that cancel. It also relies on its own inputs $e(h)$ and $M_3, M_4, M_5$. It is not a premise of the proof.



= 6. Proof of Theorem B <s6>

#policy-box[
*What this section does.* We group the pairs of §4 lattice by lattice. §6.1 defines the three parts. §6.2 proves Lemma 6.2: a lattice whose entries are all trivial has $|"Aut" L| = 4$. §6.3 identifies the lattices of $S_"gen"$ and proves them distinct (Lemma 6.3). §6.4 and §6.5 prove Propositions 6.4 and 6.5 on the groups of $S_(C_2^2)$ and $S_"CQ"$. §6.6 removes the entries of $S_"gen"$ (Lemma 6.6) and counts both parts. §6.7 obtains the total from the classified lattices and separates it from the auxiliary checks. §6.8 describes the numbering of the database.
]

== 6.1 Entries and the three parts
Fix $L in S$ and put $G="Aut" L$. The detection argument of §1.3 gives the structures needed for the classification. We use the following _entries_:
- an _odd entry_ is a conjugacy class $[h]$ with $h$ of odd prime order;
- an _involution entry_ is the set ${[sigma],[-sigma]}$ for $sigma!=plus.minus 1$, a single class if $sigma tilde -sigma$;
- a _Hermitian entry_ of order $2^a$ is the set ${[g],[-g]}$ with $g^(2^(a-1))=-1$.
For order 4 put $"U"=C(g)$ and $f=2$ or 1 according as $g$ is not or is conjugate to $-g$. Its Gaussian mass weight is $f slash |"U"|$. No assignment count $e(h)$ is needed to define these entries or the three parts below. Auxiliary assignment statistics [Aux] may still be computed inside the full groups for invariant signatures and consistency checks.

In terms of the gluings of §4.1, an involution entry of type $(d, 32 - d)$ with $d != 16$, or with $d = 16$ and $F tilde.not K$, is one gluing class $[phi]$ (seen from $sigma$; from $-sigma$ it is the class of $phi^(-1)$ in the list of type $(32 - d, d)$). For $d = 16$ and $F tilde.equiv K$ it is the pair ${[phi], [phi^(-1)]}$ of classes in the same list, which coincide exactly when $sigma tilde -sigma$.

An entry is _trivial_ if it is an involution entry with $|C(sigma)| = 4$ and $sigma tilde.not -sigma$, or a Hermitian entry of order 4 with weight $1 slash 2$. A trivial involution entry has $C(sigma)={plus.minus 1,plus.minus sigma}$. All other entries are _non-trivial_.

#thm[
*Lemma 6.1.* A Hermitian entry of order 4 has weight $1 slash 2$ if and only if $N_G (chevron.l g chevron.r) = chevron.l g chevron.r$.
]
#proofbox[
We have $chevron.l g chevron.r subset.eq "U" subset.eq N_g := N_G (chevron.l g chevron.r)$, and $N_g slash "U"$ embeds in $"Aut" chevron.l g chevron.r tilde.equiv C_2$, with $[N_g : "U"] = 2$ exactly when $g$ is conjugate to $g^(-1) = -g$, that is when $f = 1$. Hence $f slash |"U"| = 1 slash 2$ if and only if $f = 2$ and $|"U"| = 4$, if and only if $N_g = "U" = chevron.l g chevron.r$. #qed
]

We define the three parts of $S$; the subscripts name the automorphism groups that they turn out to have:
- $S_"gen"$ ("general"): $L$ has an odd entry or a non-trivial involution entry, that is, an element of odd prime order or an involution $sigma != plus.minus 1$ with $|C(sigma)| > 4$ or $sigma tilde -sigma$;
- $S_(C_2^2)$: $L in.not S_"gen"$ and $L$ has an involution $sigma != plus.minus 1$; by Proposition 6.4, $"Aut" L tilde.equiv C_2 times C_2$;
- $S_"CQ"$: $L$ has neither an element of odd prime order nor an involution other than $-1$, so $"Aut" L$ is a 2-group with $-1$ as its only involution, cyclic or generalized quaternion (Proposition 6.5).
These sets are disjoint by definition and exhaust $S$ by Proposition 3.1. If there is no odd entry or involution entry, the group is a 2-group whose only involution is $-1$. In $S_"CQ"$ the condition that $"Aut" L$ is a 2-group matters: a group such as $C_6$ has $-1$ as its only involution but contains an element of order 3, so its lattices lie in $S_"gen"$. (In the database, $"2,759"$ lattices of $S_"gen"$ have $-1$ as their only involution, $"2,308"$ of them with $"Aut" L tilde.equiv C_6$; the others have groups such as $C_12$, $"Dic"_3$, $"SL"(2, 3)$, $"SL"(2, 5)$ and $"SL"(2, 17)$ #calc("C30").)

== 6.2 Lattices with only trivial entries
#thm[
*Lemma 6.2.* Let $L in S$ and $G = "Aut" L$.
+ $|G| = 4$ if and only if all entries of $L$ are trivial. In that case $L$ has exactly one entry.
+ If $G$ is a 2-group of order at least 8, then $L$ has no trivial entry.
+ $L$ has at most one trivial entry, and if it has one, then $|G| = 4 k$ with $k$ odd.
]
#proofbox[
We use that in a finite $p$-group every proper subgroup is properly contained in its normalizer.

_Step 1._ Let $P$ be a 2-subgroup of $G$. (a) If $sigma in P$ is an involution $!= plus.minus 1$, $-1 in P$ and $H = {plus.minus 1, plus.minus sigma} != P$, then the entry of $sigma$ is non-trivial. Indeed, take $y in N_P (H) without H$; then $y sigma y^(-1)$ is an involution in $H$ other than $plus.minus 1$, so it is $sigma$ or $-sigma$. In the first case $y in C(sigma) without H$ and $|C(sigma)| > 4$; in the second $sigma tilde -sigma$. (b) If $g in P$ satisfies $g^2 = -1$ and $chevron.l g chevron.r != P$, then $N_G (chevron.l g chevron.r) supset.eq N_P (chevron.l g chevron.r) supset.neq chevron.l g chevron.r$, so the entry of $g$ is non-trivial by Lemma 6.1.

_Proof of (ii)._ Apply Step 1 with $P = G$, which has order at least 8 and contains $-1$. Involution entries and Hermitian entries of order 4 are non-trivial, odd entries do not occur, and the remaining entries are non-trivial by definition.

_Proof of (iii)._ If $sigma$ gives a trivial involution entry, let $P$ be a Sylow 2-subgroup of $G$ containing $H = {plus.minus 1, plus.minus sigma}$; by Step 1(a), $P = H$. If $g$ gives a trivial Hermitian entry, let $P$ be a Sylow 2-subgroup containing $chevron.l g chevron.r$; by Step 1(b), $P = chevron.l g chevron.r$. In both cases $|P| = 4$, so $|G| = 4 k$ with $k$ odd. If $L$ had two trivial entries, the corresponding Sylow 2-subgroups $P_1, P_2$ would be conjugate, $u P_1 u^(-1) = P_2$. They cannot be $C_2 times C_2$ and $C_4$. If $P_i = {plus.minus 1, plus.minus sigma_i}$, then $u sigma_1 u^(-1)$ is an involution of $P_2$ other than $plus.minus 1$, hence $plus.minus sigma_2$, and the two entries coincide. If $P_i = chevron.l g_i chevron.r$, then $u g_1 u^(-1) in {plus.minus g_2}$, and again the entries coincide.

_Proof of (i)._ If $|G| = 4$, then $G = {plus.minus 1, plus.minus sigma} tilde.equiv C_2 times C_2$ or $G = chevron.l g chevron.r tilde.equiv C_4$ with $g^2 = -1$. In the first case the only entry is the involution entry of $sigma$; $C(sigma) = G$ and $sigma tilde.not -sigma$ since $G$ is abelian, so it is trivial. In the second case the only entry is the Hermitian entry of $g$, and $N_G (chevron.l g chevron.r) = G = chevron.l g chevron.r$, so it is trivial by Lemma 6.1. Conversely, suppose all entries of $L$ are trivial. By the detection argument of §1.3 and $L in S$, $L$ has an entry; by (iii) it has exactly one, and $|G| = 4 k$ with $k$ odd. If $k > 1$, $G$ contains an element of odd prime order by Cauchy's theorem, whose odd entry is non-trivial. Hence $|G| = 4$. #qed
]

The lemma is a statement about finite groups containing a central involution $-1$; #calc("C21") checks it on all such groups of order at most 190.

== 6.3 The lattices of $S_"gen"$

#rem[
*Terms of this section.* They are defined in the text below; here they are collected in the order in which the argument uses them, and Figure 3 shows how the numbers arise.
- A _record_ is the result of the full-group computation for one constructed lattice: $|"Aut" L|$, the conjugacy classes with their centralizers, and one _row_ for each entry of $L$, with its _identifier_ when one was computed. A record is _successful_ if the computation finished.
- The identifier of a row is _exact_ if it is the class of the entry in a list of §4. A list type (a pair of components of an involution type, or a type $p - (z, f) - s$) is _exact_ if every row of that type in every successful record has an exact identifier.
- For an involution list with $F tilde.equiv K$, the _inverse identifier_ of a class $[phi]$ is the identifier of $[phi^(-1)]$.
- A record is _unlabelled_ if none of its rows received an exact identifier.
- Records are joined into _identification classes_ when they share an exact identifier, when an inverse identifier links them, or when an explicit integer isometry links them.
- The _minimum-shell profile_ and the _full-group signature_ are isometry invariants; they separate the unlabelled lattices.
]

#figure(kind: image, numbering: none, gap: 0.8em,
  _x-diagram(440, 196, h: 42, (
    rec: (0, 0, 100, "plain", [$"17,736"$ successful \ full-group records]),
    ids: (113, 0, 100, "plain", [$"13,396"$ classes \ shared exact identifiers]),
    inv: (226, 0, 100, "plain", [$"12,585"$ classes \ + $"8,716"$ inverse identifiers]),
    iso: (339, 0, 101, "plain", [$"12,546"$ classes \ + 39 isometries among \ the 66 unlabelled records]),
    exh: (0, 76, 136, "plain", [*none is missing* \ Lemma 6.3(ii): the rows carry \ the whole weight of every list]),
    ex: (150, 76, 140, "plain", [$"12,519"$ classes with an \ identifier of an exact type: \ distinct by Lemma 6.3(i)]),
    un: (304, 76, 136, "plain", [27 unlabelled lattices: \ distinct by their \ minimum-shell profiles]),
    res: (110, 152, 260, "result", [$S_"gen"$: $"12,546"$ pairwise non-isometric lattices, \ and every lattice of $S_"gen"$ is one of them]),
  ), (
    ("rec", "r", "ids", "l"), ("ids", "r", "inv", "l"), ("inv", "r", "iso", "l"),
    ("iso", ("b", 0.3), "ex", ("t", 0.7)), ("iso", ("b", 0.6), "un", ("t", 0.6)), ("iso", ("b", 0.1), "exh", ("t", 0.8)),
    ("exh", "b", "res", ("t", 0.15)), ("ex", "b", "res", ("t", 0.5)), ("un", "b", "res", ("t", 0.85)),
  )),
  caption: [Figure 3. How the records of §6.3 become the lattices of $S_"gen"$. Top row: records are joined into identification classes by positive identifications only. Middle row: why different classes are different lattices (Lemma 6.3(i) and the profiles), and why no lattice is missing (Lemma 6.3(ii)).],
)
Every lattice of $S_"gen"$ has a non-trivial involution entry or an odd entry, and every such entry is a double coset in the lists of §4 with $|"Stab"(phi)| > 4$, or with $|"Stab"(phi)| = 4$ and $[phi] = [phi^(-1)]$, or of odd prime order. From entries of these lists (and, starting from their published descriptions (codes, Hermitian and quaternionic structures), from 26 of the 27 lattices of Appendix A; the remaining one, $[4.L_3 (4).2^2]_(32, 2)$ of [Neb98q] with $|"Aut" L| = 2^10 3^2 5 dot 7$, is reached from an entry of type $(12, 20)$) we construct lattices $L_phi$, compute $"Aut" L$ with bliss and its conjugacy classes and centralizers with GAP, and list all entries of $L$ with their weights #calc("C1"). We call the result of this computation for one constructed lattice, that is $|"Aut" L|$, the conjugacy classes and centralizers and the rows of the entries of $L$ with their identifiers, a _record_, and the record _successful_ if the computation finished and returned these data. Two constructed lattices are identified when they share an entry (the same double coset in the same list, including the partner class $[phi^(-1)]$ for $F tilde.equiv K$); the 66 _unlabelled_ constructed records, none of whose computed entries received an exact identifier (their odd entries are of order 5, or of order 3 of the inexact types below) were grouped into 27 lattices (not related to the 27 of Appendix A) by explicit isometries, checked in integers. This gives $"12,546"$ lattices. For each of them:
- as an auxiliary check, the weights computed from all group elements add up to $1 - 2 slash |"Aut" L|$;
- every computed entry that can be matched with an entry of the lists of §4 has the same $|C|$ and the same auxiliary assignment count $e$ computed inside the full group [Aux] (no mismatch).
For the ten lattices with $|"Aut" L| > 10^7$, all identified with known lattices or with each other by explicit isometries, the conjugacy classes of two of them ($"BW"_32$ and one code lattice, with $|"Aut" L| = 2^31 3^5 5^2 7 dot 17 dot 31$ and $2^31 3^5 5^2 7$) were obtained from the structure of the group, not in a faithful permutation representation (the second paper of this repository computes them directly, on the minimal vectors, and finds the same classes of involutions and of elements of odd prime order, with the same centralizer orders): with $E = O_2 ("Aut" L)$, the classes of $"Aut" L$ are obtained from the classes of $"Aut" L slash E$ and, for a representative $t E$ of each, the orbits on the coset $t E$, under conjugation, of the preimage of the centralizer of $t E$ in $"Aut" L slash E$ #calc("C1").

For these two lattices the full group was also computed from their Gram matrices alone, with no supplied group order or generating set #calc("C22"). The input Grams are exactly those used for the known representatives in this classification. The code lattice has order $"91,321,742,131,200"$ and $"BW"_32$ has order $"48,126,558,103,142,400"$. Each is computed by the algorithm of §7.1 as twice the orbit size of a minimal line $plus.minus x$ times the order of the stabilizer of $x$; every candidate image of the line is tested, so this is the order of the full group. The groups in which the conjugacy classes were computed #calc("C1") have these orders, so they are the full automorphism groups. These full-group computations are separate from the structural conjugacy-class calculations above. The rerun uses the same algorithm as the earlier full-group calculations, not an independent mathematical method.

The positive identifications among the unlabelled records are preserved as 39 unimodular integer isometries #calc("C23"). Each satisfies $T^sans(T) G_a T=G_b$ and $|det T|=1$. The 27 resulting unlabelled lattices are pairwise non-isometric by the following invariant #calc("C24"). Let $overline(X) subset L slash 2 L$ be the image of the lines $plus.minus x$ of minimal vectors; the $"73,440"$ images are distinct, which the program checks. For $v != 0$ let $c(v)$ be the number of unordered pairs ${a, b} subset overline(X)$ of distinct elements with $a + b = v$, and attach to each $a in overline(X)$ the integers $s_j (a) = sum_(b in overline(X), b != a) c(a + b)^j$, $j = 1, 2, 3$. An isometry of lattices induces a linear isomorphism of $L slash 2 L$ that maps $overline(X)$ onto $overline(X)$, so the multiset of the triples $(s_1, s_2, s_3)(a)$ is an isometry invariant, and so is the triple chosen from it by least frequency and then lexicographic order. We call this triple the _minimum-shell profile_. The profiles of the 27 lattices are pairwise different. Different profiles prove non-isometry; equal profiles are never used to infer an isometry.

As a second, independent reason, the unlabelled records are separated from all others by their _full-group signature_: the order of $"Aut" L$ together with the multiset, over the entries of $L$, of their kind, their type and the order of their centralizer, and for Hermitian entries also $|"U"|$ and $f$. These quantities are determined by $"Aut" L$ and its conjugacy classes, so they are isometry invariants; no record with an exact identifier has the signature of an unlabelled record #calc("C24"). Of the ten lattices with $|"Aut" L| > 10^7$, nine are lattices of Appendix A, and 54 explicit isometries identify their records with them; the tenth, with $|"Aut" L| = 2^17 3^4 5$, is not among the lattices of Appendix A, and 13 explicit isometries identify its records with each other; all 67 are checked in integers #calc("C24").

The identifications above are positive: two records are put together when they share an exact identifier (defined below), when their entries are the two classes $[phi]$, $[phi^(-1)]$ of one involution entry with $F tilde.equiv K$ (§4.1.2; all $"8,716"$ such inverse identifiers ($"5,005"$ unordered pairs, $"3,711"$ of them recorded in both directions; the certificate also lists the reverses of six of them) were recomputed #calc("C24")), or by an explicit integer isometry #calc("C23", "C24"). Call the resulting classes of records _identification classes_. Step by step, the $"17,736"$ successful full-group records fall into $"13,396"$ classes by shared exact identifiers, into $"12,585"$ when the inverse identifiers are added, and into $"12,546"$ by the 39 isometries among the 66 unlabelled records; the 67 isometries for the ten lattices with $|"Aut" L| > 10^7$ join no further classes, they confirm identifications already made by identifiers #calc("C24", "C25"). Distinctness of different identification classes needs a separate argument, because not every computed entry receives an identifier.

The identifier of an entry row is _exact_ if it is the class of the entry in a list of §4. For involutions it is computed from the identified components and the canonical form of the gluing (§7.2); for four rows of type $(8, 24)$ in groups of order above $10^7$ it was completed from the canonical key of type $(8, 24)$ (§7.2), whose saved prefix matches a unique lattice $K$, and from the centralizer order. For odd prime order the components are identified and an identifier is attached only if exactly one extremal class of the list has those components and that $|"Stab"(phi)|$; otherwise the row has none. A single candidate is the class of the entry by the completeness of the list, the identification of the components and the correctness of $|C(h)| = |"Stab"(phi)|$. For each pair of components there is one ordered list, so the identifier of an involution row does not depend on whether $sigma$ or $-sigma$ represents the entry. A list type (a pair of components of an involution type, or a type $p - (z, f) - s$) is _exact_ if every row of that type in every successful full-group record has an exact identifier.

#thm[
*Lemma 6.3 (distinctness and exhaustion).* Assume that each successful record lists all conjugacy classes of the full automorphism group of its lattice, with their types and centralizer orders, and that its identifiers are exact as defined.
#set enum(numbering: "(i)")
+ If some record of an identification class $cal(I)_1$ has an identifier $x$ on a row of an exact type $T$, then the lattice of $cal(I)_1$ is not isometric to the lattice of any other identification class.
+ Let $L_1, dots, L_t$ be pairwise non-isometric lattices, one for each identification class. For a list $ell$ of §4.1–§4.2 (for $p >= 5$, the lists of the types with the same $p$ and $f$ taken together) give an extremal class $X$ the weight $w(X) = 1 slash |"Stab"(phi)|$, and give a row of type $ell$ of some $L_i$, with element $h$, the weight $1 slash |C(h)|$, or $2 slash |C(h)|$ if $ell$ is an involution list with $F tilde.equiv K$ and $h tilde.not -h$. Let $T_0 (ell)$ be a set of trivial classes of $ell$ (§6.1; empty for the lists of odd prime order, which have none) such that no class of $T_0 (ell)$, and for $F tilde.equiv K$ no partner $[phi^(-1)]$ of one, is a class of an entry of any $L_i$. If for every $ell$ the rows of type $ell$ of $L_1, dots, L_t$ have total weight $sum w(X)$ over the extremal classes $X$ of $ell$ not in $T_0 (ell)$, then every lattice of $S_"gen"$ is isometric to one of $L_1, dots, L_t$.
]
_Strategy of the proof._ Part (i) transports an entry: an isometry of lattices carries an entry with an exact identifier to an entry with the same identifier, and records that share an identifier were put into one identification class. Part (ii) reads identity (4.1) backwards. Each row of the records is an entry: one class of a list of §4, or, for an involution list with $F tilde.equiv K$, the classes $[phi]$ and $[phi^(-1)]$, which may be two. Distinct rows give disjoint sets of classes. If the rows already carry the whole weight of the list (less the trivial classes $T_0 (ell)$), no class is left over, so no lattice of $S_"gen"$ can be missing.
#proofbox[
_(i)._ Let $cal(I)_2$ be an identification class whose lattice is isometric to that of $cal(I)_1$, and let $u$ be an isometry. It maps the entry of $x$ to an entry of the lattice of $cal(I)_2$ of the same type $T$ and the same class in the list of §4, since isometric pairs have the same components and the same double coset (Proposition 4.1). A record of $cal(I)_2$ lists all entries of its full group, so it has a row for this entry; rows of odd prime order are indexed by conjugacy classes, so the transported class has its own row. As $T$ is exact, that row carries the identifier of its class, which is $x$, or for an involution type with $F tilde.equiv K$ possibly $x$ or the inverse class $x'$. Records sharing $x$ lie in one identification class, and $x$ and $x'$ are joined by an inverse identifier. Hence $cal(I)_2 = cal(I)_1$.

_(ii)._ A row of type $ell$ of $L_i$ with element $h$ gives the pair $(L_i, h)$, which is extremal; by Proposition 4.1 and the completeness of $ell$ (§4), its isometry class is an extremal class $X$ of $ell$ with $|"Stab"(phi)| = |C(h)|$. For an involution list with $F tilde.equiv K$ the row stands for the entry ${[h], [-h]}$, whose classes $[phi]$, $[phi^(-1)]$ in $ell$ coincide exactly when $h tilde -h$ (§6.1); in every case the weight of the row is the total weight of the classes of its entry. Distinct rows give disjoint sets of classes: two rows of one $L_i$ are distinct entries, and $(L_i, h) tilde.equiv (L_i, h')$ only if $h'$ is conjugate to $h$ in $"Aut" L_i$; rows of different $L_i$ have non-isometric lattices. Hence the rows of type $ell$ have the total weight of the set $A(ell)$ of classes of $ell$ that are classes of entries of $L_1, dots, L_t$, and $A(ell)$ is disjoint from $T_0 (ell)$. By hypothesis $A(ell)$ and the set of extremal classes not in $T_0 (ell)$ have the same total weight; as all weights are positive, they are equal. Now let $L in S_"gen"$. It has an odd entry or a non-trivial involution entry; a class of this entry lies in one of the lists of §4.1–§4.2, which cover all types that occur (§3.3), and is extremal, and it is not trivial, so it is not in $T_0 (ell)$. Hence it is the class of an entry of some $L_i$: for elements $h$ of $"Aut" L$ and $h_i$ of $"Aut" L_i$, $(L, h) tilde.equiv (L_i, h_i)$, so $L tilde.equiv L_i$. #qed
]

The exact types are all lists of involutions (including type $(8, 24)$ after the completion just described) and the types $3 - (12, 8) - 8$, $3 - (10, 12) - 6$, $3 - (10, 12) - 8$, $3 - (10, 12) - 10$ and $3 - (9, 14) - 7$. The types $3 - (9, 14) - 9$ and $3 - (8, 16) - 8$ have 401 and 50 rows with two to twelve candidates; the types $3 - (13, 6) - 5$, $3 - (6, 20) - 6$, $3 - (16, 0) - 0$ and all types with $p >= 5$ have no exact identifiers. Of the $"12,546"$ identification classes, $"12,519"$ contain a record with an identifier on a row of an exact type #calc("C25")\; by Lemma 6.3 each of them is distinct from every other class. The remaining 27 classes are exactly the unlabelled lattices above, which are pairwise distinct by their minimum-shell profiles. Hence the $"12,546"$ lattices of $S_"gen"$ are pairwise non-isometric. By Lemma 6.3(ii) they are all of $S_"gen"$: for every list, the rows of the $"12,546"$ lattices, one record for each, have the total weight of the extremal classes of the list less the trivial classes that §6.6 assigns to $S_(C_2^2)$, and none of these is the identifier of a row #calc("C25"). This holds for each of the 83 involution lists, for type $(8, 24)$ and for the types of odd prime order, and is checked from the saved rows (§7.4). It covers also the types whose rows carry no identifier, and the entries that were never processed because their lattices had already been reached. The full-group signatures give a second reason why the unlabelled lattices differ from the others. The ambiguous rows cannot simply be ignored: for 297 of the 451 ambiguous order-3 rows some candidate lies in another identification class #calc("C25").

== 6.4 The lattices of $S_(C_2^2)$
#thm[
*Proposition 6.4.* A lattice $L in S$ lies in $S_(C_2^2)$ if and only if it has a trivial involution entry and is not in $S_"gen"$. Then $"Aut" L = {plus.minus 1, plus.minus sigma}$, and $L$ has exactly one entry.
]
#proofbox[
If $L in S_(C_2^2)$, it has an involution $sigma != plus.minus 1$, and since $L in.not S_"gen"$ its entry is trivial. Conversely a trivial involution entry gives an involution $!= plus.minus 1$. By Lemma 6.2(iii), $|"Aut" L| = 4 k$ with $k$ odd; since $L in.not S_"gen"$, $"Aut" L$ has no element of odd prime order, so $k = 1$ by Cauchy's theorem, and Lemma 6.2(i) applies. #qed
]

== 6.5 The lattices of $S_"CQ"$
A finite 2-group with exactly one involution is cyclic or generalized quaternion $Q_(2^n) = chevron.l x, y : x^(2^(n-1)) = 1, y^2 = x^(2^(n-2)), y x y^(-1) = x^(-1) chevron.r$, $n >= 3$. Indeed, such a group has no noncyclic abelian subgroup, since a noncyclic abelian 2-group contains $C_2 times C_2$ and hence three involutions; and a $p$-group without noncyclic abelian subgroups is cyclic or, for $p = 2$, generalized quaternion of order $2^n$, $n >= 3$ [Gor80, Thm 5.4.10 (ii)]. For $L in S_"CQ"$ every element $g != plus.minus 1$ has 2-power order $2^a$ with $g^(2^(a-1)) = -1$, so all entries of $L$ are Hermitian entries.

#thm[
*Proposition 6.5.* Let $L in S$ have a Hermitian entry of order 4, $g^2 = -1$, which is not an entry of a lattice of $S_"gen"$.
+ If $g tilde.not -g$, then $L in S_"CQ"$, $"Aut" L$ is cyclic, and $|"Aut" L| = |"U"|$. This entry is the only Hermitian entry of order 4 of $L$.
+ If $g tilde -g$, then $L in S_"CQ"$, $"Aut" L tilde.equiv Q_(2^n)$, and $L$ has exactly three Hermitian entries of order 4, all with $g tilde -g$, with $|"U"| = 4, 4, 2^(n-1)$.
]
#proofbox[
$L in.not S_"gen"$ by assumption, and $L in.not S_(C_2^2)$ since $C_2 times C_2$ has no element with $g^2 = -1$. So $L in S_"CQ"$ and $G = "Aut" L$ is cyclic or generalized quaternion. If $G = C_(2^n)$, its elements of order 4 are $plus.minus g$; $G$ is abelian, so $g tilde.not -g$ and $"U" = G$. If $G = Q_(2^n)$, its elements of order 4 are $plus.minus x^(2^(n-3))$, conjugate under $y$, with centralizer $chevron.l x chevron.r$ of order $2^(n-1)$, and the elements $y x^i$, which form two conjugacy classes ($i$ even, $i$ odd), each closed under $y x^i |-> -y x^i = y x^(i + 2^(n-2))$, with centralizer $chevron.l y x^i chevron.r$ of order 4. In particular every element of order 4 of $Q_(2^n)$ is conjugate to its negative, so $g tilde.not -g$ forces $G$ to be cyclic. In both cases the claims follow. #qed
]

== 6.6 Subtraction of $S_"gen"$ and the parts $S_(C_2^2)$, $S_"CQ"$

#thm[
*Lemma 6.6 (subtraction).* Let $cal(T)$ be the set formed from the extremal classes $[phi]$ of the lists of §4.1 with $|"Stab"(phi)|=4$ and $[phi]!=[phi^(-1)]$, taking the class $[phi]$ if $F tilde.not K$ (where $[phi]!=[phi^(-1)]$ holds automatically, the inverse lying in the list with the components interchanged) and the pair ${[phi],[phi^(-1)]}$ if $F tilde.equiv K$; by Proposition 4.1(iii) and §4.1.2 these are exactly the trivial involution entries, and let $cal(G)$ be the list of Proposition 4.6. Assume:
#set enum(numbering: "(1)")
+ the lists of §4.1–§4.2 and $cal(G)$ are complete, the classes of $cal(G)$ are pairwise distinct up to isometry and complex conjugation, and for each class $X in cal(G)$ the order $|"U"(X)|$ and the value $f(X)$ are correct;
+ $L_1, dots, L_t$ are pairwise non-isometric lattices of $S_"gen"$ and every lattice of $S_"gen"$ is isometric to one of them (§6.3, Lemma 6.3(i) and (ii));
+ for each $i$, the trivial involution entry of $L_i$, if any, is identified with its element of $cal(T)$, and every Hermitian entry of order 4 of $L_i$ is identified with its class in $cal(G)$; here the entries are those of the full group $"Aut" L_i$.
Let $cal(T)_"gen" subset cal(T)$ and $cal(G)_"gen" subset cal(G)$ be the sets so identified. Then:
#set enum(numbering: "(i)")
+ $L |-> $ (its unique entry) is a bijection from $S_(C_2^2)$ onto $cal(T) without cal(T)_"gen"$, and $|"Aut" L|=4$ for these lattices;
+ $L |-> $ (its Hermitian entry of order 4) is a bijection from the lattices of $S_"CQ"$ with cyclic group onto ${X in cal(G) without cal(G)_"gen" : f(X)=2}$, and $|"Aut" L|=|"U"(X)|$;
+ the relation "same underlying lattice" partitions ${X in cal(G) without cal(G)_"gen" : f(X)=1}$ into sets of three; these sets correspond bijectively to the lattices of $S_"CQ"$ with quaternion group, and for a set with largest $|"U"|=2^(n-1)$ the group is $Q_(2^n)$.
In particular $S$ is the disjoint union of ${L_1,dots,L_t}$ and the lattices described in (i)–(iii).
]
_Strategy of the proof._ The lists of §4 contain every trivial involution entry and every Hermitian entry of order 4, of every lattice. The entries that belong to lattices of $S_"gen"$ are known from the records of §6.3. The remaining entries belong to lattices outside $S_"gen"$, and the finite-group results say how many entries each of these lattices has: one for $C_2 times C_2$ (Proposition 6.4) and for a cyclic group, three for a quaternion group (Proposition 6.5). So counting the remaining entries counts these lattices.
#proofbox[
_(i)._ Let $tau in cal(T)$ with underlying lattice $L$. If $L in S_"gen"$, then $L tilde.equiv L_i$ for some $i$ by (2), and $tau$ is a trivial involution entry of $L_i$. By Lemma 6.2(iii), $L_i$ has at most one, so $tau in cal(T)_"gen"$ by (3). Conversely, an element of $cal(T)_"gen"$ has underlying lattice $L_i in S_"gen"$. Hence $cal(T) without cal(T)_"gen"$ consists of the trivial involution entries of lattices outside $S_"gen"$, and these lattices lie in $S_(C_2^2)$ by Proposition 6.4. Every $L in S_(C_2^2)$ has a trivial involution entry by Proposition 6.4, and it occurs in $cal(T)$ by the completeness in (1). Distinct elements of $cal(T)$ are distinct entries by Proposition 4.1(ii) and §4.1.2. Since $L$ has exactly one entry by Proposition 6.4, two distinct elements of $cal(T) without cal(T)_"gen"$ have non-isometric lattices.

_(ii) and (iii)._ Let $X in cal(G)$ be the class of $(L, plus.minus g)$. As in (i), $X in cal(G)_"gen"$ if and only if $L in S_"gen"$; here the identification in (3) gives $X$ itself because the classes of $cal(G)$ are distinct by (1). If $L in.not S_"gen"$, then $L in S_"CQ"$ by Proposition 6.5. By §2.4, $f(X)=1$ exactly when $g$ is conjugate to $-g$ in $"Aut" L$. Every $L in S_"CQ"$ occurs: $"Aut" L$ is a 2-group of order greater than 2 (as $L in S$) whose only involution is $-1$, so it contains $g$ with $g^2=-1$, whose class lies in $cal(G)$ by (1). If $f(X)=2$, Proposition 6.5(i) shows that $X$ is the only Hermitian entry of order 4 of $L$ and that $|"Aut" L|=|"U"(X)|$; this gives (ii). If $f(X)=1$, Proposition 6.5(ii) shows that $L$ has exactly three such entries, all with $f=1$, and that their unitary groups have orders $4,4,2^(n-1)$ with $"Aut" L tilde.equiv Q_(2^n)$. Distinct lattices have disjoint sets of entries, which gives (iii). The last statement follows from the definitions of the three parts and Proposition 3.1. #qed
]

The hypotheses (2) and (3) are where the computations of §6.3 enter. For Hermitian entries, (3) is met through the invariance of the keys (Lemma 7.2) and the common key function of §7.2; for an involution entry with $F tilde.equiv K$, through the identifier of one of $[phi]$, $[phi^(-1)]$ and the recomputed table of $"8,716"$ inverse identifiers #calc("C28", "C24"). Hypothesis (3) requires the entries of the _full_ group: an entry missing from the computed list of some $L_i$ would leave its class in $cal(T) without cal(T)_"gen"$ or $cal(G) without cal(G)_"gen"$, and $L_i$ would be counted twice; a wrong identification would remove the class of a different lattice. By Lemma 6.2(iii), only representatives with $|"Aut" L_i| equiv 4 med (mod 8)$ can have a trivial involution entry. No odd entry and no identifier of an odd entry is used here. In the saved full-group records, every involution type is exact in the sense of §6.3 #calc("C25"), so every trivial involution entry has an exact pair identifier; no record has two trivial involution entries, and no record with a trivial involution entry has a group order not congruent to 4 modulo 8; every Hermitian entry of order 4 has a key, and all $"11,545"$ distinct keys of these $"16,163"$ entries belong to the list of Proposition 4.6, with the same $|"U"|$ and $f$ #calc("C29"). The four involution entries of type $(8, 24)$ whose identifiers were completed from the canonical key of type $(8, 24)$ (§§6.3, 7.2) belong to lattices with $|"Aut" L| > 10^7$; they are non-trivial and do not enter $cal(T)_"gen"$ in any case. These record checks confirm the bookkeeping of (3); that each record lists all conjugacy classes of the full group is the software assumption of Table 10 on bliss and GAP.

*Counting the two parts.* A lattice of $S_"gen"$ can also have a trivial involution entry; then $|"Aut" L| = 4 k$ with $k > 1$ odd by Lemma 6.2(iii), so it has an odd entry, and the computation of §6.3 lists that trivial entry among its entries. We say that the double cosets so listed are _met by $S_"gen"$_. Hence, by Lemma 6.6(i), $S_(C_2^2)$ is in bijection with the double cosets not met by $S_"gen"$ with $|"Stab"(phi)| = 4$ and $[phi] != [phi^(-1)]$, except that for $F tilde.equiv K$ the two classes $[phi]$, $[phi^(-1)]$ of one entry give one lattice. Counting them gives $"67,284"$ #calc("C26", "C28")\; for the pair with $s = 16$, $delta = 1$, for instance, $"33,006"$ classes not met by $S_"gen"$ give $"16,503"$ lattices.

The Hermitian entries of order 4 of all lattices in $cal(L)$ are the $"76,742"$ classes of Proposition 4.6. Removing those that are entries of lattices of $S_"gen"$ leaves $"65,197"$ classes, which Lemma 6.6(ii)–(iii) describes. By Proposition 6.5(i), the $"62,650"$ classes with $g tilde.not -g$ are the lattices of $S_"CQ"$ with cyclic group, with $|"Aut" L| = |"U"| = 4, 8, 16$ for $"62,468"$, $175$, $7$ of them. The values of $|"U"| = |C(g)|$ are those of the enumeration of Proposition 4.6, where for the $"62,468"$ classes with $|"U"| = 4$ the canonical form is discrete, which forces $"U" = chevron.l g chevron.r$, and for the $182$ classes with $|"U"| >= 8$ they come from a backtrack search and were all confirmed with bliss. The remaining $"2,547"$ classes, with $g tilde -g$, fall into $849$ triples with equal isometry invariants, and each triple was shown to be a single lattice by two explicit isometries checked in integers #calc("C27", "C28"). By Proposition 6.5(ii), the patterns of $|"U"|$ in the triples, $(4, 4, 4)$, $(4, 4, 8)$, $(4, 4, 16)$, $(4, 4, 32)$, occurring $793$, $48$, $6$, $2$ times, give the quaternion groups of orders $8$, $16$, $32$, $64$. The saved triples form an exact partition of these 2,547 Gaussian classes. Different triples cannot represent the same integer lattice, since Proposition 6.5 gives exactly three Gaussian entries for each quaternion lattice and Proposition 4.6 makes the entries distinct. This proves Table 8.

#ptable(
  columns: 5, align: (left, right, right, right, left),
  caption: [Table 8. The lattices of $S_"CQ"$ (Theorem B).],
  table.header([$"Aut" L$], [lattices], [$1 - 2 slash |"Aut" L|$], [total], [entries of order 4]),
  [$C_4$], [$"62,468"$], [$1 slash 2$], [$"31,234"$], [1],
  [$C_8$], [$175$], [$3 slash 4$], [$525 slash 4$], [1],
  [$C_16$], [$7$], [$7 slash 8$], [$49 slash 8$], [1],
  [$Q_8$], [$793$], [$3 slash 4$], [$"2,379" slash 4$], [3],
  [$Q_16$], [$48$], [$7 slash 8$], [$42$], [3],
  [$Q_32$], [$6$], [$15 slash 16$], [$45 slash 8$], [3],
  [$Q_64$], [$2$], [$31 slash 32$], [$31 slash 16$], [3],
  [*total*], [*$"63,499"$*], [], [*$"512,251" slash 16$*], [],
)

#rem[
*Example (one lattice of each kind).* This example is not used in the proofs.
- `32.4.112978` lies in $S_(C_2^2)$. Its only entry is an involution entry for the pair of components $sqrt(2) (D_8^2)^+$, $sqrt(2) (D_8^2)^+$ ($s = 16$, $delta = 1$, the last row of Table 6). The entry consists of two classes $[phi] != [phi^(-1)]$ of that list, both with $|"Stab"(phi)| = 4$, and neither is met by $S_"gen"$. By Lemma 6.6(i) the two classes give one lattice, with $"Aut" L = {plus.minus 1, plus.minus sigma}$; they are the pairs $(L, sigma)$ and $(L, -sigma)$. The lattice contributes $1 - 2 slash 4 = 1 slash 2$ to (5.1).
- `32.4.1` lies in $S_"CQ"$, with $"Aut" L tilde.equiv C_4$. Its only entry is one Gaussian class of Proposition 4.6, with $f = 2$ and $|"U"| = 4$. By Lemma 6.6(ii), $"Aut" L = chevron.l g chevron.r$.
- `32.8.4736` lies in $S_"CQ"$, with $"Aut" L tilde.equiv Q_8$. It has three Gaussian classes, each with $f = 1$ and $|"U"| = 4$: the pattern $(4, 4, 4)$ of §6.6. By Lemma 6.6(iii) they are one lattice with group $Q_8$, and the three classes are its three subgroups of order 4.
]

== 6.7 The total from the classified lattices

The three parts proved complete in §§6.3–6.5 are disjoint and exhaust $S$. Consequently the correction in (5.1) is obtained by adding $1-2 slash |"Aut" L|$ over their representatives. Each lattice of $S_(C_2^2)$ contributes $1 slash 2$ (its group has order 4); the cyclic and quaternionic contributions use Table 8, and the general part uses the full-group orders of §6.3.

The saved database aggregation #calc("C13") computes this correction directly from the classified rows, with no expected total as its pass condition. This is the arithmetic used in §5.1. The element-assignment computation of [Aux] shares the component lists, so agreement with it is an auxiliary consistency check, not an independent proof of completeness. Because that computation uses no full automorphism group of rank 32, it is a global check of rows 10 and 11 of Table 7 by a different arithmetic route (§5.2); it does not test the shared component and gluing lists.

For the general-part calculation, every classified general-part row is linked to an original successful full-group record with the same group order #calc("C13"). That the entries which were not processed, because their lattices had already been reached from other entries, are entries of these lattices is Lemma 6.3(ii), for all types, including those without identifiers. These provenance checks do not establish that the original lists of pairs are exhaustive. The completeness arguments of §4 and Lemma 6.3 remain necessary.

== 6.8 The database and the distribution of $|"Aut" L|$
*Numbering.* In [D] each $L in S$ has the label `32.`$|"Aut" L|$`.`$i$, where $i = 1, 2, dots$ numbers the lattices with the same $|"Aut" L|$ in the following order.
+ The isomorphism type of $"Aut" L$: groups with an identifier $[n, k]$ in the Small Groups Library by $k$, then groups identified by a structure description by that description, then the remaining $19$ groups (nine identified by the orders of a chief series and the ten of order above $10^7$) by the orders of a chief series; these $19$ lattices have $18$ distinct orders, and the two of order $"196,608"$ have chief series of different orders (recorded in `data/classified-lattices.jsonl`).
+ The part of Theorem B: $S_"gen"$, then $S_(C_2^2)$, then $S_"CQ"$. The part is recorded in the field `part` of `data/classified-lattices.jsonl`; the order of the parts is the one just given.
+ A key, which is an isometry invariant of $L$ (the first two kinds below are computed from the lists of §4): the smallest identifier of a class of the lists of §4.1 and §4.2 that is an entry of $L$ (all lattices of $S_(C_2^2)$ and $"12,519"$ of $S_"gen"$); otherwise the smallest canonical form or invariant $D$ (§7.2) of a pair $(L, plus.minus g)$ with $g^2 = -1$ (all lattices of $S_"CQ"$ and $5$ of $S_"gen"$); otherwise the hash of the configuration modulo 2 of the lines spanned by vectors of norm 4 ($22$ lattices of $S_"gen"$). Within each part and group these keys are distinct.
This order is a convention without mathematical meaning. The labels were fixed with version 1 of [D] and are kept in later versions. The Gram matrices and generators of [D] are those of the computations of §§6.3–6.5; in version 1 the generators of the 182 lattices with $"Aut" L = C_8$ or $C_16$ generated only the subgroup of order 4, and version 2 adds a generator of $"Aut" L$, found by the search of Proposition 4.6. The last label is that of $"BW"_32$, `32.48126558103142400.1`.


Table 9 gives the distribution of $|"Aut" L|$ over $cal(L)$. For the $"12,546"$ lattices of $S_"gen"$ it lists the orders computed with bliss (§6.3, Table 10); the $"67,284"$ lattices of $S_(C_2^2)$ have $|"Aut" L| = 4$, those of $S_"CQ"$ are in Table 8, and the remaining $"10,900,939"$ lattices have $"Aut" L = {plus.minus 1}$ (Corollary C). The column of numbers sums to $|S| = "143,329"$ without the first row and to $N$ with it, and the orders satisfy (5.1) as an equality of fractions #calc("C31", "C30"). In $S_"gen"$, $104$ distinct orders occur. The most frequent are $8$ ($"7,055"$ lattices), $6$ ($"2,308"$) and $16$ ($"1,233"$); the largest is $|"Aut" "BW"_32|$.

The group of each lattice of $S_"gen"$ was also identified #calc("C30"). The generators are checked in integers to be automorphisms, and the conjugacy classes, with their sizes summing to $|"Aut" L|$, were computed with GAP. For $|"Aut" L| <= 10^7$, that is, for all but ten lattices of $S_"gen"$, the group was identified as follows. For $"12,447"$ lattices it was identified by its number in the Small Groups Library of GAP, with $365$ distinct groups; the most frequent are $C_4 times C_2$ ($"3,737"$ lattices), $C_2^3$ ($"2,495"$), $C_6$ ($"2,308"$) and $D_8$ ($823$). For $80$ lattices it was described by the structure description of GAP (a description, not a proof of isomorphism), and for the remaining nine only by the orders of a chief series. Each of the $"12,536"$ groups was also given by a presentation. Its relations were checked on the generators, and the order of the finitely presented group was computed and equals $|"Aut" L|$. For nine groups of order above $10^5$ the identification run did not count this order; it was counted afterwards, by coset enumeration along a chain of point stabilizers in a permutation representation of small degree, with Reidemeister–Schreier presentations of the stabilizers #calc("C30"). None of these identifications is used in Theorems A and B.

#ptable(
  columns: 8, align: (right, right, right, right, right, right, right, right),
  caption: [Table 9. The number of $L in cal(L)$ with a given $|"Aut" L|$, read down the columns. Only the orders $4, 8, 16, 32, 64$ receive lattices of $S_(C_2^2)$ or $S_"CQ"$ (Table 8); all others except 2 (the lattices outside $S$) are orders of lattices of $S_"gen"$ #calc("C30").],
  table.header([$|"Aut" L|$], [lattices], [$|"Aut" L|$], [lattices], [$|"Aut" L|$], [lattices], [$|"Aut" L|$], [lattices]),
  [$2$], [$"10,900,939"$], [$240$], [$3$], [$"3,072"$], [$5$], [$"62,208"$], [$1$],
  [$4$], [$"129,752"$], [$256$], [$27$], [$"4,032"$], [$1$], [$"65,536"$], [$1$],
  [$6$], [$"2,308"$], [$288$], [$19$], [$"4,096"$], [$1$], [$"196,608"$], [$2$],
  [$8$], [$"8,023"$], [$320$], [$2$], [$"4,608"$], [$2$], [$"207,360"$], [$1$],
  [$10$], [$6$], [$324$], [$1$], [$"4,896"$], [$1$], [$"245,760"$], [$1$],
  [$12$], [$674$], [$384$], [$21$], [$"5,184"$], [$1$], [$"262,144"$], [$1$],
  [$16$], [$"1,288"$], [$432$], [$2$], [$"5,760"$], [$2$], [$"294,912"$], [$1$],
  [$18$], [$4$], [$480$], [$5$], [$"6,144"$], [$2$], [$"322,560"$], [$1$],
  [$20$], [$13$], [$512$], [$12$], [$"7,680"$], [$2$], [$"368,640"$], [$1$],
  [$24$], [$272$], [$576$], [$3$], [$"8,192"$], [$1$], [$"393,216"$], [$1$],
  [$32$], [$290$], [$672$], [$1$], [$"9,216"$], [$3$], [$"589,824"$], [$1$],
  [$36$], [$17$], [$768$], [$10$], [$"10,240"$], [$1$], [$"839,808"$], [$1$],
  [$40$], [$23$], [$800$], [$1$], [$"10,368"$], [$1$], [$"1,179,648"$], [$1$],
  [$48$], [$163$], [$864$], [$1$], [$"11,520"$], [$1$], [$"1,548,288"$], [$1$],
  [$54$], [$2$], [$896$], [$1$], [$"12,288"$], [$1$], [$"2,073,600"$], [$1$],
  [$64$], [$112$], [$960$], [$2$], [$"15,360"$], [$1$], [$"26,542,080"$], [$1$],
  [$72$], [$17$], [$"1,024"$], [$7$], [$"16,128"$], [$1$], [$"44,040,192"$], [$1$],
  [$80$], [$5$], [$"1,152"$], [$5$], [$"16,384"$], [$1$], [$"53,084,160"$], [$1$],
  [$96$], [$39$], [$"1,440"$], [$2$], [$"17,280"$], [$1$], [$"74,317,824"$], [$1$],
  [$108$], [$2$], [$"1,536"$], [$6$], [$"18,432"$], [$1$], [$"975,175,680"$], [$1$],
  [$120$], [$8$], [$"1,728"$], [$2$], [$"19,584"$], [$1$], [$"1,244,160,000"$], [$1$],
  [$128$], [$61$], [$"1,920"$], [$1$], [$"20,736"$], [$1$], [$"1,509,949,440"$], [$1$],
  [$144$], [$11$], [$"2,016"$], [$1$], [$"23,328"$], [$1$], [$"16,911,433,728"$], [$1$],
  [$160$], [$3$], [$"2,048"$], [$5$], [$"24,000"$], [$1$], [$"91,321,742,131,200"$], [$1$],
  [$192$], [$23$], [$"2,304"$], [$2$], [$"27,648"$], [$2$], [$"48,126,558,103,142,400"$], [$1$],
  [$216$], [$4$], [$"2,592"$], [$1$], [$"51,840"$], [$1$], [], [],
  [$224$], [$1$], [$"2,880"$], [$1$], [$"58,752"$], [$1$], [], [],
)

= 7. Computations and assumptions <s7>

#policy-box[
*What this section does.* §7.1 lists the programs. §7.2 states the two lemmas through which our lists use keys (Lemmas 7.1 and 7.2) and proves that each key is an isometry invariant. §7.3 lists everything the proofs assume (Table 10). §7.4 lists the computations with their certificates and running times (Table 11) and describes the files that accompany this paper and how they are checked. The research record [Rec] keeps the programs and runs behind the lists.
]

== 7.1 Programs
All computations used one environment: PARI/GP 2.17.3 [PARI] (through cypari2 2.2.4) for `qfminim`, `qfauto` and `qfisom`; GAP 4.15.1 [GAP] for orders of permutation groups, double cosets, conjugacy classes and centralizers; bliss 0.77 [JK07, JK11], called through its C interface, for canonical labellings and automorphism groups of graphs; SageMath 10.9 [Sage] for genus symbols and masses of genera; and our own code in Python 3.12 and C. The exact linear programs of Proposition 3.3 were solved with the Parma Polyhedra Library (PPL) through SageMath; the proof uses only their finite rational certificates (§B.1.1), which are checked without it. Calls of `qfauto` and `qfisom` go through a wrapper that retries with random changes of basis and checks every returned matrix in integers.

*Full automorphism groups.* The orders $|"Aut" L|$ of the lattices of $S_"gen"$ and of the two exceptional known lattices are computed by one algorithm #calc("C1", "C22"). Let $X$ be the set of vectors of norm 4 of $L$; the program checks that it has the $"146,880"$ elements of an extremal lattice of rank 32 (§B.3). Fix $x in X$. The _rooted local Gram graph_ of $x$ has as vertices $x$ and the vectors $y in X$ with $(x, y) = 2$, and every pair of vertices is joined by an arc coloured by their inner product (no arc is omitted); $x$ receives its own vertex colour, so that the automorphisms of the coloured graph are exactly those fixing $x$. These vectors span $L times.o QQ$, which is checked. An automorphism of $L$ fixing $x$ permutes these vertices and preserves the colours, and it is determined by that permutation. Conversely every colour-preserving automorphism of the graph that fixes $x$ extends to a unique rational isometry, and the stabilizer $"Aut"(L)_x$ consists of those extensions that preserve $L$. They are found exactly: with $M$ the sublattice spanned by the vertices, they are the elements stabilizing the subgroup $L slash M$ of the discriminant group of $M$ (a set stabilizer in GAP; when that group is not 2-elementary, a Schreier orbit of $L$ under the rational group). Next, every line $plus.minus y subset X$ is a candidate image of $plus.minus x$ unless an isometry invariant (the class profile modulo 2) excludes it. For each candidate not already in the orbit, an isometry of $L$ taking $x$ to $y$ exists if and only if some isomorphism of the rooted local Gram graphs of $x$ and $y$ extends to an isometry preserving $L$, and this is decided exactly: the isomorphisms from the graph of $x$ to that of $y$ form one coset $tau G_x$ of the automorphism group $G_x$ of the graph of $x$, where $tau$ is one of them, extended to a rational isometry taking $M_x$ to $M_y$. Then $tau^(-1)(L)$ is an overlattice of $M_x$ in $M_x^\#$, and an isometry of $L$ taking $x$ to $y$ exists if and only if some element of $G_x$ maps the subgroup $L slash M_x$ of $M_x^\# slash M_x$ to $tau^(-1)(L) slash M_x$. This transporter problem is decided by GAP in the same action of $G_x$ on $M_x^\# slash M_x$ that gives the stabilizer. Since every candidate is tested, the orbit of $plus.minus x$ is determined, not only bounded below, and $|"Aut" L| = 2 |"Aut"(L) dot (plus.minus x)| dot |"Aut"(L)_x|$. All returned generators are checked in integers. The computation assumes the completeness of the short-vector enumeration, that bliss returns all automorphisms of the graph and canonical labellings, the correctness of GAP's stabilizer and orbit computations, and that our code implements these steps (Table 10).

== 7.2 Isomorphism keys
Our lists identify objects (lattices, pairs $(L, h)$, double cosets) by _keys_: functions $kappa$ computed from a representative. A key is _invariant_ if isometric objects have equal keys. The lists are used through the following two lemmas, which need only invariance, never the converse.

#thm[
*Lemma 7.1.* Let $cal(S)$ be a finite set of isomorphism classes with weights $w(X) > 0$ and $sum_(X in cal(S)) w(X) = mu$, and let $kappa$ be an invariant key. If $X_1, dots, X_n in cal(S)$ have pairwise distinct keys and $sum_i w(X_i) = mu$, then $cal(S) = {X_1, dots, X_n}$.
]
#proofbox[
Distinct keys imply distinct classes, so $sum_i w(X_i) <= mu$, with equality only if every class of $cal(S)$ occurs. #qed
]

#thm[
*Lemma 7.2.* In the situation of Lemma 7.1, two objects of $cal(S)$ with equal keys are isomorphic.
]
#proofbox[
Let $X in cal(S)$. By Lemma 7.1, $X$ is isomorphic to some $X_i$, so $kappa(X) = kappa(X_i)$, and the $kappa(X_i)$ are pairwise distinct. So $kappa$ separates the classes of $cal(S)$. #qed
]

Lemma 7.2 is what allows the subtraction of Lemma 6.6: the Hermitian entries of order 4 of the lattices of $S_"gen"$, wherever they come from, are classes of the set of Proposition 4.6, so each is the unique class of that list with its key, and no explicit isometry is needed.

The keys are of the following kinds. For each we say how it is computed and why it is invariant.

+ _Canonical form of the generating shell._ Let $t(L)$ be the least even $t >= 2$ such that the vectors of norm at most $t$ generate $L$, and $X_t = {x in L : 0 < (x, x) <= t(L)}$ (with both signs). The graph $Gamma(L)$ has vertex set $X_t$ and an arc $x -> y$ coloured by a fixed function of $(x, y)$ (and of $(x, omega y)$ for a lattice with a structure $omega$); arcs with $|(x, y)|$ below a threshold depending only on $|X_t|$ are omitted. It is turned into a vertex-coloured graph by one layer per occurring colour, bliss computes a canonical labelling, and the key is $(t(L), |X_t|)$ and a hash of the relabelled arc set #calc("C5"). _Invariance:_ $t(L)$ and $X_t$ are defined by norms and generation only, so an isometry $u : L_1 -> L_2$ (commuting with $omega$) maps $X_t (L_1)$ onto $X_t (L_2)$; it preserves inner products, hence colours and omitted arcs, so it is an isomorphism of coloured graphs; the layer of a colour is its rank among the occurring colours, which are the same on both sides, so $u$ extends to an isomorphism of the layered graphs; a canonical labelling maps isomorphic graphs to the same graph, and the hash is a function of that graph. The choice of basis, the order of the vectors and the representatives of $plus.minus x$ only renumber vertices. A non-injective colour function or hash only makes the key coarser. With a structure $omega$, the key is an invariant of the Hermitian lattice up to isometries commuting with $omega$; it is used only for lists of Hermitian lattices up to such isometries (the lists of Appendix B.4), where a lattice and its complex conjugate are not identified.
+ _Configuration of the vectors of norm at most 4_: the same construction with $t(L)$ replaced by 4. It is invariant by the same argument; since these vectors need not generate $L$, it is used only to show non-isometry.
+ _The invariant $D$ of $(L, plus.minus g)$_, $g^2 = -1$. With $V = L slash pi.alt L tilde.equiv FF_2^16$, let $n(c)$ and $r(c)$ be the numbers of lines $chevron.l plus.minus x, plus.minus x g chevron.r$ of norm 4 and of norm 2 in the class $c in V$ ($r = 0$ for the rootless lattices of Proposition 4.6; the same function is used for the lattices with roots met on the way), form the convolutions $n * n$, $n * n * n$ and $r * n$ over the group $V$, colour each $c$ by a fixed function of the five values at $c$, and let $D$ be a hash of the multiset of colours. _Invariance:_ an isometry $phi : L -> L'$ with $phi g = plus.minus g' phi$ maps lines to lines and $pi.alt L$ to $pi.alt L'$, so it induces a group isomorphism $V -> V'$ carrying $n, r$ to $n', r'$; convolutions commute with group isomorphisms; and $x$, $-x$, $x g$, $-x g$ lie in the same class modulo $pi.alt$, so the counts do not depend on representatives. Replacing $g$ by $-g$ does not change $pi.alt L$ or the lines.
+ _Canonical Gram matrix of a pair $(L, plus.minus g)$_, $g^2=-1$ (used together with $D$ for Proposition 4.6). Select the classes of $V$ whose colour occurs once and which contain exactly one Gaussian line of norm 4, ordered by colour. Try the first $16,20,24,dots$ lines, up to 64, and choose the first such prefix that generates $L$ as a $ZZ[i]$-module. Generation means that the vectors $x_j,x_j g$ generate the integer lattice; it is checked by rank and Hermite normal form. If such a prefix exists, the pair is _discrete_. Whether a pair is discrete depends only on the colours of the classes of $V$, the numbers of lines in them and generation, which an isometry or a complex conjugation preserves as for $D$; so it is an invariant, and the key of a pair, its canonical Gram matrix if it is discrete and $D$ otherwise, is an invariant as well.
  Starting with line 0, repeatedly select the least indexed unprocessed line adjacent to a processed line under the relation $H(x_j,x_k)!=0$. Its parent is the first adjacent line in the order in which the processed lines were reached. Rotate its representative by the unique unit for which the inner product with the parent has positive real part and non-negative imaginary part. Apply the same procedure to the conjugate Hermitian form. The key is the smaller hash of the two resulting Gram matrices.
  _Connectedness:_ if the non-orthogonality graph were disconnected, its components would generate a non-trivial orthogonal decomposition of $L$. Each summand would be an even unimodular integer lattice without roots. Its rank is a multiple of 8, and a non-trivial partition of 32 has a part of rank 8 or 16, where an even unimodular lattice has a root. Thus the graph is connected.
  _Invariance and the group bound:_ changing representatives does not change the traversal or its parents. Inductively, the quadrant condition makes all normalized representatives agree up to the same unit on the first line, so the Gram matrix is unchanged. A linear isometry preserves this construction, and an anti-linear one exchanges the two constructions. Each element of $"U"$ preserves every selected line. Its units on two adjacent lines coincide, and connectedness and generation give $"U"={plus.minus 1,plus.minus i}$. If the hashes of the two Gram matrices (of the form and of its conjugate) differ, the pair has no anti-linear automorphism, since one would exchange the two constructions. Equality of the Gram matrices themselves would establish one by mapping the generating vectors; equality of hashes alone is not used for that conclusion. For all $"62,726"$ discrete classes in the final list these two hashes differ.
+ _Canonical form of a gluing_: the lexicographically smallest table of $phi$ in its double coset $overline(G)_K phi overline(G)_F$, computed by a search along stabilizer chains of $overline(G)_K$ and $overline(G)_F$ #calc("C4"). It is a function of the double coset by definition; that the program returns this minimum is a property of the program, which we assume (Table 10); it was checked on samples, and the canonical forms of all representatives were recomputed from scratch with the same result.
+ _Canonical form for the type $(8, 24)$_: the construction of the canonical form of the generating shell, applied to the 720 vectors of norm 2 of $K^\#$, with a separate call of bliss and a truncated hash #calc("C6")\; it is invariant by the same argument for all $"14,482"$ lattices, whether or not the 720 vectors generate $K^\#$.

All keys are compared only between values produced by the same program with the same version of bliss. In particular, the keys of the Hermitian entries of order 4 of the lattices of $S_"gen"$ #calc("C1") and the keys of the classes of Proposition 4.6 #calc("C12") are computed by the same function: the rule is the canonical Gram key for a discrete class and the invariant $D$ otherwise. The same version of the program was used for both. For each of the $"76,742"$ classes of Proposition 4.6, rerunning the enumeration used for its key (the same function on the saved Gram matrix of the underlying lattice of rank 32) returns exactly $"73,440"$ vectors, each of norm 4 and no two equal up to sign. That lattice is rootless (by PARI `qfminim`, Table 10), hence extremal, and so has $"146,880" = 2 dot "73,440"$ vectors of norm 4 (§B.3); every pair $plus.minus v$ of norm 4 is therefore represented, the enumeration is complete for these classes, and $D$ and the canonical Gram keys recomputed from it are the saved ones #calc("C12"). For all $"16,163"$ recorded Hermitian entries of order 4 of $S_"gen"$, the class of Proposition 4.6 with the same key has the same $|"U"|$ and the same value of $f$ #calc("C29"). The lists of §4, Appendix B and Proposition 4.6 use keys only through Lemmas 7.1 and 7.2, and every identification "equal keys, hence isometric" made while building a list concerns objects of the set $cal(S)$ to which Lemma 7.2 applies.

== 7.3 What the proofs assume

#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 10. What Theorems A and B assume. Outputs that are checked in integers (an element is an automorphism, a matrix is an isometry, a vector has the stated norm) are not assumptions; what is assumed is that the programs do not miss anything (a group is all of $"Aut"$, a list of vectors or of double cosets is complete, no further isometry exists) and that our own code implements the stated formulas.],
  table.header([assumption], [used for], [theorem]),
  [King's mass $m$ [Kin03], recomputed with our implementation of the same inclusion–exclusion #calc("C14")], [term $2m$], [A],
  [PARI `qfminim`: the enumeration of short vectors is complete], [minima, roots and norm-4 vectors of components and of glued lattices; extremality; the numbers $g(W_3)$ of §B.3 and $g(Z_+)$ of (B.2); the absence of roots in the $"76,742"$ classes of Proposition 4.6 (on which the count of §7.2 rests); the theta coefficients of the bases of the six linear programs (§B.1.1)], [A, B],
  [PARI `qfauto`: the group generated by the returned automorphisms is all of $"Aut"$], [$|"Aut"|$ and $|"U"|$ of components; orbit reduction in Lemma B.1 (a group that is too small only adds work)], [A, B],
  [PARI `qfisom`: the answer "not isometric" is correct], [merging in the enumeration of Lemma B.1 for ranks at most 18 and in the classifications of genera by neighbours; for rank 20 every "not isometric" was confirmed by different keys (§7.2)], [A, B],
  [bliss: the canonical labelling is canonical, and the automorphism group returned is all of $"Aut"$], [the canonical forms of §7.2 of the generating shell, of the configuration of vectors of norm at most 4 and for the type $(8, 24)$; $|"Aut" K|$ for $"14,462"$ lattices of type $(8, 24)$; the automorphisms of the rooted local Gram graphs of §7.1, hence $|"Aut" L|$ for the $"12,546"$ lattices of $S_"gen"$ and the two exceptional known lattices], [A, B],
  [our code: the full-group algorithm of §7.1 (extension of graph automorphisms, the stabilizer of $L slash M$, the test of every candidate line)], [$|"Aut" L|$ for $S_"gen"$ and the exceptional known lattices], [A, B],
  [GAP: orders of permutation groups, double cosets, conjugacy classes, centralizers, set stabilizers], [$|O(q)|$, $|"Stab"(phi)|$, double cosets, conjugacy classes and invariant group signatures; the $"Aut" P$-orbits of good $W$ and $W_3$ of §B.3 (complete) and their stabilizers, hence the left side of (B.1); the classes with $g^2 = -1$ and their centralizers in $W(F_4)$ and $W(E_8)$ (Proposition C.1); the stabilizers of §7.1], [A, B],
  [our canonical form of gluings returns the minimum over the double coset], [the lists of double cosets for involutions, and for odd prime order where GAP did not finish (§4.2) (the canonical form of a gluing, §7.2)], [A, B],
  [SageMath: the mass of a genus of $ZZ$-lattices (Smith–Minkowski–Siegel); dimensions of spaces of modular forms], [completeness of the classifications of genera of components $F$ (Appendix B.4) and the check in rank 16 (Appendix B.2); the bases of the six linear programs (§B.1.1)], [A, B],
  [published classifications [Hen13, Thm 1, Table 4, Rem. 3(a)], [Jür15, Satz 3.6.9]], [types $3 - (16, 0) - 0$ and $17 - (2, 0) - 0$; both are also confirmed by rootless masses (§4.2, Appendix D)], [A, B],
  [our implementations of the mass formulas of Appendices C and D, in exact arithmetic], [$M_2$; the rootless component masses of Table D1 and §D.4; the vanishing rootless Hermitian masses of reason 3 of Proposition 3.3], [A, B],
  [the weights $f slash |"U"|$ of Proposition 4.6: canonical forms, backtrack search, bliss], [$M_2$ confirmed by enumeration; $|"Aut" L|$ for the lattices of $S_"CQ"$ with cyclic group, and the $"2,547"$ classes of the quaternion lattices (Lemma 6.6(iii))], [A, B],
  [$|"Aut" K|$ for 20 lattices of type $(8, 24)$ whose 720 short vectors of $K^\#$ do not generate it: the stabilizer of $N slash N^prime$ in the automorphism group of the graph of the 720 vectors (16, §B.3), `qfauto` (4)], [the term of type $(8, 24)$], [A, B],
  [the short-vector enumeration for the Gaussian keys (floating-point Fincke–Pohst) is complete], [the keys $D$ and canonical Gram matrices of §7.2 for the entries of $S_"gen"$ (for the classes of Proposition 4.6 the enumeration is complete by the count of §7.2)], [A, B],
  [the canonical key of type $(8, 24)$ (§7.2) is invariant and its saved prefix is unique among the $"14,482"$ lattices $K$], [completion of four identifiers of type $(8, 24)$ in groups of order above $10^7$ (§6.3)], [A, B],
  [our code: exact identifiers of entries and exact types (§6.3), minimum-shell profiles, orientation and exact sums], [classification, Gaussian weights and (5.1)], [A, B],
)

== 7.4 Computations and reproduction
#calc-table()
#let release = json("release-files.json")
#let cpu = release.cpu
This paper is accompanied by the repository #link(release.repository)[#raw(release.repository.replace("https://", ""))], version #raw(release.at("version", default: "")). The directory `certificates/` holds the certificates read by the checks; Table 11 names, for each computation, its certificate and its place in the research record [Rec]. The directory `data/` holds the classified lattices of $S$ with the orders and structure of their automorphism groups. The program `verify/verify.py` needs only Python and the files of the repository. It checks the six rational certificates of §B.1.1; identity (5.1) from King's mass and the classified lattices, with the totals of Theorems A and B and Table 8; every term of (C.2) from its coefficient, unitary group order and embedding average, with every embedding average $A_16 (T)$ recomputed from the saved local statistics of §C.6, and the coefficients of $A_n$, $D_n$ and $G_m$ from the closed formulas of §C.3 (the four exceptional coefficients are inputs); in integers, the generators of the two largest groups and the 106 isometries of §6.3; the counts of §6 against the classified lattices, and from the saved rows of the full-group records the classes of Lemma 6.3(i) with a witness row for each, the identity of Lemma 6.3(ii) for every list with the distinctness of the rows that carry identifiers, and the subtraction of Lemma 6.6 as differences of sets of entries; the completeness identities from their terms, that is (4.1) for the pairs of components of §4.1 and §4.2, with $|O(q_F)|$ from the closed formulas of Lemma 2.1 (term by term for the 83 pairs of involution components and for 258 pairs of odd prime order, from the numbers of decompositions for the four pairs of $5 - (6, 8) - 6$ and from the recorded masses for the two pairs of $3 - (6, 20) - 6$ without extremal gluings, §4.2), the mass identities of Table B4, the steps of (B.2) used for Table B4 and their chains back to lists $(z, z)$ that reach the rootless mass, identity (B.1) (for $A_2^(12)$ from the terms of its $"106,624"$ orbits of good $W_3$) and the sum of Proposition 4.6 with the number of weights obtained by each method; and the finite checks of the 9-torsion reduction and of the orbits of elementary quadratic 2-forms (§D.4.6). It does not recompute the terms of these identities (automorphism groups and stabilizers, numbers of decompositions, masses of genera), the rows of the full-group records, or isometry decisions; these rest on the statements of Table 7 and the assumptions of Table 10.

The programs that produced the lists, their runs and our verification notes are kept unchanged in the research record [Rec], a release asset described by its own index. They are a record of the work, not a library; documented modules will replace them in later versions. Table 12 lists the release assets with their SHA-256 digests. Files taken from other sources (King's table [Kin03], the pages of [NS]) are not included.

The CPU times of Table 11 add up to #if cpu.lower_bound [at least] #cpu.total_text hours. The #str(cpu.excluded.len()) computations whose time was not recorded (#cpu.excluded.join(", ")) are not included.

#ptable(
  columns: (auto, auto, 1fr), align: (left, right, left),
  caption: [Table 12. Files attached to the release, with their SHA-256 digests.],
  table.header([file], [size], [content]),
  ..release.files.map(f => (
    raw(f.name),
    f.size,
    [#(if release.at("markup", default: false) { eval(f.content, mode: "markup") } else { f.content }) #linebreak() #text(size: 0.72em, raw(f.sha256))],
  )).flatten(),
)

= Appendix A. The previously known lattices <appA>

#policy-box[
*What this appendix does.* It locates in $S$ the 27 lattices of $S$ that appear in the literature. The identifications with the literature are a comparison and are not used in the proofs of Theorems A and B; the published descriptions of 26 of these lattices (codes, Hermitian and quaternionic structures; the catalogue [NS] was used only for the comparison below) are among the starting points of the records of §6.3, which the proofs use through Lemma 6.3. We state on what each identification rests, using four grades of decreasing strength.
]

The grades are:
(1) an explicit isometry to a lattice given by the source, checked in integers ($T^sans(T) G_1 T = G_2$, $det T = plus.minus 1$);
(2) the same construction: the lattice is built from a code that is shown to be equivalent to the code of the source by canonical forms, and equivalent codes give isometric lattices;
(3) forced by two complete classifications: the source classifies a family completely, our lists of §4 contain the same family completely, and an invariant that takes distinct values within the family matches value by value;
(4) agreement of invariants only.

All 27 lattices lie in $S_"gen"$ of Theorem B: 25 of them have $|"Aut" L|$ divisible by an odd prime, and the remaining two have 2-groups of orders $2^18$ and $2^14$, larger than the orders allowed in $S_(C_2^2)$ and $S_"CQ"$. Their 27 values of $|"Aut" L|$ are pairwise distinct #calc("C32", "C33"). The lattice of $S$ with $|"Aut" L| = 2^17 3^4 5$, the tenth of the ten with $|"Aut" L| > 10^7$ (§6.3), is not among them. Of the twelve rows from [KN93], the lattices $S_1$ and $S_2$ have neighbour defect 0 and are two of the five code lattices of [KV89, Satz 10]; the other ten, $S_3$, $C_1, dots, C_5$ and $G_1, dots, G_4$, are exactly the extremal lattices of neighbour defect 8 [KN93, Thm 1].

#ptable(
  columns: 3, align: (left, right, left),
  caption: [Table A1. The 27 previously known lattices of $S$, by decreasing $|"Aut" L|$, with the grade of each identification.],
  table.header([lattice], [$|"Aut" L|$], [identification (grade)]),
  [$"BW"_32$], [$2^31 3^5 5^2 7 dot 17 dot 31$], [[NS] BW32, dim32f1 (1); the code lattice of $"RM"(2, 5)$, one of the five codes of [KV89, Satz 10], by definition (2); [Neb98q, §21] (via dim32f1); [Hen13, Table 4, no. 70] (3); [Jür15, Satz 3.6.9, $L_4$] (3); [BN97, Thm 5.1, $Lambda''_32$] (4); [BCDYY25, §3.3] (4)],
  [code lattice of one of the five doubly even self-dual codes of length 32 and minimum weight 8 [KV89, Satz 10]], [$2^31 3^5 5^2 7$], [[NS] KV32F (1); [BCDYY25, §3.3] (4)],
  [[KN93] $S_1$], [$2^28 3^2 7$], [[NS] KV32U (1); [KN93, Table 1] (2); [BCDYY25, §3.3] (4)],
  [[KN93] $S_2$], [$2^25 3^2 5$], [[NS] KV32G (1); [KN93, Table 1] (2); [BCDYY25, §3.3] (4)],
  [$Lambda'_32$], [$2^13 3^5 5^4$], [[NS] BN32q, dim32f3 (1); [Neb98q, §21] (via dim32f3); [Hen13, no. 71] (3); [BN97, Thm 5.1] (3); [BCDYY25, §3.3] (4)],
  [code lattice of the extended quadratic residue code of length 32], [$2^21 3 dot 5 dot 31$], [[NS] KV32QR (1); (2) by definition; [BCDYY25, §3.3] (4)],
  [$Lambda_32$], [$2^17 3^4 7$], [[NS] BN32p (1); [Hen13, no. 72] (3); [BN97, Thm 5.1] (3)],
  [[KN93] $C_1$], [$2^21 3 dot 7$], [[NS] KV32C1 (1); [KN93, Table 1] (2)#super[†]],
  [[KN93] $S_3$], [$2^16 3^4 5$], [[NS] KV32S3 (1); [KN93, Table 1] (2); [Hen13, no. 73] (3)],
  [[Hen13] no. 74], [$2^13 3^3 7$], [(3)],
  [[KN93] $G_4$], [$2^16 3^2$], [[NS] KV32G4 (1); (2)],
  [[KN93] $C_2$], [$2^17 3$], [[NS] KV32C2 (1); (2)],
  [[Neb98q] $[4.L_3(4).2^2]_(32,2)$], [$2^10 3^2 5 dot 7$], [[NS] dim32f5 (1); [Neb98q, §21, Thm 21.1] (via dim32f5)],
  [[KN93] $G_3$], [$2^15 3^2$], [[NS] KV32G3 (1); (2)],
  [[KN93] $C_3$], [$2^18$], [[NS] KV32C3 (1); (2)],
  [[KN93] $G_1$], [$2^14 3 dot 5$], [[NS] KV32G1 (1); (2)],
  [[KN93] $C_4$], [$2^16 3$], [[NS] KV32C4 (1); (2)],
  [[Jür15] $L_3$ = [Hen13] no. 75], [$2^7 3^3 17$], [[NS] L\_32.2, dim32f6 (1); [BB92, §2, $m = 51$] (1); [Neb98q, §21] (via dim32f6); (3) for both],
  [[Hen13] no. 76], [$2^8 3^4$], [(3)],
  [[Jür15] $L_2$], [$2^7 3^2 17$], [(3)],
  [[KN93] $C_5$], [$2^14$], [[NS] KV32C5 (1); (2)#super[†]],
  [[Hen13] no. 77], [$2^8 3^2 7$], [(3)],
  [[KN93] $G_2$], [$2^12 3$], [[NS] KV32G2 (1); (2); [BCDYY25, §3.3] (4)],
  [[Hen13] no. 78], [$2^8 3^2 5$], [(3)],
  [[Hen13] no. 79], [$2^10 3^2$], [(3)],
  [[Hen13] no. 80], [$2^6 3^4$], [(3)],
  [[Jür15] $L_1$], [$2^5 3^2 17$], [(3)],
)

#small[
Notes. [Hen13, Thm 1, Table 4, Rem. 3(a)] shows that there are exactly eleven extremal lattices with an automorphism of order 3 without fixed points, and gives the orders $|"U"|$ of their unitary groups; our list of type $3 - (16, 0) - 0$ has eleven lattices, its completeness is proved by the rootless Hermitian mass (§4.2), and $|C(h)| = |"U"|$ holds for all eleven with eleven distinct values. Our representatives of this type were constructed without the Gram matrices of [Hen13], by neighbour methods and from codes; [Hen13, Table 4] was used only to compare the orders $|"U"|$. [Hen13, Table 4] also names rows 70, 71, 72 as $"BW"_32 = Lambda''_32$, $Lambda'_32$, $Lambda_32$, which agrees with the identification through [BN97].

[Jür15, Satz 3.6.9] classifies a genus of Hermitian $ZZ[zeta_17]$-lattices of rank 2, shows that six Hermitian classes have extremal trace lattices, representing four integer-lattice classes, and gives the orders $|O(L_i)|$, which are distinct. Our four lattices of type $17 - (2, 0) - 0$ were obtained by Kneser neighbours inside that genus, starting from its lattice with minimum 2 and preserving the $ZZ[zeta_17]$-structure at each step #calc("C32")\; they are therefore lattices of the genus, our list is complete because the six Hermitian classes they carry have $sum 1 slash |"U"(Z)| = 1 slash 12$, the rootless mass of the genus (§4.2), and they are identified with $L_1, dots, L_4$ by $|"Aut" L| = |O(L_i)|$.

For [KN93], the twelve codes were reconstructed from Table 1 there and matched by canonical forms with the codes from which our lattices are built; the lattice is determined by the code [KN93, Thm 2(b)]. #super[†]For $C_1$ and $C_5$ the printed word $c$ does not give a code with the stated weight enumerator; after changing two coordinates it does, and the matching uses the corrected word. The column $|"Aut" C|$ of [KN93, Table 1] differs from our recomputation in six rows.

[BN97, Thm 5.1] classifies the eleven classes of rank-8 lattices over the Hurwitz order in the relevant genus, complete by the mass formula [BN97, Rem. 5.5]; our representatives of $Lambda_32$ and $Lambda'_32$ were constructed by neighbours in that genus, so they carry a structure over the Hurwitz order by construction, and the orders of their unitary groups are those in [BN97, Table 4]. This gives grade (3). For $"BW"_32$ we use only that it is one of the three extremal lattices with such a structure [BN97, Rem. 5.2].

[BCDYY25, §3.3] gives seven orders $|"Aut" L|$ of the 19 extremal lattices of rank 32 in the Magma database: those of the five code lattices, without saying which order belongs to which code; $2^13 3^5 5^4$, the only other order above $2^21 3 dot 5 dot 31$; and the least order, $2^12 dot 3$. In Table A1 they are the orders of the five code lattices, of $Lambda'_32$ and of $G_2$ of [KN93]; this is why the comparison with [BCDYY25] is only at grade (4).

An entry "[NS] _name_ (1)" means that the Gram matrix on the page _name_ of the catalogue [NS] was mapped onto our lattice by an isometry found with bliss and checked in integers #calc("C33")\; the 25 extremal pages of rank 32 of the catalogue give 19 distinct lattices, all in Table A1, which names 22 of the pages; the other three are further pages of $"BW"_32$, and eight lattices of the table have no page there. The catalogue names the pages dim32f1, dim32f3, dim32f5, dim32f6 by the groups $[2^(1+10)_+ . O^+_(10)(2)]_32$, $[((SL_2(5) compose SL_2(5)) times.square (SL_2(5) compose SL_2(5))) : S_4]_(32,2)$, $[4.L_3(4).2^2]_(32,2)$, $[SL_2(17) square.stroked S_3]_(32,1)$ of [Neb98q, §21]; there these are listed among the rational irreducible maximal finite subgroups of $GL_32(QQ)$ with an invariant lattice $L$ of determinant 1 and minimum 4, and $G = "Aut" L$ [Neb98q, Thm 21.1 and its proof]. Of these, $[4.L_3(4).2^2]_(32,2)$ appears in no other source we know.

[BB92, §2] gives, for $m = 51$, the Gram matrix of a unimodular trace form lattice over $ZZ[zeta_51]$ by its first row; it is isometric to our lattice $L_3$ of [Jür15] #calc("C33").
]

= Appendix B. Component lists <appB>

#policy-box[
*What this appendix does.* §B.1 lists the quadruples excluded in Proposition 3.3. §B.2 proves that the components of the involution types of ranks 8 to 20 are reached by chains of index-2 overlattices, which makes their enumeration complete. §B.3 reduces the type $(8, 24)$ to three Niemeier lattices and proves the list of $"14,482"$ lattices complete. §B.4 treats the components of the types of odd prime order.
]

== B.1 Excluded types of odd prime order
#ptable(
  columns: 2, align: (left, left),
  caption: [Table B1. The 74 quadruples $p - (z, f) - s$ excluded in the proof of Proposition 3.3, with the reason. The density exclusions are grouped by $p$ #calc("C2").],
  table.header([type], [reason]),
  [$3 - (2, 28) - 0$], [no genus],
  [$3 - (6, 20) - 0$], [no genus],
  [$3 - (7, 18) - 5$], [no rootless $U$],
  [$3 - (7, 18) - 7$], [no rootless $U$],
  [$3 - (8, 16) - 6$], [linear programming],
  [$3 - (9, 14) - 5$], [linear programming],
  [$3 - (10, 12) - 0$], [no genus],
  [$3 - (11, 10) - 7$], [linear programming],
  [$3 - (11, 10) - 9$], [linear programming],
  [$3 - (12, 8) - 6$], [linear programming],
  [$3 - (14, 4) - 0$], [no genus],
  [$3 - (14, 4) - 4$], [no genus],
  [$7 - (2, 20) - 0$], [no genus],
  [$7 - (3, 14) - 3$], [linear programming],
  [$11 - (2, 12) - 0$], [no genus],
  [$p = 3$: $(1, 30) - 1$, $(2, 28) - 2$, $(3, 26) - 1$, $(3, 26) - 3$, $(4, 24) - 0$, $(4, 24) - 2$, $(4, 24) - 4$, $(5, 22) - 1$, $(5, 22) - 3$, $(5, 22) - 5$, $(6, 20) - 2$, $(6, 20) - 4$, $(7, 18) - 1$, $(7, 18) - 3$, $(8, 16) - 0$, $(8, 16) - 2$, $(8, 16) - 4$, $(9, 14) - 1$, $(9, 14) - 3$, $(10, 12) - 2$, $(10, 12) - 4$, $(11, 10) - 1$, $(11, 10) - 3$, $(11, 10) - 5$, $(12, 8) - 0$, $(12, 8) - 2$, $(12, 8) - 4$, $(13, 6) - 1$, $(13, 6) - 3$, $(14, 4) - 2$, $(15, 2) - 1$], [density],
  [$p = 5$: $(1, 28) - 1$, $(2, 24) - 0$, $(2, 24) - 2$, $(3, 20) - 1$, $(3, 20) - 3$, $(4, 16) - 0$, $(4, 16) - 2$, $(5, 12) - 1$, $(5, 12) - 3$, $(6, 8) - 0$, $(6, 8) - 2$, $(7, 4) - 1$], [density],
  [$p = 7$: $(1, 26) - 1$, $(2, 20) - 2$, $(3, 14) - 1$, $(4, 8) - 0$, $(4, 8) - 2$, $(5, 2) - 1$], [density],
  [$p = 11$: $(1, 22) - 1$, $(2, 12) - 2$, $(3, 2) - 1$], [density],
  [$p = 13$: $(1, 20) - 1$, $(2, 8) - 0$, $(2, 8) - 2$], [density],
  [$p = 17$: $(1, 16) - 1$], [density],
  [$p = 19$: $(1, 14) - 1$], [density],
  [$p = 23$: $(1, 10) - 1$], [density],
  [$p = 29$: $(1, 4) - 1$], [density],
)

#include "finite-type-exclusions.typ"

== B.2 Components of involutions

Let $cal(M)(d)$ be the set of isometry classes of even 2-elementary lattices $M$ of rank $d$ with $min M >= 4$.

#thm[
*Lemma B.1.* Let $M in cal(M)(d)$ with $M^\# slash M tilde.equiv FF_2^s$. There are a unimodular lattice $U$ of rank $d$ without vectors of norm 1 and a chain $sqrt(2) U = Y_0 subset Y_1 subset dots subset Y_m = M$, $m = (d - s) slash 2$, with $[Y_(j+1) : Y_j] = 2$, such that every $Y_j$ lies in $cal(M)(d)$ and $Y_(j+1) = Y_j + ZZ v_j$ for some $v_j in Y_j^\#$ with $(v_j, v_j) in 2 ZZ$.
]
#proofbox[
The proof of Proposition 3.2 gives $T = sqrt(2) U subset M$ with $U$ unimodular and $min U >= 2$, and $M slash T$ is elementary abelian of order $2^m$, since $2 M subset 2 M^\# subset T$. Take any chain of subgroups of $M slash T$ with successive quotients of order 2, and let $Y_j$ be their preimages. Each $Y_j$ is a sublattice of $M$, so it is even with minimum at least 4. Since $Y_j supset T$, $Y_j^\# subset T^\# = U slash sqrt(2)$ and $2 Y_j^\# subset sqrt(2) U = T subset Y_j$, so $Y_j$ is 2-elementary. Finally $Y_(j+1) subset M subset M^\# subset Y_j^\#$, so $Y_(j+1) = Y_j + ZZ v_j$ with $v_j in Y_j^\#$, and $(v_j, v_j) in 2 ZZ$ because $M$ is even. #qed
]

So $cal(M)(d)$, or its part with $s >= s_0$, is obtained as follows: start from the lattices $sqrt(2) U$, where $U$ runs over the unimodular lattices of rank $d$ without vectors of norm 1 [SPLAG, Table 16.7], and at each step replace every $Y$ by the lattices $Y + ZZ v$, where $v + Y$ runs over representatives of the $"Aut" Y$-orbits of elements of order 2 of $Y^\# slash Y$ with $(v, v) in 2 ZZ$, keeping those with minimum at least 4 and merging isometric ones. Since $s$ drops by 2 at each step, the part with $s >= s_0$ is reached by chains staying in $s >= s_0$. The generators of $"Aut" Y$ used for the orbits are checked to be automorphisms in integers, and merging uses explicit isometries checked in integers; two lattices are kept apart only if `qfisom` finds no isometry between them; in rank 20 this answer was also confirmed by different isometry invariants (§7.2). Hence no class is lost, and the list has no repetition #calc("C5", "C7").

#ptable(
  columns: 5, align: (left, left, right, right, left),
  caption: [Table B2. The components of involution types. "Classes" is the number of classes in $cal(M)(d)$ with the given $(s, delta)$ #calc("C5", "C6", "C7").],
  table.header([rank $d$], [$(s, delta)$], [classes], [$|"Aut" M|$], [used in type]),
  [8], [$(8, 0)$], [1], [$"696,729,600"$ ($sqrt(2) E_8$)], [$(8, 24)$],
  [12], [$(12, 1)$], [1], [$"980,995,276,800"$ ($sqrt(2) D_12^+$)], [$(12, 20)$],
  [14], [$(14, 1)$], [1], [$"16,855,282,483,200"$ ($sqrt(2) (E_7^2)^+$)], [$(14, 18)$],
  [15], [$(15, 1)$], [1], [$"41,845,579,776,000"$ ($sqrt(2) A_15^+$)], [$(15, 17)$],
  [16], [all], [10], [Table 6], [$(16, 16)$],
  [17], [$(15, 1)$], [1], [$"8,957,952,000"$], [$(15, 17)$],
  [18], [$(14, 1)$], [5], [$"7,077,888"$, $"11,197,440"$, $"477,757,440"$, $"566,231,040"$, $"2,717,908,992"$], [$(14, 18)$],
  [20], [$(12, 1)$], [66], [from $64$ to $"79,272,345,600"$], [$(12, 20)$],
  [24], [$(8, 0)$, $min M^\# >= 2$], [$"14,482"$], [from $2$ to $"253,671,505,920"$], [$(8, 24)$],
)

For ranks 8, 12, 14, 15 and 16 the enumeration was run for all $s$; it shows that $sqrt(2) E_8$, $sqrt(2) D_12^+$, $sqrt(2) (E_7^2)^+$ and $sqrt(2) A_15^+$ are the only lattices in $cal(M)(d)$ for $d = 8, 12, 14, 15$, and gives the ten classes of rank 16. Since $F$ and $K$ have the same $(s, delta)$, this determines which $(s, delta)$ are needed in ranks 17, 18, 20 and 24; in ranks 17, 18 and 20 the enumeration was stopped at $s_0 = 15, 14, 12$ respectively. In rank 20 it gives $297$ classes with $s >= 12$, of which 66 have $s = 12$, $delta = 1$. Rank 24 is treated in §B.3.

Two checks, not used in the proofs, compare these lists with computations that do not use chains. In rank 16, a classification of all seventeen genera of even 2-elementary lattices by Kneser 3-neighbours, completed when $sum 1 slash |"Aut"|$ reached the mass of the genus, contains exactly the ten classes of minimum at least 4 #calc("C7"). For the 57 genera of rank at most 18 and for $(20, 12, 1)$, $sum 1 slash |"Aut" M|$ over the list equals, as a fraction, the rootless mass of the genus computed by King's inclusion–exclusion over root systems [Kin03] with the local densities at 2 for 2-elementary lattices #calc("C5"). For $(20, 12, 1)$ this second computation found the same 66 classes, matched one-to-one with the chain enumeration by explicit isometries.

== B.3 The type $(8, 24)$

Let $sigma$ be an involution of type $(8, 24)$ of an even unimodular lattice $L$ of rank 32, with $F = L^sigma tilde.equiv sqrt(2) E_8$ (§B.2) and $K = L^(-sigma)$. Then $K$ is even, of rank 24, 2-elementary with $s = 8$ and $delta = 0$, so $(y, y) in ZZ$ for $y in K^\#$.

#thm[
*Proposition B.2.* $L$ is extremal if and only if $min K >= 4$ and $min K^\# >= 2$.
]
#proofbox[
If $L$ is extremal, $min K >= 4$. By the criterion of §4.1.3, $L$ has a root if and only if some class $a in F^\# slash F$ and its image $phi(a)$ both contain a vector of norm 1. Since $F^\# = E_8 slash sqrt(2)$, a class of $F^\# slash F = E_8 slash 2 E_8$ contains a vector of norm 1 if and only if it is the class of a root of $E_8$, that is, if and only if $q_F (a) equiv 1$. A vector $y in K^\#$ of norm 1 lies in a class $b$ with $q_K (b) equiv 1$, so $phi^(-1)(b)$ has $q_F equiv -1 equiv 1$ and contains a vector of norm 1. Hence $L$ has a root if and only if $min K < 4$ or $K^\#$ has a vector of norm 1; and $K^\#$ has no non-zero vector of norm less than 1 (§4.1.3). #qed
]

Let $cal(K)$ be the set of isometry classes of such $K$. Let $T' subset K^\# slash K tilde.equiv FF_2^8$ be a maximal totally singular subspace for $q_K$; it has dimension 4, and $P = K + T'$ is an even unimodular lattice of rank 24, a Niemeier lattice. There are $270$ such $T'$, and each of the $135$ non-zero singular elements lies in $270 dot 15 slash 135 = 30$ of them, since $O(q_K)$ is transitive on these elements.

#thm[
*Theorem B.3.* For every $K in cal(K)$, $K^\#$ has exactly $720$ vectors of norm 2, the $270$ Niemeier lattices $K + T'$ have on average exactly $80$ roots, and hence at least one of them is the Leech lattice, $A_1^(24)$ or $A_2^(12)$.
]
_Strategy of the proof._ The lattices $K$ of rank 24 cannot be listed by a direct search, so we look at the 270 Niemeier lattices $K + T'$ that contain them. The vectors of norm 4 of $L$ form a spherical 7-design, and this fixes the number of vectors of norm 2 of $K^\#$. That number fixes the average number of roots of the Niemeier lattices $K + T'$. A Niemeier lattice with at most this many roots is the Leech lattice, $A_1^(24)$ or $A_2^(12)$, so every $K$ lies in one of these three lattices, where the search is finite.
#proofbox[
By Proposition B.2, $L$ is extremal. Let $X$ be the set of its vectors of norm 4; there are $"146,880"$ of them, the coefficient of $q^2$ in the theta series of an extremal lattice of rank 32, which is the unique modular form of weight 16 for $"SL"_2(ZZ)$ of the form $1 + O(q^2)$. By [BV01, Cor. 4.1], $X$ is a spherical 7-design. For $x in X$ write $x = x_+ + x_-$ with $x_+ in F^\#$, $x_- in K^\#$. For a 7-design on the sphere of radius 2 in $RR^32$ and the 8-dimensional space $F times.o RR$, the averages of the polynomials $|x_+|^(2 j)$, $j <= 3$, are those of the uniform measure:
$ sum_(x in X) (x_+, x_+)^j = |X| dot 4^j product_(i = 0)^(j - 1) (8 + 2 i) / (32 + 2 i). $
Since $(x_+, x_+) in {0, 1, 2, 3, 4}$ and the 240 vectors with $(x_+, x_+) = 4$ are those of $F$, the four equations $j = 0, dots, 3$ determine the numbers $N_a$ of $x$ with $(x_+, x_+) = a$: $N_0 = "12,240"$, $N_1 = "122,880"$, $N_2 = "11,520"$, $N_3 = 0$. A vector $y in K^\#$ of norm 2 is not in $K$ and lies in a class $b$ with $q_K (b) = 0$, so $phi^(-1)(b)$ is a class of $E_8 slash 2 E_8$ of vectors of norm 4, which contains exactly 16 vectors of norm 2 of $F^\#$ (in $E_8$, two vectors of norm 4 that are congruent modulo $2 E_8$ are equal up to sign or orthogonal, since their difference lies in $2 E_8$ and so has norm 0 or at least 8; hence a class contains at most 16 of them. The same comparison of norms shows that two congruent roots are equal up to sign and that no vector of norm 4 is congruent to 0 or to a root. The zero class and the $120$ classes $plus.minus r$ of roots leave $135$ of the $256$ classes, and $"2,160" = 135 dot 16$ forces exactly 16 vectors of norm 4 in each); each gives a vector $x_+ + y in X$. So $N_2 = 16 n$, where $n$ is the number of vectors of norm 2 in $K^\#$, and $n = 720$.

The roots of $P = K + T'$ are the vectors of norm 2 of $K^\#$ lying in the classes of $T'$, since $K$ has none. Each of the 720 vectors lies in a non-zero singular class and is counted in 30 of the $T'$, so the average number of roots is $720 dot 30 slash 270 = 80$. The number of roots of a Niemeier lattice is $24 h$, where $h$ is the Coxeter number of its root system, and the only values at most 80 are $0$, $48$ and $72$, attained only by the Leech lattice, $A_1^(24)$ and $A_2^(12)$ [SPLAG, Table 16.1]. #qed
]

Conversely, for a Niemeier lattice $P$ and a 4-dimensional subspace $W subset P slash 2 P$, totally singular for $x |-> (x, x) slash 2 mod 2$, with preimage $tilde(W) subset P$, put $K_W = {x in P : (x, tilde(W)) subset 2 ZZ}$. Then $K_W^\# = P + tilde(W) slash 2$, so $K_W$ is even, 2-elementary with $s = 8$ and $delta = 0$, and every $K in cal(K)$ is $K_W$ for $P = K + T'$ and $W = 2 K^\# slash 2 P$. Moreover $min K_W >= 4$ if and only if no root of $P$ is orthogonal modulo 2 to $tilde(W)$, and $min K_W^\# >= 2$ if and only if every non-zero class of $W$ has minimum at least 8 in $P$. We call such $W$ _good_.

*Enumeration.* For $P$ the Leech lattice and $A_1^(24)$, all $"Aut" P$-orbits of good $W$ were enumerated, subspace by subspace of dimension $1, dots, 4$, with an orbit-canonical choice at each step; there are $"1,244"$ and $"376,306"$ orbits, and they give $"1,216"$ and $"13,231"$ classes of $K$ with a $T'$ of $0$ or $48$ roots #calc("C6"). By Theorem B.3 every other $K in cal(K)$ has a $T'$ with $K + T' tilde.equiv A_2^(12)$. For $P = A_2^(12)$ we use the identity
#tagged("B.1", $ sum_([W]) 1 / (|"Stab"_("Aut" P) (W)|) = sum_(K in cal(K)) (a(K, P)) / (|"Aut" K|), $)
where the left side runs over the $"Aut" P$-orbits of good $W$ and $a(K, P)$ is the number of $T'$ with $K + T' tilde.equiv P$. It holds because $(P, W) |-> (K_W, P slash K_W)$ is a bijection between the isometry classes of the two kinds of pairs which preserves automorphism groups, and both sides are the orbit sums of this set of pairs. Call a 3-dimensional subspace $W_3 subset P slash 2 P$ good if it is totally singular and every non-zero class of $W_3$ has minimum at least 8 in $P$; the condition on roots is not imposed, since it does not pass to subspaces. The left side was computed from the $"106,624"$ orbits of good $W_3$, enumerated as for the Leech lattice: for each, $g(W_3)$ is the number of non-zero singular $x in W_3^perp slash W_3$ such that every class of $x + W_3$ has minimum at least 8 and every root of $P$ orthogonal modulo 2 to $tilde(W)_3$ has odd inner product with $x$, that is, the number of good $W supset W_3$. Since every hyperplane of a good $W$ is a good $W_3$, counting the pairs $W_3 subset W$ gives $sum_([W]) 1 slash |"Stab"(W)| = 1/15 sum_([W_3]) g(W_3) slash |"Stab"(W_3)| = "57,341,069,029" slash "69,120"$. The right side over the $"14,447"$ lattices already found, with $a(K, A_2^(12))$ counted over all 270 subspaces $T'$, falls short by $"429,583" slash 768$. Every lattice not yet found has $a(K, A_2^(12)) >= 1$, so it contributes a positive term to the difference. A search adapted to $A_2^(12)$ (a good $W$ maps onto each of the twelve summands $A_2 slash 2 A_2$ of $A_2^(12) slash 2 A_2^(12)$, so it is described by twelve lines in $W$) found 35 further lattices, and their terms fill the difference exactly. Hence no lattice is missing, and $cal(K)$ has $"14,482"$ elements. How the 35 lattices were found plays no role in this argument. Identity (B.1) holds exactly also for the Leech lattice and $A_1^(24)$, where both sides were computed independently #calc("C6").

#ptable(
  columns: 4, align: (left, right, right, left),
  caption: [Table B3. The $"14,482"$ lattices $K$ of type $(8, 24)$, by the smallest number of roots of a Niemeier lattice $K + T'$ #calc("C6").],
  table.header([$P$], [good $W$ (orbits)], [classes $K$], [completeness]),
  [Leech], [$"1,244"$], [$"1,216"$], [all orbits],
  [$A_1^(24)$], [$"376,306"$], [$"13,231"$], [all orbits],
  [$A_2^(12)$], [—], [$35$], [identity (B.1)],
  [total], [], [$"14,482"$], [],
)

The classes of $K$ are separated by the canonical form (§7.2) of the graph on the 720 vectors of norm 2 of $K^\#$, with inner products as colours. For $"14,462"$ of the lattices these vectors generate $K^\#$, so $"Aut" K$ is the automorphism group of the graph and the canonical form is a complete invariant. For the remaining 20, let $J^prime subset J = sqrt(2) K^\#$ be the sublattice spanned by the 720 vectors (of norm 4 in $J$); the automorphism group $A$ of the graph acts on $J^(prime\#) slash J^prime$, $"Aut" K$ is the stabilizer in $A$ of the subgroup $J slash J^prime$, so $|"Aut" K| = |A| slash |A dot (J slash J^prime)|$, and the orbit, written in the canonical coordinates, completes the invariant (16 lattices); for 4 lattices the orbit was too large and $"Aut" K$ was computed with `qfauto`; each of these four has an order $|"Aut" K|$ that no other of the $"14,482"$ lattices has, which separates them from each other and from the rest #calc("C6").

== B.4 Components of elements of odd prime order

Let $h$ have type $p - (z, f) - s$ with $p$ odd. We need all $F$ of minimum at least 4 in the genus of even $p$-elementary lattices of rank $f$ and determinant $p^s$, and all $Z$ of minimum at least 4 in the genus of Hermitian $cal(O)$-lattices of rank $z$ with $Z^\# slash Z tilde.equiv FF_p^s$ (Proposition D.1). Table B4 gives the result.

*The lattices $F$.* Except in two genera, the whole genus was classified by Kneser $q$-neighbours for a prime $q != p$, which stay in the genus; classes were identified by explicit isometries checked in integers and separated by invariants, and the classification is complete when $sum 1 slash |"Aut" F|$ over the classes found equals the mass of the genus given by the Smith–Minkowski–Siegel mass formula [SPLAG, Ch. 16] #calc("C2"). The lattices of minimum at least 4 were then selected. For the genera $(20, 6)$ ($p = 3$) and $(16, 4)$ ($p = 5$) only rootless lattices were searched for, and the list is complete when its total weight equals the rootless mass of the genus (Appendix D).

*The lattices $Z$.* If $s = z$, $Z$ is unimodular for $tilde(H)$ (Lemma 4.4), and the list of rootless $Z$ is complete when its total weight equals the rootless Hermitian mass of rank $z$, obtained from the mass of the genus and (4.2) (Appendix D). If $s < z$, let $ell$ be a line of the symplectic space $Z slash Z^*$ of Lemma 4.4 and $Z_+ subset Z$ the preimage of its orthogonal $ell^perp$. Then $[Z : Z_+] = p$, $Z_+ slash Z_+^* tilde.equiv ell^perp slash ell$, so $Z_+$ lies in the genus $(z, s + 2)$, and $Z = Z_+ + cal(O) y$ for an isotropic element $y + Z_+$ of the discriminant group of $Z_+$; if $Z$ has minimum at least 4, so has $Z_+$. Each $Z$ contains exactly $n = (p^(z - s) - 1) slash (p - 1)$ such $Z_+$. Let $g(Z_+)$ be the number of isotropic lines of the discriminant group of $Z_+$ that give a lattice of minimum at least 4. Counting pairs $Z_+ subset Z$ in two ways gives
#tagged("B.2", $ sum_Z n / (|"U"(Z)|) = sum_(Z_+) g(Z_+) / (|"U"(Z_+)|), $)
the sums over the lattices of minimum at least 4 in the two genera. Given the complete list for $(z, s + 2)$, the right side is computed from it, and the list for $(z, s)$ is complete when its left side reaches this value. For $3 - (13, 6) - 5$ and for $p = 5, 7$ the lists were completed instead against the rootless mass of the genus $(z, s)$ (Appendix D). For $p = 5$, (B.2) was also checked at every step of the chains that lead to the genera of Table B4, from $(6, 6)$ to $(6, 4)$ and from $(7, 7)$ through $(7, 5)$ to $(7, 3)$ (for $p = 7$ the genus $(4, 4)$ is itself the start of a chain, and the list for $5 - (8, 0) - 0$ was found directly, §4.2). For $z = 13$ ($p = 3$) it was checked only from $s = 13$ to $s = 11$ #calc("C9", "C2"). As a check not used in the proofs, for $p = 3$ each list of a genus $(z, s)$ with $s < z$ obtained by (B.2) also has $sum 1 slash |"U"(Z)|$ equal to the rootless Hermitian mass of the genus computed by the method of Appendix D.4 #calc("C16").

#ptable(
  columns: 5, align: (left, left, left, left, left),
  caption: [Table B4. Components of the 19 types of odd prime order: the genus $(f, s)$ of $F$ and the number of its classes of minimum at least 4 (in brackets, the number of classes of the whole genus when it was classified); the genus $(z, s)$ of $Z$ and the number of its classes of minimum at least 4. "Mass" means the mass of the genus, "rootless" the rootless mass of Appendix D, "(B.2)" the identity above, starting from the rootless unimodular lattices #calc("C8", "C9", "C2").],
  table.header([type], [$F$], [$Z$], [$F$ complete by], [$Z$ complete by]),
  [$3 - (10, 12) - 6$], [$(12, 6)$: 1 [10], $K_12$], [$(10, 6)$: 1], [mass], [(B.2)],
  [$3 - (9, 14) - 7$], [$(14, 7)$: 1 [29]], [$(9, 7)$: 1], [mass], [(B.2)],
  [$3 - (10, 12) - 8$], [$(12, 8)$: 1 [6]], [$(10, 8)$: 1], [mass], [(B.2)],
  [$3 - (8, 16) - 8$], [$(16, 8)$: 6 [163]], [$(8, 8)$: 1], [mass], [rootless],
  [$3 - (12, 8) - 8$], [$(8, 8)$: 1, $sqrt(3) E_8$], [$(12, 8)$: 51], [$s = f$], [(B.2)],
  [$3 - (9, 14) - 9$], [$(14, 9)$: 2 [21]], [$(9, 9)$: 1], [mass], [rootless],
  [$3 - (10, 12) - 10$], [$(12, 10)$: 1 [3]], [$(10, 10)$: 2], [mass], [rootless],
  [$3 - (13, 6) - 5$], [$(6, 5)$: 1 [1], $sqrt(3) E_6^\#$], [$(13, 5)$: 105], [mass], [rootless],
  [$3 - (6, 20) - 6$], [$(20, 6)$: 24], [$(6, 6)$: 1, $K_12$], [rootless], [rootless],
  [$3 - (16, 0) - 0$], [—], [$(16, 0)$: 11], [—], [[Hen13]; rootless],
  [$5 - (4, 16) - 4$], [$(16, 4)$: 1], [$(4, 4)$: 1], [rootless], [rootless],
  [$5 - (5, 12) - 5$], [$(12, 5)$: 1 [41]], [$(5, 5)$: 1], [mass], [rootless],
  [$5 - (6, 8) - 4$], [$(8, 4)$: 1 [5]], [$(6, 4)$: 5], [mass], [rootless; (B.2)],
  [$5 - (6, 8) - 6$], [$(8, 6)$: 1 [2]], [$(6, 6)$: 4], [mass], [rootless],
  [$5 - (7, 4) - 3$], [$(4, 3)$: 1 [1], $sqrt(5) A_4^\#$], [$(7, 3)$: 23], [mass], [rootless; (B.2)],
  [$5 - (8, 0) - 0$], [—], [$(8, 0)$: 9], [—], [rootless],
  [$7 - (4, 8) - 4$], [$(8, 4)$: 1 [8]], [$(4, 4)$: 3], [mass], [rootless],
  [$17 - (2, 0) - 0$], [—], [$(2, 0)$: 6 Hermitian classes; 4 trace lattices], [—], [[Jür15]; rootless],
  [$31 - (1, 2) - 1$], [$(2, 1)$: 1], [$(1, 1)$: 8], [direct], [class number],
)

For $s = f$ the lattice $F$ is $sqrt(p) M$ with $M$ even unimodular, so $F = sqrt(3) E_8$ for $3 - (12, 8) - 8$.

= Appendix C. The rootless mass over $ZZ[i]$ <appC>

#policy-box[
*What this appendix does.* It proves Theorem 4.5. §C.1 sets up the problem and states the strategy (King's inclusion–exclusion over root systems). §C.2 classifies the root lattices over $ZZ[i]$ and §C.3 their root sublattices, which gives the coefficients of the inclusion–exclusion. §C.4 proves a Siegel-type formula for the average number of embeddings of a lattice $T$, and §C.5 computes its local factors, the essential new point being the ramified prime 2. §C.6 explains how the 2-adic sums were computed and §C.7 lists the checks. The appendix can be read without the rest of the paper except §2.4.
]

== C.1 Setting and strategy

Let $cal(O) = ZZ[i]$, $pi.alt = 1 + i$, and $chi = chi_(-4)$ the character of $QQ(i)$. A _theta lattice_ of rank $r$ is a free $cal(O)$-module $L$ of rank $r$ with a positive definite Hermitian form $h$ (linear in the first variable) such that $h(L, L) subset cal(O)$, $h(x, x) in 2 ZZ$ for all $x$, and $L = {x in L times.o QQ : h(x, L) subset cal(O)}$ [HKN12]. With $(x, y) = "Re" h(x, y)$, theta lattices of rank $r$ are exactly the even unimodular lattices of rank $2 r$ with an automorphism $g$ ($x |-> i x$) satisfying $g^2 = -1$; this is the normalization $H' = d H$, $d = 2$, of the proof of Lemma 2.2. A _root_ is a vector with $h(x, x) = 2$. The mass $mu_r = sum_L 1 slash |"U"(L)|$ over isometry classes is known [HKN12, Thm 2]; in the notation of (C.6) below, $mu_r = S'_r$ for $r equiv 0 med (mod 4)$, the case of [HKN12, Thm 2] and the only one we use. We want the rootless mass $mu_r^emptyset$; then $M_2 = mu_16^emptyset$.

In the equivalent Möbius form of [Kin03, Prop. 1], for a lattice $T$ generated by roots let $r(T, L)$ be the number of embeddings $T -> L$ (isometric $cal(O)$-linear maps) and
#tagged("C.1", $ A_r (T) = sum_L r(T, L) / (|"U"(L)|). $)
§C.4–§C.5 compute $A_r (T)$ from local data without knowing the classes $L$, and §C.3 expresses $mu_r^emptyset$ as a finite linear combination of the $A_r (T)$.

== C.2 Root lattices over $ZZ[i]$

The roots of a theta lattice form a root system in the trace form, preserved by $g$. Since $g^2 = -1$ preserves every irreducible component, $g$ acts on the set of components as an involution.

#thm[
*Proposition C.1.* The irreducible root lattices over $cal(O)$ are the following, where $C[i] = C times.o_ZZ cal(O)$ and $G_m = {z in cal(O)^m : sum z_k in pi.alt cal(O)}$ with the standard form.
#align(center, table(columns: 4, align: (left, left, right, left),
  table.header([name], [root system over $ZZ$], [rank over $cal(O)$], [$|"U"|$]),
  [$A_n [i]$, $n >= 1$], [$A_n^2$], [$n$], [$4$ for $n = 1$, else $4 (n + 1)!$],
  [$D_n [i]$, $n >= 4$], [$D_n^2$], [$n$], [$"2,304"$ for $n = 4$, else $2^(n + 1) n!$],
  [$E_n [i]$, $n = 6, 7, 8$], [$E_n^2$], [$n$], [$2 |"Aut" E_n|$],
  [$G_m$, $m >= 2$], [$D_(2 m)$], [$m$], [$96$ for $m = 2$, else $4^m m!$],
  [$E_8^g$], [$E_8$], [4], [$"46,080"$],
))
]
#proofbox[
If $g$ interchanges two components $C$ and $C g$, identifying them by $g$ gives $C[i]$. If $g$ preserves a component $C$, then $g in "Aut" C$ with $g^2 = -1$. For $A_n$ ($n >= 2$) and $E_6$, $"Aut" C = W times {plus.minus 1}$ with $-1 in.not W$, so $g^2 = w^2 = -1$ is impossible; $A_1$, $E_7$ and $D_n$ with $n$ odd have odd rank. For $D_(2 m)$, $m >= 3$, $"Aut" C$ is the group of signed permutations, and $g^2 = -1$ forces a fixed-point-free involutive permutation with each transposition carrying sign product $-1$; these form one conjugacy class, giving $G_m$. For $D_4$ the 12 elements with $g^2 = -1$ in $W(F_4)$ form one class ($|W(F_4)| slash 96 = 12$), giving $G_2$. For $E_8$, the elements of $W(E_8)$ with characteristic polynomial $Phi_4^4$ form one conjugacy class, of size $"15,120"$ and with centralizer of order $"46,080"$ #calc("C11", "C20")\; this gives $E_8^g$. The orders $|"U"|$ in the table are centralizers of $g$. If $g$ interchanges $C$ and $C g$, the centralizer of $g$ in $"Aut"(C perp C g) = "Aut" C wreath C_2$ consists of the elements acting by the same $a in "Aut" C$ on both components, together with their products with $g$; hence $|"U"(C[i])| = 2 |"Aut" C|$, which gives the first three rows ($|"Aut" A_1| = 2$, $|"Aut" A_n| = 2 (n+1)!$ for $n >= 2$, $|"Aut" D_4| = "1,152"$, $|"Aut" D_n| = 2^n n!$ for $n >= 5$). For $G_m$, $m >= 3$, the centralizer of $g$ in the group of signed permutations is $(ZZ slash 4)^m times.r S_m$, of order $4^m m!$; the values for $G_2$ and $E_8^g$ are the centralizer orders above. #qed
]

Different ways of letting $g$ act on the same root system over $ZZ$ are different types over $cal(O)$; for instance $E_8^3$ gives the types $(E_8^g)^3$ and $E_8 [i] perp E_8^g$. The numbers of multisets of irreducible root lattices over $cal(O)$ of total rank exactly $0, 4, 8, 12, 16$ are $1, 12, 113, 802, "4,663"$; over all total ranks at most 16 there are $"12,934"$.

== C.3 Inclusion–exclusion

#thm[
*Proposition C.2.* For an irreducible root lattice $Y$ define $c(Y) = -sum_(Q subset.neq Y) c(Q)$, the sum over the sublattices $Q$ of $Y$ generated by roots (as $cal(O)$-modules), with $c(0) = 1$, and extend $c$ multiplicatively to orthogonal sums. Then
#tagged("C.2", $ mu_r^emptyset = sum_T c(T) (A_r (T)) / (|"U"(T)|), $)
the sum over isomorphism classes of root lattices $T$ over $cal(O)$ of rank at most $r$.
]
#proofbox[
Let $R$ be the root lattice of $L$. Every sublattice $Q subset L$ generated by roots lies in $R$, and is the orthogonal sum of its intersections with the irreducible components $Y_k$ of $R$, since every root lies in one component: each nonzero component in an orthogonal sum of even positive definite lattices has norm at least 2. Hence $sum_(Q subset R) c(Q) = product_k sum_(Q_k subset Y_k) c(Q_k)$, which is 1 if $R = 0$ and 0 otherwise by the definition of $c$. The number of sublattices of $L$ isomorphic to $T$ and generated by roots is $r(T, L) slash |"U"(T)|$. Weight by $1 slash |"U"(L)|$ and sum over $L$. #qed
]

To compute $c$ we need all root sublattices of each irreducible $Y$. For $0<=kappa<=2$ let $cal(O)_kappa = ZZ[zeta_(2^kappa)]$ ($cal(O)_0 = cal(O)_1 = ZZ$), $mu_(2^i)$ the group of $2^i$-th roots of unity, and $pi.alt_i = 1 - zeta_(2^i)$ for $i >= 2$, $pi.alt_1 = 2$, $pi.alt_0 = 0$. For $kappa >= 1$ put $G_k^((kappa)) = {z in cal(O)_kappa^k : sum z_p in pi.alt_kappa cal(O)_kappa}$; with the standard form, $G_k^((kappa))$ is $D_(2^(kappa - 1) k)$, $G_k^((1)) = D_k$ and $G_k^((2)) = G_k$; put $G_k^((0)) = {z in ZZ^k : sum z_p = 0} = A_(k - 1)$. Up to roots of unity, the roots of $G_k^((kappa))$ are the _edges_ $e_p - epsilon e_q$ ($p != q$, $epsilon in mu_(2^kappa)$) and the _half-edges_ $(1 - eta) e_p$ ($eta in mu_(2^kappa)$, $eta != plus.minus 1$); for $kappa = 0$ only the edges with $epsilon = 1$.

#thm[
*Proposition C.3.* Let $0<=kappa<=2$ and let $Q$ be an $cal(O)_kappa$-submodule of $G_k^((kappa))$ generated by roots. There are a partition of the $k$ coordinates into blocks $V$ and, for each block, roots of unity $u_p$ ($p in V$) and $0 <= i <= kappa$, such that $Q$ is the orthogonal sum of
$ Q_V = {sum_(p in V) z_p u_p e_p : z_p in cal(O)_kappa, sum_p z_p in pi.alt_i cal(O)_kappa}, $
where for $|V| = 1$ either $Q_V = 0$ or $i >= 2$. The pair $(V, i)$ is determined by $Q$, and for fixed $(V, i)$ there are $(2^(kappa - i))^(|V| - 1)$ such $Q_V$. The type of $Q_V$ is the base change to $cal(O)_kappa$ of $A_(|V| - 1)$ for $i = 0$ and of $D_(|V|)$ for $i = 1$ (for $kappa = 2$ these are $A_(|V| - 1)[i]$ and $D_(|V|)[i]$), and $G_(|V|)^((i)) times.o cal(O)_kappa$ for $i >= 2$.
]
#proofbox[
Let $V$ be a connected component of the graph on the coordinates whose edges are the edges of $Q$. Choose a spanning tree and rescale the basis vectors by roots of unity ($e_q |-> u_q e_q$, an isometry of $G_k^((kappa))$) so that the tree edges become $e_p - e_q$; then $e_p - e_q in Q$ for $p, q in V$. Let $Gamma = {eta in mu_(2^kappa) : (1 - eta) e_p in Q}$. Since $(1 - eta) e_q = (1 - eta) e_p - (1 - eta)(e_p - e_q)$, $Gamma$ does not depend on $p in V$, and $1 - eta eta' = (1 - eta) + eta (1 - eta')$ shows that it is a subgroup, $Gamma = mu_(2^i)$. An edge $e_p - epsilon e_q = (e_p - e_q) + (1 - epsilon) e_q$ lies in $Q$ only if $epsilon in Gamma$. So $Q_V$ is generated by the $e_p - e_q$ and by $(1 - eta_i) e_p$ for a generator $eta_i$ of $Gamma$; since $(1 - eta_i) cal(O)_kappa = pi.alt_i cal(O)_kappa$, $Q_V$ has the stated form (for $i = 1$, $2 e_p$ is not a root but lies in the span of $e_p plus.minus e_q$). Conversely, for $eta$ of order $2^j$ the valuation of $1 - eta$ at the prime of $cal(O)_kappa$ above 2 is $2^(kappa - j)$, so $1 - eta in pi.alt_i cal(O)_kappa$ if and only if $eta in mu_(2^i)$; hence the roots of $Q_V$ are those with $eta, epsilon in mu_(2^i)$, $i$ is determined by $Q_V$, and no new roots arise. The count is the number of classes of the ratios $u_q slash u_p$ along the tree modulo $mu_(2^i)$. #qed
]

For $1<=kappa<=2$, if $C$ is a lattice over $cal(O)_b$, $0<=b<=kappa$, the roots of $C times.o_(cal(O)_b) cal(O)_kappa$ split into $2^(kappa-max(b,1))$ orthogonal copies of the roots of $C$, permuted cyclically by $g$; so the root sublattices of $C times.o cal(O)_kappa$ correspond bijectively to those of $C$. In particular those of $A_k [i]$ and $D_k [i]$ correspond to the root sublattices of $A_k$ and $D_k$ over $ZZ$, which are given by Proposition C.3 with $kappa = 0, 1$. For $kappa = 2$: a root sublattice of $G_m$ is given by a partition of the coordinates and, on each block $V$, one of $A_(|V| - 1) [i]$ ($4^(|V| - 1)$ ways), $D_(|V|) [i]$ ($2^(|V| - 1)$ ways) or $G_(|V|)$ (one way), singletons being unused or $G_1 = A_1 [i]$. If $w(n)$ denotes the total coefficient of a block of $n$ coordinates (for instance $w(n) = 2^(n-1) c(A_(n-1)) + c(D_n)$ for $kappa = 1$), Proposition C.3 says that the sum of $product_(V in cal(P)) w(|V|)$ over the set partitions $cal(P)$ of $k$ coordinates is the sum of $c(Q)$ over all root sublattices $Q$ of $G_k^((kappa))$, which is $0$ for $k >= 2$ (and for every $k$ in the $G$ family) and $1$ for $k = 1$ otherwise; so $exp sum_n w(n) x^n slash n! = 1 + x$ for the $A, D$ families and $= 1$ for the $G$ family, that is $w(n) = (-1)^(n-1) (n-1)!$ or $0$. Solving for the last term of $w(n)$ gives
$ c(A_n)=(-1)^n n!, quad c(D_n)=(-1)^n (n-1)!(2^(n-1)-1), quad c(G_m)=(-1)^m 2^(m-1) (m-1)!. $
The $D$ formula includes $D_2=A_1^2$, $D_3=A_3$; the $G$ formula includes $G_1=A_1[i]$. Base change gives $c(C[i])=c_(ZZ)(C)$. For $E_n [i]$ and $E_8^g$ the root sublattices were counted from explicit roots ($c(E_8 [i]) = -"28,183,680"$, $c(E_8^g) = -"5,184"$). For ADE systems a symmetric root subset closed under root sums is exactly the set of roots in its integer span, so it represents a root-generated submodule after imposing $g$-stability. The closed formulas agree with explicit enumeration of these subsets for small $(kappa,k)$ #calc("C11").

*The exceptional coefficients.* For $E_6[i],E_7[i],E_8[i]$, enumerate all ADE root types of rank at most that of the real root system. Count the ordered root tuples with each type's positive definite simple-root Gram and divide by the full automorphism-group order of the source root lattice. This counts each root sublattice once. At each recursive step the program uses reflections in roots orthogonal to the already selected tuple, splits the candidates into their orbits, and multiplies one continuation count by the orbit size. These reflections fix the selected tuple pointwise, so they give bijections between all continuations in an orbit; the subgroup need not be the full stabilizer.

For $E_8^g$, take the 240 explicit real roots and the rotation $g$ of successive coordinate pairs. Test both real and imaginary parts of the Hermitian Gram. The Hermitian reflection in a root $x$ is the product of the real reflections in the orthogonal roots $x,g x$; it preserves the root set, commutes with $g$, and fixes the preceding Hermitian-orthogonal tuple. Enumerating the possible Gaussian root types of rank at most 4 gives the same orbit-counting argument. Omitting the ambient exceptional lattice itself and taking the negative coefficient sum gives $c(E_6[i])="7,680"$, $c(E_7[i])=-"233,280"$, $c(E_8[i])=-"28,183,680"$, and $c(E_8^g)=-"5,184"$ #calc("C11", "C19").

== C.4 A Siegel-type formula

Let $Lambda_n$ be the lattice of even Hermitian $n times n$ matrices over $cal(O)$ (diagonal in $2 ZZ$). For a theta lattice with Gram matrix $S$ over $cal(O)_p = cal(O) times.o ZZ_p$ and $T in Lambda_n$ positive definite, the local density is
$ alpha_p (T) = lim_(a -> oo) p^(-a (2 r n - n^2)) \#{X in M_(r times n)(cal(O) slash p^a) : X^* S X equiv T med (mod p^a Lambda_n)}. $
All theta lattices of rank $r$ lie in one genus (by Lemma C.9(iii) at $p = 2$ and the classical theory at odd $p$), so $alpha_p (T)$ does not depend on $L$.

#thm[
*Theorem C.4.* For $n < r$,
#tagged("C.3", $ A_r (T) = mu_r dot 2^n (pi^(n (2 r - n + 1) slash 2)) / (product_(j = 0)^(n - 1) Gamma(r - j)) det(T)^(r - n) product_p alpha_p (T), $)
and the product converges absolutely.
]
_Strategy of the proof._ This is Siegel's formula in the form "an average number of representations is a product of local densities", proved by unfolding an adelic integral. By Witt's theorem the unitary group acts transitively on the solutions of $X^* S X = T$, so the sum over embeddings becomes an integral over the solutions, which factors into local volumes. The Tamagawa numbers of the unitary group and of the stabilizer are both 2 and cancel.
#proofbox[
Let $V = L times.o QQ$, $G = "U"(V)$ and $Omega = {X in V^n : X^* S X = T}$. $Omega(QQ)$ is non-empty: with $delta = det S det T > 0$, the form $T perp chevron.l delta, 1, dots, 1 chevron.r$ has the same dimension as $V$, is positive definite, and has determinant $det S (det T)^2$, which lies in the class of $det S$ modulo norms since $det T in QQ^times$; so it is isometric to $V$, because a Hermitian space over a number field is determined by its dimension, its indices at the real places and its determinant modulo norms [Shi08, Thm 2.2(i)] (Landherr). By Witt's theorem $G$ acts transitively on $Omega$, over $QQ$ and over every completion, with stabilizer $H = "U"(W)$, $W = T^perp$, of rank $r - n >= 1$. Let $K_f subset G(AA_f)$ be the stabilizer of $L times.o hat(ZZ)$ and $phi$ the characteristic function of $(L times.o hat(ZZ))^n$. Unfolding,
$ integral_(G(QQ) backslash G(AA)) sum_(x in Omega(QQ)) phi(g^(-1) x) dif g = "vol"(H(QQ) backslash H(AA)) integral_(Omega(AA)) phi dif mu, $
which converges since $G(RR)$ is compact. The left side is $"vol"(K_f G(RR)) A_r (T)$, and for $T = emptyset$ it gives $"vol"(K_f G(RR)) mu_r = "vol"(G(QQ) backslash G(AA))$. With Tamagawa measures $omega_G = omega_Omega and omega_H$, the convergence factors $L_p (1, chi)$ of $G$ and $H$ cancel, and $tau(H) = tau(G)$ [Wei65, n°53–54] ($tau("U"(m)) = 2$ for all $m >= 1$, including the norm-one torus $m = 1$). The local volumes of $omega_Omega = dif X slash dif T$, normalized so that $cal(O)_p^(r n)$ and $Lambda_n times.o ZZ_p$ have volume 1, are the $alpha_p (T)$, and at infinity $alpha_infinity (T) = 2^n c det(T)^(r - n)$ with $c$ from Lemma C.5: the covolume of $cal(O)^(r n)$ in $CC^(r n)$ is 1 since $det S = 1$, and that of $Lambda_n$ is $2^n$. #qed
]

#thm[
*Lemma C.5.* Let $dif X$ be Lebesgue measure on $M_(r times n)(CC) = RR^(2 r n)$ and $dif T$ Lebesgue measure on Hermitian matrices in the coordinates $t_(k k)$, $"Re" t_(k l)$, $"Im" t_(k l)$ ($k < l$). Under $X |-> T = X^* X$ ($n <= r$), $dif X$ pushes forward to $c det(T)^(r - n) dif T$ on positive definite $T$, with $c = pi^(r n - n (n - 1) slash 2) slash product_(j = 0)^(n - 1) Gamma(r - j)$.
]
#proofbox[
For $g in "GL"_n (CC)$, $X |-> X g$ multiplies $dif X$ by $|det g|^(2 r)$ and $T |-> g^* T g$ multiplies $dif T$ by $|det g|^(2 n)$; since $"GL"_n (CC)$ acts transitively on positive definite matrices, the density has the form $c det(T)^(r - n)$. To find $c$, compare $pi^(r n) = integral e^(-"tr" X^* X) dif X = c integral_(T > 0) e^(-"tr" T) det(T)^(r - n) dif T$. Write $T = B^* B$ with $B$ upper triangular with positive diagonal; then $dif T = 2^n product_k b_(k k)^(2 (n - k) + 1) dif B$ and $det T = product_k b_(k k)^2$, and the integral is $pi^(n (n - 1) slash 2) product_(k = 1)^n Gamma(r - k + 1)$. #qed
]

For $n = r$ we do not use (C.3). In that case $r(T, L) slash |"U"(L)|$ summed over $L$ counts the even unimodular overlattices $T' supset T$ in $T times.o QQ$, because an embedding of a lattice of full rank is determined by its image and the only automorphism of $T'$ fixing $T$ pointwise is the identity. Overlattices correspond to submodules $C subset T^* slash T$ with $C^perp = C$ that are isotropic for $q(x) = h(x, x) mod 2$, and $T^* slash T$ is the direct sum of its $p$-parts; hence
#tagged("C.4", $ A_r (T) = product_p \#{C subset D_p (T) : C "isotropic", C^perp = C} quad (n = r), $)
with $D_p (T) = T^* slash T times.o ZZ_p$.

== C.5 Local densities

An embedding $X$ is _primitive_ at $p$ if $X$ has rank $n$ modulo every prime of $cal(O)$ above $p$.

#thm[
*Lemma C.6.* If $X$ is primitive, the differential of $X |-> X^* S X$ maps onto $Lambda_n times.o ZZ_p$. Consequently the primitive part $alpha_(p, "prim")(T)$ of $alpha_p (T)$ is determined by solutions modulo $p$ for odd $p$ and modulo 2 for $p = 2$.
]
#proofbox[
The differential is $A |-> X^* S A + (X^* S A)^*$. $X^* S$ maps $cal(O)_p^r$ onto $cal(O)_p^n$, and $Y |-> Y + Y^*$ maps $M_n (cal(O)_p)$ onto $Lambda_n times.o ZZ_p$, because $"Tr"_(cal(O) slash ZZ) (cal(O)) = 2 ZZ$ covers the even diagonal. Replacing $X$ by $X + 2 A$ changes $X^* S X$ by $2 dif F_X (A) + 4 A^* S A$, and $4 A^* S A in 4 Lambda_n$ since $L$ is even; so each solution modulo $2 Lambda_n$ lifts to exactly $2^(2 r n - n^2)$ solutions at every higher level. The odd case is the usual Hensel argument. #qed
]

On $overline(V) = L slash pi.alt L tilde.equiv FF_2^r$ the form $omega(x, y) = h(x, y) mod pi.alt$ is non-degenerate and alternating (proof of Lemma 2.2). For $T$ let $2 a$ be the rank of $omega$ on $T slash pi.alt T$ and $d = n - 2 a$ the dimension of its radical; $d$ is the number of generators of the $cal(O)$-module $D_2 (T)$.

#thm[
*Proposition C.7.* $alpha_(2, "prim")(T) = product_(k = 0)^((n + d) slash 2 - 1) (1 - 4^k 2^(-r))$. For odd $p$, with $d$ the number of generators of $D_p (T)$, $alpha_(p, "prim")(T) = product_(j = r - n - d + 1)^r (1 - chi(p)^j p^(-j))$.
]
#proofbox[
Write the columns of $X$ modulo 2 as $x_j = hat(w)_j + pi.alt u_j$ with $w_j in overline(V)$ independent (primitivity) and $u_j in overline(V)$. Since $overline(pi.alt) equiv pi.alt$ modulo 2,
$ h(x_j, x_k) equiv h(hat(w)_j, hat(w)_k) + pi.alt (omega(w_j, u_k) + omega(u_j, w_k)), quad h(x_j, x_j) slash 2 equiv h(hat(w)_j, hat(w)_j) slash 2 + omega(w_j, u_j) quad (mod 2). $
For fixed $w$ these are linear conditions on $u$: the map $u |-> B = (omega(w_k, u_j))_(k, j) in M_n (FF_2)$ is onto, and the conditions prescribe $B_(j j)$ and $B_(j k) + B_(k j)$, which are $n (n + 1) slash 2$ independent linear forms. So there are $2^(r n - n (n + 1) slash 2)$ choices of $u$ for each admissible $w$. The admissible $w$ are the isometric injections of $(T slash pi.alt T, omega)$ into $(overline(V), omega)$, whose number is $product_(k = 0)^(a - 1) (2^(r - 2 k) - 1) 2^(r - 2 k - 1) dot product_(j = 0)^(d - 1) (2^(r - 2 a - j) - 2^j)$. Multiplying and dividing by $2^(2 r n - n^2)$ gives the formula. For odd $p$ the same count over $cal(O) slash p$ (a Hermitian space over $FF_(p^2)$ if $p equiv 3 med (mod 4)$, a pair of dual spaces over $FF_p$ if $p equiv 1$) gives the second formula. #qed
]

#thm[
*Theorem C.8.* Let $t = r - n$. For each isotropic $cal(O)$-submodule $C subset D_p (T)$ let $T_C supset T$ be the corresponding even overlattice and $d(C)$ the number of generators of $C^perp slash C$. Then
#tagged("C.5", $ alpha_p (T) = sum_C |C|^(-t) alpha_(p, "prim")(T_C), $)
and, with $S_t = 2^(1 - t) product_(j = 1)^t |L(1 - j, chi^j)|$ ($chi^j = chi$ for odd $j$, the trivial character for even $j$, so that $L(1 - j, chi^j) = zeta(1 - j)$), $S'_t = 2^(-t [t "odd"]) S_t$ and $S'_0 = 1$,
#tagged("C.6", $ A_r (T) = S'_t det(T)^t product_p f_p (T, t), quad f_p (T, t) = sum_C |C|^(-t) product_(j in J_p (t, d(C))) (1 - chi(p)^j p^(-j)), $)
where $J_p (t, d) = {t - d + 1, dots, t}$ for odd $p$, $J_2 (t, d) = {j "even" : t - d + 2 <= j <= t}$, and $chi(2)^j$ is read as 1.
]
_Strategy of the proof._ An embedding of $T$ need not be primitive: the closure of its image is an overlattice $T_C supset T$, given by an isotropic submodule $C$. (C.5) splits the local density according to this overlattice, and Proposition C.7 gives each primitive density in closed form. For (C.6), the factors that do not depend on $T$ are divided out and, with the functional equations of the $L$-functions, combined into the constant $S'_t$.
#proofbox[
(C.5): the closure of the image of $X$ is some $T_C$; the embeddings with closure $T_C$ are the primitive embeddings of $T_C$ composed with the inclusion, and the change of basis multiplies the volume on the $X$ side by $|C|^(-r)$ and on the $Lambda_n$ side by $|C|^(-n)$.

(C.6) for $n < r$: if $p divides.not det T$ then $D_p (T) = 0$ and $alpha_p (T) = product_(j = t + 1)^r (1 - chi(p)^j p^(-j))$ by Proposition C.7 (for $p = 2$, the factors with $j = r - 2 k$ even). Dividing $alpha_p (T)$ by this product gives $f_p (T, t)$ for every $p$, again by Proposition C.7. The product of the divisors over all $p$ is $product_(j = t + 1)^r L(j, chi^j)^(-1)$. By the functional equations, $|L(1 - j, chi^j)| = 2 (j - 1)! dot 2^((2 j - 1) [j "odd"]) (2 pi)^(-j) L(j, chi^j)$. Inserting this and $mu_r = S'_r$ into (C.3), the factors $pi$ and $(j - 1)!$ cancel against those of (C.3), and the remaining power of 2 is $2^E$ with $E = sum_(j = t + 1)^r (1 - j + (2 j - 1)[j "odd"]) - r [r "odd"] + t [t "odd"]$; the summand is $j$ for odd $j$ and $1 - j$ for even $j$, so consecutive pairs (odd $j$, $j + 1$) cancel, and the parity terms cancel the unpaired end terms. Hence $E = 0$ and (C.6) holds.

(C.6) for $n = r$: $t = 0$, $S'_0 = 1$, and in $f_p (T, 0)$ the product is 0 unless $d(C) = 0$ (it contains the factor $j = 0$, and $d(C)$ is even at $p = 2$), so $f_p (T, 0)$ counts the $C$ with $C^perp = C$, which is (C.4). #qed
]

The formula expresses the averaged embedding number through the complement rank $t=r-n$ and local gluing data. It need not be the ordinary mass of a single complement genus.

== C.6 Computing the 2-adic sums

For odd $p$, $D_p (T)$ is small and the isotropic submodules are counted in closed form. At $p = 2$, $D_2 (T)$ can be large ($D_2 (A_1 [i]^12)$ has about $2.2 times 10^14$ isotropic submodules), but only the distribution of $(|C|, d(C))$ is needed. We compute it by adding the irreducible components of $T$ one at a time. The state is the isometry type of the discriminant module of an overlattice of the part already added, recorded by its Jordan decomposition:

#thm[
*Lemma C.9.* Let $L$ be an even lattice over $cal(O)_2 = ZZ_2 [i]$.
+ $L$ is an orthogonal sum of Jordan components of rank 1 and 2.
+ A component $chevron.l 2^e u chevron.r$ of rank 1 is determined by $u mod 4$.
+ A plane of scale $pi.alt^(2 e)$ is $2^e H$, $H$ the hyperbolic plane; in particular an even unimodular lattice is $H^(r slash 2)$.
+ A plane of scale $pi.alt^(2 e + 1)$ is one of the two planes $P(e, t)$, $t = 0, 1$, with Gram matrix $mat(2^(e + 1) t, pi.alt^(2 e + 1); overline(pi.alt)^(2 e + 1), 2^(e + 1) t)$. $P(e, 0)$ has a primitive isotropic vector and $P(e, 1)$ has none, so they are not isometric.
+ $chevron.l 2^e chevron.r^2 tilde.equiv chevron.l 3 dot 2^e chevron.r^2$ and $P(e, 1)^2 tilde.equiv P(e, 0)^2$.
]
#proofbox[
We use that the norms of units are $N(cal(O)_2^times) = 1 + 4 ZZ_2$: the inclusion holds modulo 4 since $u^2 + v^2 equiv 1$ for $u + v i$ a unit, and $5 = N(1 + 2 i)$ generates $1 + 4 ZZ_2$ topologically while the norm group is closed. Let $s$ be the $pi.alt$-adic scale of $L$.

(i) If a diagonal entry attains the scale, it splits off a component of rank 1; otherwise only off-diagonal entries attain it, and a plane spanned by two basis vectors splits off.

(ii) $chevron.l 2^e u chevron.r tilde.equiv chevron.l 2^e u' chevron.r$ if and only if $u' slash u$ is a norm of a unit, that is $u' equiv u$ modulo 4.

(iii) Let the plane be $2^e mat(2 a, epsilon; overline(epsilon), 2 b)$ with $epsilon$ a unit; rescaling the basis we may take $epsilon = 1$. For $x = e_1 + lambda e_2$, $lambda = u + v i$, half the norm divided by $2^e$ is $a + b (u^2 + v^2) + u$, whose derivative in $u$ is odd, so Hensel's lemma applies; a root modulo 2 always exists ($(0, 0)$ if $a$ is even, $(1, 1)$ if $a$ and $b$ are odd, $(1, 0)$ if $b$ is even). So there is a primitive isotropic $x$ with $h(x, e_2)$ a unit, and $e_2$ can be moved by a multiple of $x$ to an isotropic vector, since $"Tr"(cal(O)) = 2 ZZ$. Hence the plane is $2^e H$.

(iv) Let the plane be $mat(2^(e + 1) a, pi.alt^s epsilon; overline(pi.alt^s epsilon), 2^(e + 1) b)$ with $s = 2 e + 1$. For $x = e_1 + lambda e_2$ the norm is $2^(e + 1)(a + b (u^2 + v^2) + c u + d v)$ with $c, d$ odd. If $a$ is even, Hensel's lemma at $(0, 0)$ gives an isotropic $x$ and, as in (iii), the plane is $P(e, 0)$ (the plane with an isotropic primitive vector); if $b$ is even, exchange $e_1$ and $e_2$. If $a$ and $b$ are odd, Hensel's lemma from the root $(0, 0)$ of $a - 1$ modulo 2 gives $x$ of norm $2^(e + 1)$; rescale $e_2$ so that $h(x, e_2) = pi.alt^s$. With $rho = 2^(e + 1) pi.alt^(-s)$ and $y = eta e_2 + mu x$, $eta = 1 - mu overline(rho)$, we keep $h(x, y) = pi.alt^s$, and the norm of $y$ divided by $2^(e + 1)$ is a polynomial in $mu$ that is odd at $mu = 0$ with a derivative having odd coefficients; Hensel's lemma gives $y$ of norm $2^(e + 1)$. This plane is $P(e, 1)$.

(v) In $chevron.l 2^e chevron.r^2$ the vector $(1, pi.alt)$ has norm $3 dot 2^e$, equal to the scale, so it splits off and its complement is determined by the determinant. In $P(e, 1)^2$ with basis $e_1, f_1, e_2, f_2$, the vectors $e_1 + e_2$ and $f_1$ span a plane with an isotropic primitive vector, hence $P(e, 0)$; with $rho$ as above, $z = -rho e_1 + f_1 + rho e_2 - f_2$ is a primitive vector of its orthogonal complement of norm $2 dot 2^(e + 1)$, so the complement is also $P(e, 0)$. For $P(e, 1)$, $h(x, x) slash 2^(e + 1) equiv 1 med (mod 2)$ for every primitive $x$, by the computation above with $a$, $b$ odd. #qed
]

#include "m2-odd-factors.typ"

#ptable(
  columns: 4, align: (left, left, left, left),
  caption: [Table C1. The 2-adic Jordan blocks of the finite transition. $u in {1, 3}$ modulo 4; the discriminant module is that of the block, with its generator(s).],
  table.header([block], [Gram matrix], [discriminant module], [occurs in the root lattices]),
  [$chevron.l 2^e u chevron.r$, $e >= 1$], [$(2^e u)$], [$cal(O) slash 2^e$, generator $e_1 slash 2^e$, $q = u slash 2^e$], [$A_n [i]$ ($n + 1 = 2^e m$, $m$ odd, $u equiv m n$), $D_n [i]$ ($n$ odd: $e = 2$, $u equiv n$; $n equiv 2$ mod 4: two blocks with $e = 1$, $u equiv n slash 2$), $E_7 [i]$ ($e = 1$, $u = 3$), $G_m$ ($m$ odd, $e = 1$, $u equiv m$)],
  [$2^e H$, $e >= 1$], [$2^e mat(0, 1; 1, 0)$], [$(cal(O) slash 2^e)^2$], [$D_n [i]$, $n equiv 0$ mod 4 ($e = 1$)],
  [$P(e, t)$, $e >= 0$, $t = 0, 1$], [$mat(2^(e+1) t, pi.alt^(2e+1); overline(pi.alt)^(2e+1), 2^(e+1) t)$], [$(cal(O) slash pi.alt^(2 e + 1))^2$], [$G_m$, $m$ even ($e = 0$; $t = 1$ iff $m equiv 2$ mod 4)],
)
The state of the programme is a pair $(ell, E)$: $ell = log_2 |C|$, and $E$ is the multiset of non-unimodular Jordan blocks of the current overlattice, normalized by Lemma C.9(v) ($P(e, 1)^2 -> P(e, 0)^2$, $chevron.l 3 dot 2^e chevron.r^2 -> chevron.l 2^e chevron.r^2$); unimodular blocks are dropped. The blocks $P(0, t)$ are processed first, in closed form (below), then the planes $2^e H$, then the blocks of rank 1. A transition adds $c = log_2 |Z|$ to $ell$ and replaces $E$ by the normalized Jordan decomposition of $Gamma_phi^perp slash Gamma_phi$; at the end the weights are aggregated by $(ell, d(E))$, $d(E)$ the sum of the ranks of the blocks of $E$. $E_6 [i]$, $E_8 [i]$ and $E_8^g$ contribute no block at 2.
For $G_m$ the discriminant group of $D_(2 m)$ is ${0, v, s, c}$ with $q(v) = 1$, $q(s) = q(c) = m slash 2$; $g$ has $m$ sign changes, so for even $m$ it lies in $W(D_(2 m))$, fixes $s$ and $c$, and the discriminant module is $(cal(O) slash pi.alt)^2$ with $t = 0$ or $1$ according as $m equiv 0$ or $2$ modulo 4, while for odd $m$ it exchanges $s$ and $c$ and the module is $cal(O) slash 2$.

*The finite transition.* Write $cal(O)_2=ZZ_2[i]$ and $pi.alt=1+i$. A state after adding blocks $B_1,dots,B_a$ is $(ell,E)$, with integer weight $W_(a)(ell,E)$ equal to the number of isotropic submodules $C$ of their direct sum for which $|C|=2^ell$ and $C^perp slash C tilde.equiv E$. The initial state is $W_(0)(0,0)=1$.

#thm[
*Lemma C.10 (graph transition).* Let $A,B$ be finite discriminant modules. Isotropic submodules $C subset A perp B$ are in bijection with triples consisting of an isotropic $C_0 subset A$, a submodule $Z subset B$, and an $cal(O)_2$-linear map $phi:Z -> E=C_0^perp slash C_0$ whose graph $Gamma_phi subset E perp B$ is isotropic. The map need not be injective. Then $|C|=|C_0| dot |Z|$ and $C^perp slash C tilde.equiv Gamma_phi^perp slash Gamma_phi$.
]
#proofbox[
Take $C_0=C inter A$ and project $C$ onto $B$. Orthogonality to $C_0$ puts $C$ inside $C_0^perp perp B$. Modulo $C_0$ the projection onto its image $Z$ is injective, so it is the graph of a unique linear map. Conversely, the inverse image of any such isotropic graph is the unique corresponding $C$. The order and residual-form identities follow from the same quotient. #qed
]

For each actual submodule $Z$ of the new block, fix one generating set and count each allowed map once. There is no division by $|"Aut" Z|$. Put $c=log_2|Z|$ and $E'=Gamma_phi^perp slash Gamma_phi$; this pair adds $W_(a)(ell,E)$ to $W_(a+1)(ell+c,E')$. Isometric representatives of $E$ have bijective sets of these transitions. It is enough that each normalization preserves isometry; a complete canonical classification of state keys is not required.

*Conditions on maps.* A rank-1 block has $B tilde.equiv cal(O)_2 slash pi.alt^s$, and all its submodules are $pi.alt^j B$, $0<=j<=s$. With fixed generator $z=pi.alt^(j-s)$, its image $y$ must satisfy $pi.alt^(s-j)y=0$ and $q_(E)(y)+q_(B)(z)=0$. A rank-2 block has $B tilde.equiv (cal(O)_2 slash pi.alt^s)^2$. Its submodules are at most two-generated; enumerate cyclic submodules, then their pairwise sums, retaining each resulting subset once. For fixed generators $z_1,z_2$ the candidate images $y_1,y_2$ must satisfy every relation
$ a_1 z_1+a_2 z_2=0 quad arrow.r.double quad a_1 y_1+a_2 y_2=0, $
as well as $"ann"(z_i) y_i = 0$, $q_(E)(y_i)+q_(B)(z_i)=0$ in $QQ_2 slash 2 ZZ_2$, and $h_(E)(y_1,y_2)+h_(B)(z_1,z_2)=0$ in $QQ_(2)(i) slash cal(O)_2$. These conditions are necessary and sufficient for an isotropic graph. Relations with coefficients in $cal(O)_2 slash pi.alt^s$ then suffice, since $pi.alt^s y_i = 0$. Individual orders or norms alone do not replace the relations or the mixed Hermitian condition.

*Orbit compression.* In each Jordan block use a subgroup of its actual isometry group. For a single image use its orbits on elements; for two images use its diagonal orbits on ordered pairs, preserving their relative position. If a block occurs $m$ times, an orbit pattern with $k_j$ copies of an orbit of size $v_j$ represents
$ m! slash (product_j k_j!) product_j v_j^(k_j), quad sum_j k_j=m $
candidate images. Permuting equal blocks and applying their isometries takes any two candidates of this pattern to each other, preserving relations, isotropy and the residual form. Thus one representative has this multiplicity. A subgroup in place of the whole block isometry group refines the orbits and increases work; it does not omit candidates. The implementation retains joint ordered-pair orbits, not just the two marginal element orbits.

*New states and endpoints.* Add the graph generators to the local lattice and compute its Jordan decomposition. In a one-generator step, replace the basis vector in the coordinate of least $pi.alt$-adic valuation by the new vector. Every old basis vector is still integral over this new local basis; its index is the annihilator exponent of the added class. Repeating this gives the index $2^c$ and the residual form. Every transition has $|E'|=|E| dot |B| slash 2^(2c)$.
The $P(0,t)$ blocks are processed first. Their discriminant forms are the binary planes $q(x,y)=x y$ for $t=0$ and $q(x,y)=x+y+x y$ for $t=1$. For $m$ such planes with sign $epsilon=(-1)^(sum t)$, the number of totally singular $k$-subspaces is
$ (product_(j=0)^(k-1) (2^(m-j)-epsilon)(2^(m-j-1)+epsilon)) slash (product_(j=1)^k (2^j-1)). $
It is zero for $k>m$ and for $epsilon=-1$, $k=m$, and one for $k=0$. The residual space has dimension $2(m-k)$ and sign $epsilon$. Count ordered bases: given $V_j = chevron.l v_1, dots, v_j chevron.r$ totally singular, $v_(j+1)$ is a lift of a non-zero singular vector of $V_j^perp slash V_j$, a non-degenerate space of dimension $2(m - j)$ and sign $epsilon$ with $(2^(m-j) - epsilon)(2^(m-j-1) + epsilon)$ such vectors, each with $2^j$ lifts, all singular; dividing by $|GL_k (FF_2)| = 2^(k(k-1) slash 2) product_(j=1)^k (2^j - 1)$ cancels the factors $2^j$. These are the closed initial weights before the $H$ and rank-1 blocks are added.
Finally aggregate by $(ell,d(E))$, where $d(E)$ is the sum of the ranks of the remaining non-unimodular blocks, and substitute the counts into (C.6). In the present even ambient rank, $t$ and $d(E)$ have the same parity. If $d(E)>t$, the product contains its zero factor at $j=0$; for $t=0$ only $d(E)=0$ survives. This includes the full-rank, Lagrangian endpoint. The finite computation and its local basis operations are part of #calc("C11")\; the small checks support their implementation and do not replace this transition argument.

The full root-type index and all terms of (C.2) are also recorded in #calc("C19"). This reconstruction regenerates the exceptional root coefficients, checks the index against coefficient extraction from $product_Y (1-x^("rank" Y))^(-1)$, and sums in exact arithmetic. It reuses fixed local DP statistics for 1,498 block multisets; those transition counts are not independently recomputed by this aggregation.

== C.7 Checks and the value

For ranks $r = 4, 8, 12$ the rootless masses are known independently, and (C.2) reproduces them exactly (Table C2). At rank 12 the only rootless theta lattice is the Leech lattice with its unique class of automorphisms $g$ with $g^2 = -1$ in $"Co"_0$, whose centralizer has order $"2,012,774,400"$ (computed from the character table of $2."Co"_1$ in GAP); this value is obtained independently of (C.2). At rank 8 three root systems have non-zero mass, and $1 slash m(R)$ equals the order of the unitary group of the corresponding lattice for all three; at rank 12 the 25 root systems with non-zero mass map to exactly the 19 Niemeier root systems that admit a structure over $ZZ[i]$ [HKN12, Remark 4(d)].

#ptable(
  columns: 4, align: (right, right, right, left),
  caption: [Table C2. Rootless masses of theta lattices computed by (C.2) #calc("C11"), with the independent values.],
  table.header([$r$], [root types], [$mu_r^emptyset$ from (C.2)], [independent value]),
  [4], [22], [0], [0 (the only theta lattice is $E_8$ [HKN12, Cor. 1])],
  [8], [254], [0], [0 (no rootless even unimodular lattice of rank 16)],
  [12], [2,035], [$1 slash "2,012,774,400"$], [$1 slash |C_("Co"_0)(g)|$ (Leech lattice)],
  [16], [$"12,934"$], [$"4,128,748,922,231,417" slash "120,324,096,000"$], [Proposition 4.6],
)

At rank 16 the value is $mu_16^emptyset approx "34,313.567"$, with denominator $2^21 dot 3^3 dot 5^3 dot 17$. Of the $"12,934"$ root types $R$, the mass $m(R)$ of the theta lattices whose root system is exactly $R$, obtained from the values $A_16 (T)$ by solving $A_16 (T) = sum_R m(R) r(T, R)$, is non-zero for $"1,154"$ types, exactly 0 for the other $"11,780"$, and never negative; an error in a single value $A_r (T)$ would make some mass that is truly 0 non-zero, except for 240 of the $"12,933"$ non-empty types $T$. These are consistency checks; the value rests on the proofs above and is confirmed by Proposition 4.6.



= Appendix D. Rootless masses for elements of odd prime order <appD>

#policy-box[
*What this appendix does.* It proves the values of the rootless masses used in §4.2 and Appendix B.4. §D.1 shows that the Hermitian lattices $Z$ of a given type form one genus and recalls its mass. §D.2 derives the rootless masses of unimodular Hermitian lattices from (4.2). §D.3 computes the rootless masses of the genera $(z, s)$ over $ZZ[zeta_p]$, $p = 5, 7, 17$, by inclusion–exclusion over the root lattices $A_(p-1)$ and $E_8$, with the method of Appendix C. §D.4 treats the $ZZ$-lattices of the genera $(20, 6)$ for $p = 3$ and $(16, 4)$ for $p = 5$, and the Eisenstein lattices of the genus $(13, 5)$.
]

Throughout, $p$ is an odd prime, $KK = QQ(zeta_p)$, $KK^+$ its real subfield, $cal(O) = ZZ[zeta_p]$, $cal(O)^+ = cal(O) inter KK^+$, $pi.alt = 1 - zeta_p$ and $nu = pi.alt overline(pi.alt)$. An _$cal(O)$-lattice of type $(z, s)$_ is a torsion-free $cal(O)$-module $Z$ of rank $z$ with a positive definite Hermitian form whose trace form $(x, y) = "Tr"_(KK slash QQ) H (x, y)$ is integral, with $Z^\# slash Z tilde.equiv FF_p^s$ and $zeta_p$ acting trivially on $Z^\# slash Z$; we use the rescaled form $tilde(H) = nu^((p-3) slash 2) H$ of Lemma 4.4, written $h$ in this appendix, so that $Z^* = pi.alt Z^\#$ and $Z^* subset Z subset pi.alt^(-1) Z^*$.

== D.1 The genus and its mass

#thm[
*Proposition D.1.* The $cal(O)$-lattices of type $(z, s)$ form at most one genus, and their trace forms are even. Here $0 <= s <= z$, since $Z^* = pi.alt Z^\# subset Z$ gives $[Z^\# : Z] <= [Z^\# : pi.alt Z^\#] = p^z$. The genus exists if and only if $z equiv s med (mod 2)$, and $s > 0$ or $(-1 slash p)^(z slash 2) = 1$. For $s = z$ the quadratic space $U slash pi.alt U$ with the form $h mod pi.alt$ has square discriminant, so for even $z$ it is hyperbolic if and only if $(-1 slash p)^(z slash 2) = 1$.
]
#proofbox[
_Evenness._ $KK slash KK^+$ has different $pi.alt cal(O)$, so $KK^+ slash QQ$ has different $pi.alt^(p - 3) cal(O) inter cal(O)^+$. For $x in Z$, $h(x, x) in KK^+ inter pi.alt^(-1) cal(O) = cal(O)^+$ (elements of $KK^+$ have even $pi.alt$-adic valuation), so $H (x, x) = nu^(-(p-3) slash 2) h(x, x)$ lies in the inverse different of $KK^+ slash QQ$, and $(x, x) = 2 "Tr"_(KK^+ slash QQ) H (x, x) in 2 ZZ$.

_Places other than $pi.alt$._ There $Z_v = Z_v^*$. If $v$ splits in $KK slash KK^+$, a unimodular Hermitian lattice has a basis in which $h$ is the identity. If $v$ is inert with residue extension $FF_(q^2) slash FF_q$, $h(x, x)$ is not identically 0 modulo $v$ (otherwise $h(x, y) + overline(h(x, y)) equiv 0$ for all $x, y$, and replacing $y$ by $c y$ gives $h equiv 0$), so some $h(x, x)$ is a unit, $chevron.l h(x, x) chevron.r$ splits off, and since the extension is unramified every unit of $cal(O)^+_v$ is a norm. Hence $Z_v tilde.equiv chevron.l 1, dots, 1 chevron.r$.

_The place $pi.alt$._ Let $Lambda = Z_pi.alt$. The form $a(x, y) = pi.alt h(x, y) mod pi.alt$ on $Lambda slash Lambda^*$ is non-degenerate and alternating (Lemma 4.4), so $z - s = 2 m$. Take $x, y$ with $h(x, y) = pi.alt^(-1)$; then $cal(O) x + cal(O) y$ has Gram matrix $pi.alt^(-1)$ times a unimodular matrix and is an orthogonal summand. Replacing $x$ by $x + r pi.alt y$ with $r in cal(O)^+_pi.alt$ chosen by Hensel's lemma (the derivative of the norm in $r$ is $-2$ modulo $pi.alt$), we get $h(x, x) = 0$, and similarly $h(y, y) = 0$. So $Lambda = R perp P^m$ with $P$ the plane with Gram matrix $mat(0, pi.alt^(-1); overline(pi.alt)^(-1), 0)$ and $R$ unimodular of rank $s$. The form $h mod pi.alt$ on $R slash pi.alt R$ is a non-degenerate symmetric form over $FF_p$; it can be diagonalized, $chevron.l u chevron.r$ depends only on the square class of $u mod pi.alt$ (units $equiv 1 mod pi.alt$ are norms), and $chevron.l u, v chevron.r tilde.equiv chevron.l 1, u v chevron.r$. So $R tilde.equiv chevron.l 1, dots, 1, d chevron.r$, and $Lambda$ is determined by $z$, $s$ and $chi_pi.alt (d)$, where $chi_pi.alt$ is the norm residue symbol of $KK_pi.alt slash KK^+_pi.alt$ on units, $chi_pi.alt (u) = (u mod pi.alt slash p)$.

_Product formula._ Let $delta = det h in KK^+$. At the real places $delta > 0$ is a norm, and at the finite places $v != pi.alt$, $chi_v (delta) = 1$ by the above. By the product formula $chi_pi.alt (delta) = 1$. Since $det P = -nu^(-1)$ and $nu$ is a norm, $chi_pi.alt (det P) = chi_pi.alt (-1) = (-1 slash p)$, hence $chi_pi.alt (d) = (-1 slash p)^m$. If $s > 0$ this determines $Lambda$; if $s = 0$ it requires $(-1 slash p)^(z slash 2) = 1$. So all completions are determined by $(z, s)$, and two such lattices lie in one genus. (Their spaces are isometric: they are locally isometric everywhere, so their determinants agree modulo local norms at every place, hence modulo global norms by the Hasse norm theorem for the cyclic extension $KK slash KK^+$, and [Shi08, Thm 2.2(i)] applies.) Conversely, if the conditions hold, the local norm classes of the determinants satisfy the product formula, so by global class field theory (the norm residue symbols of an element of $KK^+$ can be prescribed at all places, equal to 1 at almost all, subject only to the product formula) there is $delta in KK^+$, totally positive, with these local classes; [Shi08, Thm 2.2(ii)] gives a positive definite Hermitian space with determinant class $delta$, and a lattice with the prescribed completions exists in it. For $s = z$, $m = 0$ and $chi_pi.alt (d) = 1$, so the discriminant of $U slash pi.alt U$ is a square. #qed
]

The standard mass of [Jür15, §3.3.2] is
#tagged("D.1", $ "std"_z = 2^(1 - z r) product_(j = 1)^z |L(1 - j, chi^j)| dot cases(1 & (z "even"), 2^(-1) & (z "odd")), quad r = (p - 1) / 2, $)
where $chi$ is the quadratic character of $KK slash KK^+$, $chi^j$ is the trivial character for even $j$, and $L(1 - j, chi^j)$ is the $L$-function of $KK^+$. By [Jür15, Satz 3.3.15] it is the mass of the genus $G_z^KK$ of [Jür15, §3.3.1], whose local shape at the ramified place is hyperbolic for even $z$ and hyperbolic plus one rank-1 summand for odd $z$; in our notation this is the genus of type $(z, z mod 2)$, whenever that genus exists.

For other $s$ we count pairs. Let $U$ be of type $(z, z)$ and $Z$ of type $(z, s)$, $z - s = 2 m$. Then $U subset Z$ with $Z^* subset U$ if and only if $W = pi.alt Z slash pi.alt U$ is an $m$-dimensional totally isotropic subspace of the quadratic space $U slash pi.alt U$: indeed $h(pi.alt^(-1) u, pi.alt^(-1) u') = h(u, u') slash nu$ lies in $pi.alt^(-1) cal(O)$ if and only if $h(u, u') equiv 0 mod pi.alt$, and then $Z^* = {x in U : h(x, u) in pi.alt cal(O) "for all" u in W}$ has index $p^m$ in $U$, so $[Z : Z^*] = p^(2 m)$ and $Z$ has type $(z, s)$. Conversely every $Z$ contains $product_(i = 1)^m (p^i + 1)$ such $U$, one for each Lagrangian subspace of $Z slash Z^*$ (Lemma 4.4). By Proposition D.1 all these $Z$ lie in one genus, so with $I(z, m)$ the number of $m$-dimensional totally isotropic subspaces of $U slash pi.alt U$, which depends only on $z$, $m$ and the type of $U slash pi.alt U$ given by Proposition D.1,
#tagged("D.2", $ "mass"(z, s) dot product_(i = 1)^m (p^i + 1) = "mass"(z, z) dot I(z, m). $)
For a non-degenerate quadratic space of dimension $n$ over $FF_p$, $I(n, m) = product_(i = 0)^(m - 1) (p^(2(k - i)) - 1) slash (p^(i+1) - 1)$ if $n = 2 k + 1$, and $I(n, m) = product_(i = 0)^(m - 1) (p^(k - i) - epsilon)(p^(k - i - 1) + epsilon) slash (p^(i+1) - 1)$ if $n = 2 k$ and the space has type $epsilon = plus.minus 1$. Taking $s = z mod 2$ in (D.2) gives $"mass"(z, z)$ from (D.1) whenever the genus $(z, z mod 2)$ exists, that is, unless $z equiv 2 med (mod 4)$ and $p equiv 3 med (mod 4)$; and then (D.2) gives every $"mass"(z, s)$. For $p = 5, 17, 31$ this covers all genera we need. For $p = 7$ the only genus we need is $(4, 4)$ (Table D1), and (D.1)–(D.2) give $"mass"(z, s)$ for $z = 1, 3, 4$. But the triangular system (D.3) for $z = 3, 4$ needs $"mass"^emptyset (2)$, which (D.3) for $z = 2$ takes from $"mass"(2, 2)$, and (D.5) for $(4, 4)$ uses $"mass"(2, 2)$ for the complement of $A_6^2$; the genus $(2, 0)$ is empty. We determine the genus $(2, 2)$ directly. Its trace lattices are even, of rank 12 and determinant $7^2$. A rootless one would have minimum at least 4 and so centre density at least $7^(-1) approx 0.143$, more than the bound $0.06279$ of [CE03, Table 3] in dimension 12. So every $U$ of type $(2, 2)$ has a root, and by (4.2) $U = A_6^k perp U'$ with $k >= 1$ and $U'$ rootless of type $(2 - k, 2 - k)$. If $k = 1$, the trace lattice of $U'$ is even, of rank 6 and determinant 7, with minimum at least 4, so its centre density is at least $7^(-1 slash 2) approx 0.378$, more than the bound $0.08084$ of [CE03, Table 3] in dimension 6. Hence $k = 2$, the genus $(2, 2)$ consists of $A_6^2$ alone, and by (4.2) $"mass"(2, 2) = 1 slash ((2 dot 7)^2 dot 2!) = 1 slash 392$. The computation #calc("C2", "C9") takes $"mass"(z, z)$ for even $z$ from the closed form $2^(-z r) product_(j = 1)^z |L(1 - j, chi^j)| dot (p^(z slash 2) + (frac(-1, p))^(z slash 2))$; for $p >= 5$ it agrees with (D.1)–(D.2) in every case we use in which the genus $(z, 0)$ exists, and for $p = 7$, $z = 2$ it gives $1 slash 392$, the value just proved. For $p = 3$ and $z equiv 2 med (mod 4)$ we use instead Theorem D.5 of §D.4, which computes the masses of all genera over the Eisenstein integers from the masses of maximal lattices [GHY01]; for $z = 6$ and $z = 10$ it gives $"mass"(z, z) = 13 slash "235,146,240"$ and $"8,899" slash "31,345,665,638,400"$. As an independent check, these values agree with the mass formula for unimodular Hermitian lattices over $ZZ[zeta_3]$ of [HK89, Thm 5.6] (also quoted in [Jür15, §3.7.1]), and [Jür15, Table 4.8], reproducing Feit's classification, has exactly one rootless class for $z = 6$, the Coxeter–Todd lattice, and exactly two for $z = 10$, as (D.3) predicts (the misprints in the mass column of that table, at $N = 2$ and $N = 11$, do not concern $z = 6, 10$).

== D.2 Rootless unimodular lattices

By (4.2) every $U$ of type $(z, z)$ is uniquely $A_(p-1)^k perp U'$ with $U'$ rootless of type $(z - k, z - k)$, and conversely every such sum is of type $(z, z)$ and so lies in the genus. Since $|"U"(U)| = (2 p)^k k! |"U"(U')|$, summing over the genus gives
#tagged("D.3", $ "mass"(z, z) = sum_(k = 0)^z ("mass"^emptyset (z - k)) / ((2 p)^k k!), $)
where $"mass"^emptyset (n)$ is the sum of $1 slash |"U"(U')|$ over the rootless $U'$ of type $(n, n)$ ($"mass"^emptyset (0) = 1$). This triangular system gives $"mass"^emptyset (z)$ from the masses of §D.1 #calc("C2"). In particular $"mass"^emptyset (z) = 0$ for the ranks used in reason 3 of Proposition 3.3, and for $p = 31$, $z = 1$, $"mass"(1, 1) = 9 slash 62$ and $"mass"^emptyset (1) = 9 slash 62 - 1 slash 62 = 4 slash 31$.

== D.3 Rootless masses for $p = 5, 7, 17$

Let $Z$ be an $cal(O)$-lattice of type $(z, s)$ with $z (p - 1) <= 32$, and $p >= 5$. The roots of $Z$ (vectors of trace norm 2) span a root lattice on which $zeta_p$ acts without fixed vectors.

#thm[
*Lemma D.2.* The irreducible components of this root lattice are $A_(p-1)$ (of rank 1 over $cal(O)$) and, for $p = 5$ only, $E_8$ (of rank 2).
]
#proofbox[
If $zeta_p$ permuted $p$ components cyclically, the sum of the orbit of a root would be a non-zero fixed vector. So each component $C$ is preserved, and $zeta_p |_C$ is an automorphism of order $p$ without fixed vector; its characteristic polynomial is a power of $Phi_p$, so $p - 1$ divides the rank $n$ of $C$. For $C = A_n$, $"Aut" C = S_(n+1) times {plus.minus 1}$, an element of order $p$ is a permutation, and its fixed space in $A_n$ has dimension one less than its number of cycles, so it is a single $p$-cycle and $n = p - 1$. For $C = D_n$ an element of order $p$ is a signed permutation ($|"Aut" D_4| = 2^7 3^2$ has no element of order $p >= 5$) which, after conjugation, has all signs $+$ on its $p$-cycles, and it fixes the sum of the coordinates of a cycle. For $E_6$, $p - 1 | 6$ gives $p = 7$, but $7$ does not divide $|"Aut" E_6|$; $p - 1 | 7$ is impossible for $E_7$; for $E_8$, $p - 1 | 8$ gives $p = 5$, and such elements exist (characteristic polynomial $Phi_5^2$). #qed
]

On $A = A_(p-1)$ every root generates $A$ as an $cal(O)$-module, so its root sublattices are $0$ and $A$, and $c(A) = -1$ in the notation of Proposition C.2; $|"U"(A)| = 2 p$. In $E_8$ with an element of order 5 of characteristic polynomial $Phi_5^2$, the $cal(O)$-submodules generated by roots are $0$, twelve copies of $A_4$, six copies of $A_4^2$ (the orthogonal complement of $A_4$ in $E_8$ is again $A_4$), and $E_8$, so $c(E_8) = -(1 - 12 + 6) = 5$, and $|"U"(E_8)| = 600$. The proof of Proposition C.2 applies verbatim over $cal(O)$ and gives
#tagged("D.4", $ "mass"^emptyset (z, s) = sum_(a, b) (-1)^a 5^b (A(A^a E_8^b)) / ((2 p)^a a! dot 600^b b!), $)
where $A(T) = sum_Z r(T, Z) slash |"U"(Z)|$ over the genus $(z, s)$, $r(T, Z)$ is the number of $cal(O)$-linear isometric embeddings $T -> Z$, and $b = 0$ unless $p = 5$.

*Computing $A(T)$.* An embedding $iota : T -> Z$ extends to $T times.o QQ$, and $iota^(-1)(Z)$ is an overlattice $T_C$ of $T$ given by an isotropic subspace $C$ of $D(T) = FF_p^a$ [Nik79, Prop. 1.4.1] (all subspaces are $cal(O)$-stable, since $cal(O)$ acts on $D(T)$ through $cal(O) slash pi.alt = FF_p$). So $A(T) = sum_C A_"prim" (T_C)$ with $A_"prim"$ counting primitive embeddings. For a primitive $T' subset Z$ of rank $n$ and $D(T') = FF_p^(s')$, the complement $Y = T'^perp$ is an $cal(O)$-lattice of type $(z - n, s - s' + 2 k)$ (by the argument given in §D.4.4 for $p = 3$, which applies to every odd $p$) and $Z$ is the gluing of $T' perp Y$ along an anti-isometry $gamma$ from a $k$-dimensional subspace $cal(H) subset D(T')$ onto a subspace of $D(Y)$ [Nik79, §1.5, Prop. 1.5.1]; conversely every such gluing is of type $(z, s)$, hence in the genus (Proposition D.1). Counting the data $(cal(H), gamma, Y)$ with weights $1 slash |"U"(Y)|$ as in Proposition 4.1,
#tagged("D.5", $ A_"prim" (T') = sum_cal(H) ("emb"(cal(H), D(T')) dot "emb"(cal(H)^-, D(Y_cal(H)))) / (|O(cal(H))|) dot "mass"(z - n, s - s' + 2 dim cal(H)), $)
where $cal(H)$ runs over the isometry types of (possibly degenerate) quadratic spaces over $FF_p$, $"emb"(cal(H), D)$ is the number of isometric embeddings of $cal(H)$ into $D$, computed by Witt's theorem, $cal(H)^-$ is $cal(H)$ with the form negated, and $D(Y_cal(H))$ is the discriminant form of the genus of $Y$, determined by Proposition D.1. With the masses of §D.1 this gives every $"mass"^emptyset (z, s)$ in exact arithmetic #calc("C9").

#ptable(
  columns: 4, align: (left, right, right, left),
  caption: [Table D1. Rootless masses of the genera $(z, s)$ over $ZZ[zeta_p]$ used in §4.2 and Appendix B.4 #calc("C9", "C2").],
  table.header([$p$, $(z, s)$], [$"mass"^emptyset (z, s)$], [approx.], [used for]),
  [5, $(4, 4)$], [$1 slash "72,000"$], [$1.39 times 10^(-5)$], [$5 - (4, 16) - 4$],
  [5, $(5, 5)$], [$1 slash "150,000"$], [$6.67 times 10^(-6)$], [$5 - (5, 12) - 5$],
  [5, $(6, 6)$], [$191 slash "7,200,000"$], [$2.65 times 10^(-5)$], [$5 - (6, 8) - 6$],
  [5, $(6, 4)$], [$"4,567" slash "1,440,000"$], [$3.17 times 10^(-3)$], [$5 - (6, 8) - 4$],
  [5, $(7, 7)$], [$"10,358,729" slash "11,340,000,000"$], [$9.13 times 10^(-4)$], [$5 - (7, 4) - 3$ (first step)],
  [5, $(7, 5)$], [$"16,927" slash "28,800"$], [$0.588$], [$5 - (7, 4) - 3$ (second step)],
  [5, $(7, 3)$], [$"799,867" slash "1,440,000"$], [$0.5555$], [$5 - (7, 4) - 3$],
  [5, $(8, 0)$], [$"168,947" slash "138,240,000"$], [$1.22 times 10^(-3)$], [$5 - (8, 0) - 0$],
  [7, $(4, 4)$], [$439 slash "395,136"$], [$1.11 times 10^(-3)$], [$7 - (4, 8) - 4$],
  [17, $(2, 0)$], [$1 slash 12$], [$0.0833$], [$17 - (2, 0) - 0$],
  [31, $(1, 1)$], [$4 slash 31$], [$0.129$], [$31 - (1, 2) - 1$],
)

*Checks.* For $p = 5$, $(2, 0)$ and $(4, 0)$, (D.4) gives that $E_8$ and $E_8^2$ are the only lattices, with masses $1 slash 600$ and $1 slash "720,000" = 1 slash (600^2 dot 2)$. For $s = z$, (D.4) agrees with (D.3). Three values are compared with the Leech lattice $Lambda$: an element $g in "Co"_0$ with $g^p = 1$ and no fixed vector makes $Lambda$ a rootless lattice of type $(24 slash (p - 1), 0)$, and the genus contains no other rootless lattice since $Lambda$ is the only rootless even unimodular lattice of rank 24; so the rootless mass is $\#{"such" g "up to conjugacy, weighted"} = sum 1 slash |C_("Co"_0)(g)|$, read from the character table of $2."Co"_1$ in GAP. The values $1 slash "6,048,000"$ ($p = 5$, $(6, 0)$), $1 slash "35,280"$ ($p = 7$, $(4, 0)$) and $1 slash 312$ ($p = 13$, $(2, 0)$) agree with (D.4) #calc("C9"). Finally, for $17 - (2, 0) - 0$ the value $1 slash 12$ is the sum of $1 slash |"U"(Z)|$ over the six Hermitian classes whose trace lattices are the four lattices of [Jür15, Satz 3.6.9] (§4.2).

== D.4 The case $p = 3$ and the genera of $F$

#policy-box[
*What this section does.* For $p = 3$ the root lattices over $cal(O) = ZZ[omega]$ include $D_4$, whose discriminant group has a 2-part, so the orthogonal complements of root sublattices of $Z$ leave the genera of §D.1 and acquire a 2-adic level. §D.4.1 describes these genera $(N, s, b)$ and §D.4.2 proves a closed formula for their masses (Theorem D.5), anchored at the maximal lattices of [GHY01]; in particular it gives the masses of the genera $(z, z)$ for $z = 6, 10$, where the genus $(z, 0)$ is empty. §D.4.3 lists the root lattices over $ZZ[omega]$ and §D.4.4 computes the rootless masses as in §D.3; §D.4.5 gives the value for $(13, 5)$ and the checks. §D.4.6 treats the $ZZ$-lattices $F$ of the genera $(20, 6)$ for $p = 3$ and $(16, 4)$ for $p = 5$ by the same inclusion–exclusion over root lattices, with the masses of the complements given by Siegel's mass formula; the new point there is a local formula (Theorem D.10) that expresses the local factor of a complement at a prime where $F$ is unimodular by counting on $F slash p F$.
]

Throughout §D.4.1–§D.4.5, $p = 3$, $cal(O) = ZZ[omega]$ with $omega = zeta_3$, $pi.alt = 1 - omega$ and $chi = chi_(-3)$; for $p = 3$ the form $tilde(H)$ of Lemma 4.4 is $H$, written $h$ below, and $Z^\# = {x : h(x, Z) subset pi.alt^(-1) cal(O)}$. The prime 2 is inert in $QQ(omega)$ and $cal(O) slash 2 = FF_4$. A _root_ is a vector of trace norm 2, that is, with $h(x, x) = 1$.

=== D.4.1 Genera with a 2-adic level

An _$cal(O)$-lattice of type $(N, s, b)$_ is a free $cal(O)$-module $Lambda$ of rank $N$ with a positive definite Hermitian form $h$ such that $h(Lambda, Lambda) subset pi.alt^(-1) cal(O)$, and, with $Lambda^\# = {y : h(y, Lambda) subset pi.alt^(-1) cal(O)}$:
- at 3: $pi.alt Lambda_3^\# subset Lambda_3$ and $Lambda_3^\# slash Lambda_3 tilde.equiv FF_3^s$;
- at 2: $2 Lambda_2^\# subset Lambda_2$ and $Lambda_2^\# slash Lambda_2 tilde.equiv FF_4^b$;
- $Lambda_v^\# = Lambda_v$ at all other finite places.
Type $(N, s, 0)$ is type $(N, s)$ of §D.1. As there, the trace form is even. On $D_2 (Lambda) = Lambda_2^\# slash Lambda_2$ the form $beta_2 (x, y) = 2 h(x, y) mod 2 cal(O)_2$ is a non-degenerate Hermitian form over $FF_4$ ($beta_2 (x, Lambda_2^\#) = 0$ means $h(x, Lambda_2^\#) subset cal(O)_2$, that is $x in Lambda_2$); on $D_3 (Lambda) = Lambda_3^\# slash Lambda_3$ the form $beta_3 (x, y) = 3 h(x, y) mod pi.alt$ is a non-degenerate symmetric form over $FF_3$. Since the trace form is $(x, y) = 2 "Re" h(x, y)$, the discriminant form of the trace lattice is $q(x) = 2 h(x, x) mod 2 ZZ$, and $beta_3 (x, x) = (3 slash 2) q(x) mod 3$ is the value $Q(x)$ of the convention of §B.1.1 (Table 3).

#thm[
*Lemma D.3.* Let $Lambda$ be a lattice over $cal(O)_2 = ZZ_2 [omega]$ with a Hermitian form, $2 Lambda^\# subset Lambda subset Lambda^\#$ and $Lambda^\# slash Lambda tilde.equiv FF_4^b$. Then $Lambda tilde.equiv chevron.l 1 chevron.r^(N - b) perp chevron.l 2 chevron.r^b$.
]
#proofbox[
On $W_2 = Lambda slash 2 Lambda^\#$ the form $h mod 2$ is a non-degenerate Hermitian form over $FF_4$. If $W_2 != 0$, some $h(x, x)$ is a unit: otherwise $h(x, y) + h(y, x) equiv 0$ for all $x, y$, and replacing $y$ by $omega y$ gives $h equiv 0$. Since $QQ_2 (omega) slash QQ_2$ is unramified, every unit of $ZZ_2$ is a norm of a unit, so $chevron.l 1 chevron.r$ splits off. When $W_2 = 0$, $Lambda = 2 Lambda^\#$ and $h slash 2$ is unimodular on $Lambda$; the same argument gives $chevron.l 2 chevron.r^b$. #qed
]

#thm[
*Proposition D.4.* The $cal(O)$-lattices of type $(N, s, b)$ form at most one genus. The genus exists if and only if $0 <= s <= N$, $N equiv s med (mod 2)$, $0 <= b <= N$, and $s > 0$ or $(-1)^(N slash 2 + b) = 1$. Write $m = (N - s) slash 2$. The form $beta_3$ on $D_3 (Lambda)$ has determinant $(-1)^(m + b)$ modulo squares; for even $s$ its type (hyperbolic $+1$, anisotropic part of dimension 2 $-1$) is $tau = (-1)^(N slash 2 + b)$.
]
#proofbox[
At 3, the proof of Proposition D.1 gives $Lambda_3 tilde.equiv R perp P^m$ with $R tilde.equiv chevron.l 1, dots, 1, d chevron.r$ unimodular of rank $s$, determined by $N$, $s$ and $chi_pi.alt (d)$; the form $h mod pi.alt$ on $R slash pi.alt R$ is isometric to $beta_3$ (send $x in R$ to $x slash pi.alt in R^\#$; then $3 h(x slash pi.alt, y slash pi.alt) = h(x, y)$). At 2, Lemma D.3 determines $Lambda_2$ by $b$. At the other places $Lambda_v$ is unimodular, hence unique (proof of Proposition D.1). For the product formula, let $delta = det h$. At 2, $delta in 2^b ZZ_2^times$, and since $QQ_2 (omega)$ is unramified the norms are the elements of even valuation, so $chi_2 (delta) = (-1)^b$. The other places $v != 3$ contribute 1, and at 3, $chi_pi.alt (delta) = chi_pi.alt (d) dot (-1 slash 3)^m$ as in the proof of Proposition D.1. Hence $chi_pi.alt (d) = (-1)^(m + b)$, which is the determinant of $beta_3$ modulo squares; for $s = 2 l$, $tau = (-1)^l det beta_3 = (-1)^(l + m + b) = (-1)^(N slash 2 + b)$. For $s = 0$ this is the condition $(-1)^(N slash 2 + b) = 1$. Uniqueness of the genus and existence follow as in Proposition D.1, with [Shi08, Thm 2.2]. #qed
]

=== D.4.2 The masses of the genera $(N, s, b)$

For $m >= 0$ let $u_m = product_(i = 1)^m (2^i - (-1)^i)$ and $g_m = product_(i = 1)^m (9^i - 1)$. Put
$ lambda_2 (N, b) = u_N / (u_(N - b) u_b), quad lambda_3 (N, s) = (2 g_(floor(N slash 2))) / (g_m o_s), quad o_0 = 1, quad o_(2 l + 1) = 2 g_l, quad o_(2 l) = 2 (3^l - tau) g_(l - 1) quad (l >= 1), $
with $m = (N - s) slash 2$ and $tau$ as in Proposition D.4. Here $|"U"_b (FF_4)| = 2^(b (b - 1) slash 2) u_b$, and $o_s$ is the order of the orthogonal group of $(D_3, beta_3)$ divided by its power of 3.

#thm[
*Theorem D.5.* Let $N >= 1$. If the genus of type $(N, s, b)$ exists, its mass is
#tagged("D.6", $ "mass"(N, s, b) = 2^(-N) product_(j = 1)^N |L(1 - j, chi^j)| dot lambda_3 (N, s) dot lambda_2 (N, b). $)
]

At rank zero we instead set $"mass"(0,0,0)=1$, the mass of the zero lattice, and give non-existing types mass zero. The empty-product expression (D.6) is not applied at rank zero. These conventions also apply to the zero-rank complements in (D.9).

For $b = 0$ and $s = N mod 2$ this is (D.1): $lambda_3 = 2$ for $s = 0$ and $lambda_3 = 1$ for $s = 1$. Theorem D.5 covers every existing type $(N, s, b)$ of positive rank, including the genera $(z, z)$ with $z equiv 2 med (mod 4)$, where it is anchored at the maximal lattices with $s = 2$, $tau = -1$; for instance $"mass"(6, 6, 0) = 13 slash "235,146,240"$ and $"mass"(10, 10, 0) = "8,899" slash "31,345,665,638,400"$. Being the true masses, the values (D.6) with $b = 0$ satisfy (D.2).

The proof has two steps. The genera with $b <= 1$ and $s <= 2$ that contain maximal lattices are given by [GHY01]; the others are reached from them by counting pairs of lattices of index 3 (Lemma D.6) and of index 4 (Lemma D.7).

*The maximal lattices.* [GHY01] computes the mass of a genus of lattices that are maximal among lattices with $h(x, x) in ZZ$ (for $k = QQ$, $d = 1$, absolute rank $N$):
#tagged("D.7", $ "Mass"(Lambda) = 2^(-N) |L(M)| dot tau(G) dot product_(v in S) lambda_v, quad L(M) = product_(j = 1)^N L(1 - j, chi^j), quad tau(G) = 2, $)
[GHY01, (2.12), Prop. 2.13, (4.2), (4.3), Prop. 4.4, 4.5], with $lambda_v = q^(-N(overline(H)^0)) |overline(H)^0 (FF_q)| slash q^(-N(overline(G))) |overline(G)(FF_q)|$ and $S$ the set of places where the reduction $overline(G)$ of the stabilizer differs from that of the quasi-split model. The mass is positive while $L(M) < 0$ for $N equiv 2 med (mod 4)$, so the product of $L$-values is read in absolute value. The local factors needed are [GHY01, Tables 1, 2]:

#ptable(
  columns: 4, align: (left, left, left, left),
  caption: [Table D2. The maximal lattices and their local factors in [GHY01]. The last column is the factor $lambda_3 (N, s)$, $lambda_2 (N, b)$ of (D.6).],
  table.header([local type], [$overline(G)$], [$lambda_v$ of [GHY01]], [factor in (D.6)]),
  [$p = 2$, $b = 0$], [$"U"_N$], [$1$ ($2 in.not S$)], [$lambda_2 (N, 0) = 1$],
  [$p = 2$, $b = 1$, $N$ odd], [$"U"_(N - 1) times "U"_1$], [$(2^N + 1) slash 3$], [$lambda_2 (N, 1)$],
  [$p = 2$, $b = 1$, $N$ even], [$"U"_(N - 1) times "U"_1$], [$(2^N - 1) slash 3$], [$lambda_2 (N, 1)$],
  [$p = 3$, $N$ odd, $s = 1$], [$ZZ slash 2 times "Sp"_(N - 1)$], [$1 slash 2$], [$lambda_3 = 2 dot 1 slash 2 = 1$],
  [$p = 3$, $N$ even, $s = 0$, $tau = +1$], [$"Sp"_N$], [$1$ ($3 in.not S$)], [$lambda_3 = 2$],
  [$p = 3$, $N$ even, $s = 2$, $tau = -1$], [$"Sp"_(N - 2) times O_2^-$], [$(3^N - 1) slash 8$], [$lambda_3 = (3^N - 1) slash 4$],
)

The factor 2 in $lambda_3$ is $tau(G)$. To identify these lattices with our types it suffices, by Lemma D.3 and Proposition D.4, to compare the numerical invariants. The lattice $Delta = A_E e plus.o frak(D)^(-1) f$ of [GHY01, (3.2)] is the plane $P$ at 3 (after replacing $f$ by $-f$) and is unimodular at 2. The lattice $A_D$ at 3 with $alpha = -1$ is $cal(O) plus.o cal(O) z$ with the form $chevron.l 1, 1 chevron.r$: the unit parts of $N(x_1)$ and $N(x_2)$ are $1 mod 3$, so their sum does not cancel the denominator; $R slash pi.alt R$ is anisotropic, so $s = 2$, $tau = -1$. At 2 with $alpha = 2$ it is $cal(O)_2 plus.o cal(O)_2 z$ with the form $chevron.l 1, -2 chevron.r tilde.equiv chevron.l 1, 2 chevron.r$, so $b = 1$. Conversely the lattices with $b >= 2$ ($chevron.l 2, 2 chevron.r$ lies in a unimodular lattice), $s >= 3$, or $s = 2$, $tau = +1$ (an isotropic vector of $D_3$ gives an overlattice) are not maximal; [GHY01] treats only these extremal vertices (p. 102). In each row of Table D2 the last column equals the corresponding factor of (D.6), so (D.6) holds for the genera $(N, s_0, b_0)$ with $b_0 in {0, 1}$ and $s_0 = 1$ ($N$ odd), $s_0 = 0$ ($N$ even, $tau = +1$), $s_0 = 2$ ($N$ even, $tau = -1$).

For a non-degenerate quadratic space over $FF_3$ of dimension $n$ (and type $tau$ if $n = 2 l$) let $I_O (n)$ be its number of isotropic lines, $(3^(n - 1) - 1) slash 2$ for odd $n$ and $(3^l - tau)(3^(l - 1) + tau) slash 2$ for $n = 2 l$; for a non-degenerate Hermitian space over $FF_4$ of dimension $n$ let $I_"U" (n) = (2^n - (-1)^n)(2^(n - 1) - (-1)^(n - 1)) slash 3$ be its number of isotropic lines.

#thm[
*Lemma D.6.* If $N - s >= 2$, then $"mass"(N, s + 2, b) dot I_O (s + 2) = "mass"(N, s, b) dot (3^(N - s) - 1) slash 2$, where $I_O (s + 2)$ is taken for the form of type $(N, s + 2, b)$.
]
#proofbox[
Count the pairs $Lambda' subset Lambda$ of $cal(O)$-lattices of index 3 with $Lambda$ of type $(N, s, b)$ and $Lambda'$ of type $(N, s + 2, b)$, weighted by $1 slash |"U"(Lambda)|$ and $1 slash |"U"(Lambda')|$; by the orbit–stabilizer relation both weighted counts are equal. Only the completions at 3 change.

_Sublattices._ A sublattice of index 3 is the kernel of a non-zero linear form $f$ on $Lambda slash pi.alt Lambda$. If $Lambda'$ satisfies $pi.alt Lambda'^\# subset Lambda'$, then $pi.alt Lambda^\# subset pi.alt Lambda'^\# subset Lambda'$, so $f$ vanishes on $pi.alt Lambda^\# slash pi.alt Lambda$ and is a form on $W = Lambda slash pi.alt Lambda^\#$. Conversely, for such $f$, non-degeneracy of $a$ on $W$ (Lemma 4.4) gives $f = a(dot, v)$ with $v in Lambda without pi.alt Lambda^\#$, and $Lambda' = {x in Lambda : h(x, v) in cal(O)}$. Then $v slash pi.alt in Lambda'^\# without Lambda^\#$, comparing indices $Lambda'^\# = Lambda^\# + ZZ_3 v slash pi.alt$, and $pi.alt Lambda'^\# = pi.alt Lambda^\# + ZZ_3 v subset Lambda'$ because $h(v, v) in ZZ_3$; so $D_3 (Lambda') tilde.equiv FF_3^(s + 2)$. The forms $f$ correspond to the hyperplanes of $W$, of which there are $(3^(N - s) - 1) slash 2$.

_Overlattices._ If $Lambda' subset Lambda$, then $Lambda subset Lambda^\# subset Lambda'^\#$ and $I = Lambda slash Lambda'$ is a line in $D' = D_3 (Lambda')$. For $x in Lambda'^\#$, $Lambda' + ZZ_3 x$ is an $cal(O)$-module since $(omega - 1) x in pi.alt Lambda'^\# subset Lambda'$, and it is integral if and only if $h(x, x) in ZZ_3$, that is $beta_3 (x, x) = 0$. Then $pi.alt Lambda^\# subset Lambda$ holds automatically, and $D_3 (Lambda) = I^perp slash I$ has dimension $s$. So the overlattices correspond to the isotropic lines of $D'$.

_Types._ $D' tilde.equiv D_3 (Lambda) perp$ (hyperbolic plane), so for even $s$ the type $tau$ is unchanged; $b$ is unchanged. By Proposition D.4 the lattices obtained lie in the genera $(N, s, b)$ and $(N, s + 2, b)$. #qed
]

#thm[
*Lemma D.7.* If $N - b >= 2$, then $"mass"(N, s, b + 2) dot I_"U" (b + 2) = "mass"(N, s, b) dot I_"U" (N - b)$.
]
#proofbox[
Count pairs of index 4 at 2 as in Lemma D.6. A sublattice $Lambda'$ of index 4 is the kernel of a form $f = h(dot, v) mod 2$ on $W_2 = Lambda slash 2 Lambda^\#$ ($v in Lambda without 2 Lambda^\#$); the condition $2 Lambda'^\# = 2 Lambda^\# + cal(O)_2 v subset Lambda'$ holds if and only if $h(v, v) in 2 ZZ_2$, that is the line of $v$ is isotropic in $W_2$, which has dimension $N - b$: there are $I_"U" (N - b)$ of them. The overlattices of $Lambda'$ of index 4 with $2 Lambda^\# subset Lambda$ correspond to the isotropic lines of $D_2 (Lambda')$, of dimension $b + 2$. By Lemma D.3 and Proposition D.4 the lattices lie in the stated genera; $s$ and $tau$ are unchanged since the Hermitian space is the same. #qed
]

#proofbox[
(of Theorem D.5) The right side of (D.6) satisfies the relations of Lemmas D.6 and D.7: $lambda_3 (N, s + 2) slash lambda_3 (N, s) = (g_m o_s) slash (g_(m - 1) o_(s + 2))$, where $g_m slash g_(m - 1) = 9^m - 1 = 3^(N - s) - 1$ and $o_(s + 2) slash o_s$ equals $9^(l + 1) - 1$ for $s = 2 l + 1$, $(3^(l + 1) - tau)(3^l + tau)$ for $s = 2 l >= 2$ and $4$ for $s = 0$ ($tau = +1$), in each case $2 I_O (s + 2)$; and $lambda_2 (N, b + 2) slash lambda_2 (N, b) = u_(N - b) u_b slash (u_(N - b - 2) u_(b + 2)) = I_"U" (N - b) slash I_"U" (b + 2)$. It agrees with the maximal lattices (Table D2). Every existing genus is reached from one of these: from $(N, s_0, b_0)$ with $b_0 = b mod 2$, Lemma D.6 reaches every $s <= N$ with $N - s$ even (each step needs $N - s >= 2$, and the case $s = 2$, $tau = +1$ is one step from $s = 0$ with $I_O (2) = 2$; the case $s = 0$, $tau = -1$ is empty, consistent with $I_O (2) = 0$ for $tau = -1$), and Lemma D.7 then reaches every $b <= N$ with $b equiv b_0$ (each step needs $N - b >= 2$). In every step the genus on the far side is non-empty, since $I_O (n) > 0$ for $n >= 3$ and $I_"U" (n) > 0$ for $n >= 2$. The genus $(N, s_0, b_0)$ itself is non-empty: take maximal lattices at every place in $V = chevron.l 1 chevron.r^(N - 1) perp chevron.l c chevron.r$ with $c = 1$ for $b_0 = 0$ and $c = 2$ for $b_0 = 1$. #qed
]

=== D.4.3 Root lattices over $ZZ[omega]$

#thm[
*Lemma D.8.* The irreducible components of the root lattice of an $cal(O)$-lattice are $A_2$, $D_4$, $E_6$, $E_8$, each with a unique $cal(O)$-structure, with the data of Table D3.
]
#proofbox[
For a root $r$, $(r, omega r) = "Tr" overline(omega) = -1$, so $r$ and $omega r$ lie in the same component and $omega$ preserves each component $Y$, acting without fixed vectors ($1 + omega + omega^2 = 0$). For $A_n$, $n >= 2$, $"Aut" Y = W times {plus.minus 1}$ and an element of order 3 is a permutation in $S_(n + 1)$ whose fixed space has dimension one less than the number of cycles; so it is one 3-cycle and $n = 2$. $A_1$, $E_7$ and $D_n$ with $n$ odd have odd rank. For $D_n$, $n >= 5$, $"Aut" Y = W(B_n)$, and an element of order 3 is a signed permutation whose cycles have sign product $+1$, so it fixes a vector on each cycle. The four remaining lattices carry $cal(O)$-structures (explicitly inside $E_8 = {x in cal(O)^4 : x mod pi.alt in "tetracode"}$ with $h(x, y) = 1 / 3 sum x_i overline(y_i)$). Each such structure is an $cal(O)$-lattice of the type in Table D3, and by Theorem D.5 the masses of these genera are $1 slash 6$, $1 slash 72$, $1 slash "1,296"$, $1 slash "155,520"$, equal to $1 slash |"U"(Y)|$; so each genus consists of one class and the structure is unique. #qed
]

#ptable(
  columns: 6, align: (left, right, right, left, right, left),
  caption: [Table D3. Irreducible root lattices over $ZZ[omega]$ and their root sublattices (as $cal(O)$-modules), counted from explicit roots #calc("C16").],
  table.header([$Y$], [rank], [$|"U"(Y)|$], [type $(N, s, b)$], [$c(Y)$], [root sublattices]),
  [$A_2$], [1], [6], [$(1, 1, 0)$], [$-1$], [$emptyset$: 1, $A_2$: 1],
  [$D_4$], [2], [72], [$(2, 0, 1)$], [$3$], [$emptyset$: 1, $A_2$: 4, $D_4$: 1],
  [$E_6$], [3], [1,296], [$(3, 1, 0)$], [$-24$], [$emptyset$: 1, $A_2$: 12, $A_2^2$: 12, $D_4$: 9, $A_2^3$: 4, $E_6$: 1],
  [$E_8$], [4], [$"155,520"$], [$(4, 0, 0)$], [$324$], [$emptyset$: 1, $A_2$: 40, $A_2^2$: 240, $D_4$: 90, $A_2^3$: 160, $D_4 A_2$: 360, $E_6$: 40, $A_2^4$: 40, $E_6 A_2$: 40, $D_4^2$: 45, $E_8$: 1],
)

The discriminant modules are $D_3 (A_2) tilde.equiv FF_3$ with $beta_3 = chevron.l 1 chevron.r$, $D_3 (E_6) tilde.equiv FF_3$ with $beta_3 = chevron.l -1 chevron.r$, and $D_2 (D_4) tilde.equiv FF_4$; so for $T = A_2^a D_4^b E_6^c E_8^d$, $D_3 (T)$ has dimension $a + c$ and determinant $(-1)^c$, and $D_2 (T) tilde.equiv FF_4^b$. The proof of Proposition C.2 applies verbatim over $cal(O)$ and gives
#tagged("D.8", $ "mass"^emptyset (z, s) = sum_(a, b, c, d) ((-1)^a 3^b (-24)^c 324^d) / (6^a a! dot 72^b b! dot "1,296"^c c! dot "155,520"^d d!) A(A_2^a D_4^b E_6^c E_8^d), $)
the sum over $a + 2 b + 3 c + 4 d <= z$ (194 terms for $z = 13$), with $A(T) = sum_Z r(T, Z) slash |"U"(Z)|$ over the genus $(z, s)$.

=== D.4.4 The terms $A(T)$

As in §D.3, $A(T) = sum_C A_"prim" (T_C)$, now with $C = C_2 plus.o C_3$, $C_3$ a totally singular subspace of $(D_3 (T), beta_3)$ and $C_2$ a totally isotropic subspace of $(D_2 (T), beta_2)$; both are automatically $cal(O)$-stable, and every such $C$ gives an integral overlattice. For a primitive $T' subset Z$ of rank $n$ with $D_3 (T') tilde.equiv FF_3^(s')$ and $D_2 (T') tilde.equiv FF_4^(b')$ let $Y = T'^perp$. Then:
- $D_3 (Y)$ is killed by $pi.alt$: for $y in Y^\#$ choose $w = y + v in Z^\#$ with $v in T'^\#$; then $pi.alt w in Z$ and $pi.alt v in T' subset Z$, so $pi.alt y in Z inter QQ Y = Y$;
- at 2, $Z_2$ is unimodular, so $Z$ is glued along an anti-isometry of all of $D_2 (T')$ onto all of $D_2 (Y)$, and $Y$ has type $(z - n, s - s' + 2 dim cal(H), b')$ for the subspace $cal(H) subset D_3 (T')$ glued at 3;
- conversely every such gluing is an $cal(O)$-lattice of type $(z, s, 0)$ (subspaces of $FF_3$- and $FF_4$-spaces are $cal(O)$-stable since $omega equiv 1 mod pi.alt$ and $cal(O) slash 2 = FF_4$), hence lies in the genus.
The anti-isometries $D_2 (T') -> D_2 (Y)$ form a torsor under $"U"(D_2 (T')) tilde.equiv "U"_(b') (FF_4)$. Counting as in (D.5),
#tagged("D.9", $ A_"prim" (T') = |"U"_(b') (FF_4)| sum_cal(H) ("emb"(cal(H), D_3 (T')) dot "emb"(cal(H)^-, D_3 (Y_cal(H)))) / (|O(cal(H))|) dot "mass"(z - n, s - s' + 2 dim cal(H), b'), $)
with the masses of Theorem D.5 and $D_3 (Y_cal(H))$ the form of that genus (Proposition D.4). The numbers of totally singular subspaces of $D_3 (T)$ and of totally isotropic subspaces of $D_2 (T)$ and the embedding numbers in (D.9) are given by Witt's theorem and the orders of the finite orthogonal and unitary groups; all terms are exact rationals #calc("C16").

=== D.4.5 The genus $(13, 5)$ and checks

#thm[
*Theorem D.9.* $"mass"^emptyset (13, 5) = "27,547,211" slash "6,635,520" approx 4.1515$.
]

The mass of the genus $(13, 5)$ is $"mass"(13, 5, 0) approx 6.9221$, so the rootless part is $0.59974$ of it. Of the 194 root types, 51 have non-zero mass (the masses of the classes with a given root system are obtained from the $A(T)$ by the triangular inversion of (D.8) over all root systems), and none is negative #calc("C16"). The 105 classes of type $(13, 5)$ listed in Appendix B.4 are pairwise non-isometric and the sum of $1 slash |"U"(Z)|$ over them equals Theorem D.9 #calc("C16")\; as for Proposition 4.6, this shows that the list is complete. For type $3-(13,6)-5$, this component mass also determines the gluing-pair mass: $F=sqrt(3)E_6^\#$ and $"Aut" F$ maps onto $O(q_F)$.

#ptable(
  columns: 2, align: (left, left),
  caption: [Table D4. Checks of (D.6)–(D.9), in exact arithmetic #calc("C16"). They are consistency checks; the values rest on the proofs above.],
  table.header([check], [result]),
  [(D.6) with $b = 0$ against (D.1)–(D.2), $N <= 20$], [equal],
  [$"mass"^emptyset (z, z)$ by (D.8) against (D.3), $z <= 16$], [equal ($1 slash "39,191,040"$ for $z = 6$, $109 slash "451,480,780,800"$ for $z = 10$)],
  [$A(T)$ against direct counts inside $A_2$, $A_2^2$, $A_2^3$, $E_6$, $E_6 A_2$, $E_8$], [44 values equal],
  [$(9, 5)$, $(9, 7)$, $(10, 6)$, $(10, 8)$, $(11, 9)$, $(12, 8)$, $(12, 10)$, $(13, 13)$: sum over the complete lists of Appendix B.4], [equal (e.g. $"65,919,403,483" slash "2,612,138,803,200"$ for $(12, 8)$, 51 classes)],
  [$(12, 0)$: the 5 root systems with non-zero mass], [$emptyset$: $1 slash |6."Suz"|$; $E_8^3$, $E_6^4$, $D_4^6$, $A_2^(12)$: $1 slash |"U"|$ of the lattices $E_8^3$, $E_6^4$, $D_4^6$, $A_2^(12)$],
  [$(16, 0)$: rootless mass $"373,348,469" slash "360,972,288,000"$], [equal to the sum over the eleven lattices of [Hen13, Table 4]],
  [masses of the classes with a given root system, all $(z, s)$ with $z <= 16$], [none negative],
)

For $(12, 0)$ the only rootless lattice is the Leech lattice with an automorphism $g$ of order 3 without fixed vectors; in $2."Co"_1$ these form the single class $3 a$ (character value $-12$ on the 24-dimensional representation), with centralizer of order $"2,690,072,985,600" = 6 |"Suz"|$ (character table library of GAP), and (D.8) gives exactly $1 slash "2,690,072,985,600"$. A random walk on 2-neighbours of the genus $(13, 5)$, which does not use (D.8), gives the rootless fraction $0.5975 plus.minus 0.0047$ against $0.59974$; this is a statistical check only.

=== D.4.6 The $p$-elementary $ZZ$-lattices $F$

Let $p_0 in {3, 5}$, and let $F$ be an even $ZZ$-lattice of rank $d$ with $F^\# slash F tilde.equiv FF_(p_0)^s$; these lattices form one genus [Nik79, Cor. 1.9.4]. We need the rootless mass $mu_F^emptyset (d, s) = sum_(min F >= 4) 1 slash |"Aut" F|$ for $(p_0; d, s) = (3; 20, 6)$ and $(5; 16, 4)$.

*Inclusion–exclusion.* With the coefficients $c$ of King's inclusion–exclusion over the root lattices of type $A D E$ [Kin03, Prop. 1] (Proposition C.2 over $ZZ$),
#tagged("D.10", $ mu_F^emptyset (d, s) = sum_T c(T) (A(T)) / (|"Aut" T|), quad A(T) = sum_F (r(T, F)) / (|"Aut" F|), $)
the sum over the $"10,632"$ root lattices $T$ of rank at most 20 (for $d = 20$). As before $A(T) = sum_C A_"prim" (T_C)$ over the isotropic subgroups $C$ of the discriminant form of $T$, and for a primitive $T' subset F$ with complement $Y$ of rank $N = d - "rank" T'$, $F$ is glued from $T' perp Y$ [Nik79, §1.5, Prop. 1.5.1]. At a prime $p != p_0$, $F_p$ is unimodular, so the gluing at $p$ is an anti-isometry of all of $D_p (T')$ onto all of $D_p (Y)$, and $Y_p$ is determined by $T'_p$ [Nik79, Cor. 1.9.3]; at $p_0$ the gluing is along a subgroup $cal(H) subset D_(p_0) (T')$ as in (D.5). By Siegel's mass formula in the form [CS88, (2), (3)], the mass of the genus of $Y$ is $c_infinity (N) product_p 2 m_p (Y)$ with $c_infinity (N) = 2 pi^(-N (N + 1) slash 4) product_(j = 1)^N Gamma(j slash 2)$ and $m_p$ depending only on $Y_p$. Since both the gluing data and the local factors are local, the sum over $C$ and over the genera of $Y$ factors over the primes; each choice of local data comes from one global genus, since the rational complement exists globally by the Hasse–Minkowski theorem and Witt cancellation:
#tagged("D.11", $ A(T) = c_infinity (N) product_p Psi_p (T), quad Psi_p (T) = sum_(C_p) sum_((cal(H), gamma) "at" p) 2 m_p (Y_(C_p, cal(H))), $)
where for $p in.not S = {2, p_0} union {p | det T}$ the only term is $C_p = 0$ with $Y_p$ unimodular, and the product of these factors with $c_infinity (N)$ is the standard mass of [CS88, (6)–(8)] (written there with $zeta_D$; we use the primitive $L(s, chi_(D_0))$ and remove the factors at $S$, which gives the same value). For $N = 1$ the leading 2 of [CS88, (2)] is replaced by 1 [CS88, §6], and $N = 0$ is counted directly by overlattices.

At $p = p_0$, $Psi_p$ is computed from [CS88, (3), (5), Table 2], since $Y_p$ is a Jordan sum of components of scales $1$, $p$ and (when $T$ has a component $A_8$ or $A_(17)$) $p^2$. At $p != p_0$ we avoid the Jordan decomposition of $Y_p$, which at $p = 2$ involves the type I and bound components of [CS88, Table 1], by the following local formula.

#thm[
*Theorem D.10.* Let $L_p$ be a unimodular $ZZ_p$-lattice of rank $d$ (even if $p = 2$), $T' subset L_p$ a primitive sublattice of rank $n$ with orthogonal complement $Y$ of rank $N >= 1$, and $v = v_p (det T')$. With $m_p$ as in [CS88, §12],
#tagged("D.12", $ 2 m_p (Y) dot |O(D_p (T'))| = p^(v (N - 1) slash 2) dot 2 m_p (L_p) dot alpha_"prim" (T', L_p), $)
where, with $Q(x) = (x, x) slash 2$,
$ alpha_"prim" (T', L_p) = epsilon_p dot p^(-(d n - n (n + 1) slash 2)) dot \#{"isometric injections" (T' slash p T', Q) -> (L_p slash p L_p, Q)}, $
with $epsilon_p = 1$ for odd $p$ and $epsilon_p = 2^n$ for $p = 2$.
]
_Strategy of the proof._ Both sides of (D.12) are local volumes. The primitive embeddings of $T'$ into $L_p$ form finitely many orbits of $O(L_p)$, one for each anti-isometry $D_p (T') -> D_p (Y)$ modulo the image of $O(Y)$. Computing the volume of the set of primitive embeddings in two ways, by these orbits and by counting solutions modulo $p$ and lifting them, gives (D.12).
#proofbox[
_Measures._ For a symmetric matrix $S$ of size $k$ let $omega_S = dif X slash dif S$ be the measure on $O(S)$ given by the fibres of $X |-> X^t S X$, with $dif S$ normalized by $"vol"("Sym"_k (ZZ_p)) = 1$. By the definition of the density, $"vol"_(omega_S) O(Lambda) = alpha_p (Lambda)$, and [CS88, §12] defines $m_p (Lambda) = p^(v_p (det Lambda)(k + 1) slash 2) slash alpha_p (Lambda)$. If $S' = C^t S C$ then $dif S' = |det C|^(k + 1) dif S$, so $omega^0 = |det S|^((k + 1) slash 2) omega_S$ is independent of the basis and $m_p (Lambda) = 1 slash "vol"_(omega^0) O(Lambda)$.

_Decomposition._ Let $X = {A : A^t S A = S_(T')}$ with $omega_X = dif A slash dif S_(T')$ and $G_Y = O(Y times.o QQ_p)$. In a basis adapted to $S_0 = S_(T') perp S_Y$ write $g = (A, B)$; near $A = (I, 0)^t$, with $B = (B_1, B_2)^t$, $A^t S_0 B = S_(T') B_1$, so $dif(A^t S_0 B) = |det S_(T')|^N dif B_1$ and $omega_(S_0) = omega_(X, S_0) and |det S_(T')|^(-N) omega_(S_Y)$. Returning to a basis of $L_p$, $S_0 = C^t S C$ with $|det C| = |det S_(T')| = |det S_Y|$ (as $L_p$ is unimodular and $T'$ primitive), $omega_S = |det C|^(d + 1) omega_(S_0)$ and $omega_(X, S) = |det C|^n omega_(X, S_0)$. Collecting the factors, $omega^0_(O(L)) = |det S_(T')|^((1 - N) slash 2) omega_X and omega^0_(G_Y)$.

_Orbits._ Let $X_"prim" subset X$ be the primitive embeddings. The complement of each is isometric to $Y$, since it has the same rank, determinant class and discriminant form, which determine a $p$-adic (even, if $p = 2$) lattice [Nik79, Cor. 1.9.3]. The $O(L_p)$-orbits on $X_"prim"$ correspond to the classes of anti-isometries $D_p (T') -> D_p (Y)$ modulo the image of $O(Y)$, and the stabilizer of an embedding is the kernel of $O(Y) -> O(D_p (Y))$. Hence
$ "vol"_(omega_X) (X_"prim") = |det S_(T')|^((N - 1) slash 2) |O(D_p (T'))| ("vol"^0 O(L_p)) / ("vol"^0 O(Y)) = p^(-v (N - 1) slash 2) |O(D_p (T'))| (m_p (Y)) / (m_p (L_p)). $
The left side is the primitive density.

_Lifting._ For odd $p$, $B |-> A^t S B + B^t S A$ maps onto all symmetric matrices when $A$ is primitive and $S$ unimodular, and Hensel's lemma gives the count modulo $p$. For $p = 2$ the image is the set $"Sym"^e$ of even symmetric matrices, whose volume for $dif S$ is $2^(-n)$; measuring with $"vol"("Sym"^e) = 1$ instead multiplies by $2^n$. In that normalization the solutions are counted modulo $2^r$ (diagonal modulo $2^(r + 1)$); for $r = 1$ they are exactly the injections preserving $Q$, and since $B |-> A^t S B + B^t S A$ maps onto $"Sym"^e$, each lifts by the same factor $2^(d n - n (n + 1) slash 2)$ at every level. #qed
]

For a unimodular $ZZ_p$-lattice $L_p$ of rank $d$ (even if $p = 2$) put $epsilon = ((-1)^(d slash 2) det L_p slash p)$ (Kronecker symbol) when $d$ is even; for $p = 2$, $epsilon = 1$ exactly when $L_2$ is an orthogonal sum of hyperbolic planes. Put
$ hat(s)(p, d, epsilon) = product_(i = 1)^(ceil(d slash 2) - 1) (1 - p^(-2 i))^(-1) dot cases((1 - epsilon p^(-d slash 2))^(-1) & (d "even"), 1 & (d "odd")), $
so that $hat(s)(p, d, epsilon) = 2 p^(dim "GO"_d) slash |"GO"_d^epsilon (FF_p)|$ for odd $p$, with $dim "GO"_d = d (d - 1) slash 2$; this is twice the standard $p$-mass of [CS88, (6)–(8)]. Formally, for $n = d$ (so $T' = L_p$, outside the hypothesis $N >= 1$) the formula gives $2 m_p (L_p) = hat(s)(p, d, epsilon)$ for odd $p$ and $hat(s)(2, d, epsilon) slash 2^d$ for $p = 2$; and for $N = 1$, $2 m_p (chevron.l a chevron.r) = 1$ ($p$ odd) or $1 slash 2$ ($p = 2$), in agreement with [CS88, §6]. The numbers of isometric injections are given by Witt's theorem in terms of the radical of $T' slash p T'$ and the type of its non-degenerate part; for $p = 2$ the invariants are the dimension of the radical, whether $Q$ vanishes on it, and if so the Arf invariant of the non-degenerate part.

*The sums over $C$ at $p = 2$.* As in §C.6, the distribution of the isotropic subgroups $C_2 subset D_2 (T)$ by the invariants entering (D.12) is computed by adding the irreducible components of $T$ one at a time. Write $D_2 (T) = E perp C_"cyc"$, where $E$ is 2-elementary and $C_"cyc"$ is a sum of cyclic blocks $ZZ slash 2^e$, $e = 2, 3, 4$, each with $q = u slash 2^e$ ($u$ odd, modulo $2^(e + 1)$). The 2-elementary part is recorded by $(rho, delta, sigma)$: $rho = dim_(FF_2) E$, $delta$ as in §2.3, and $sigma in ZZ slash 8 ZZ$ its Brown invariant, additive under orthogonal sums, with $sigma = 1, 7, 0, 4$ for $chevron.l 1 slash 2 chevron.r$, $chevron.l 3 slash 2 chevron.r$, the hyperbolic plane $q(x, y) = x y$ and the anisotropic plane $q(x, y) = x^2 + x y + y^2$. For $delta = 0$, $sigma in {0, 4}$ is the type ($+1$ or $-1$) of the quadratic space; for $delta = 1$, $E tilde.equiv chevron.l 1 slash 2 chevron.r^a perp chevron.l 3 slash 2 chevron.r^b$ with $a + b = rho$ and $a - b equiv sigma$ (mod 8). The isotropic subgroups of $E$ are counted in closed form ($delta = 0$: totally singular subspaces of a quadratic space; $delta = 1$: split by the subspace $D_0$ where $q$ is integral and the characteristic vector $c$); the cyclic blocks are then added one at a time by a dynamic programme on the graph of a map from a subgroup of the new block, counting images by orbits of the block automorphisms. The state of the programme is the pair $((rho, delta, sigma), C_"cyc")$ of the residual form $C^perp slash C$, with integer weights indexed by $log_2 |C|$. At the end only three invariants of $T' = T_C$ enter (D.12) at $p = 2$: the number $rho'$ of generators of $D_2 (T')$ (the dimension of the radical of $(T' slash 2 T', Q)$), $delta' = delta$ (whether $Q$ vanishes on that radical), and, if $delta' = 0$, the type $eta = ((w slash 2)) s_E product epsilon(b)$ of the non-degenerate part, where $w$ is the odd part of $det T$, $s_E = plus.minus 1$ the type of $E$, $epsilon(b) = ((u slash 2))$ for a cyclic block and $plus.minus 1$ for a hyperbolic or anisotropic plane, and $((dot slash 2))$ is the Kronecker symbol. The programme agrees with brute-force enumeration for 42 small discriminant forms and for the 515 forms with $|D_2| <= 2^9$ that occur #calc("C17").

#include "integer2-orbits.typ"

#include "exponent9.typ"

#thm[
*Theorem D.11.* $mu_F^emptyset (20, 6) = "3,275,143" slash "7,464,960" approx 0.43874$ for $p_0 = 3$, and $mu_F^emptyset (16, 4) = 1 slash "288,000"$ for $p_0 = 5$.
]

In the first case 811 of the $"10,632"$ root types have non-zero mass and none is negative; in the second 329 #calc("C17"). The 24 classes $F$ of the genus $(20, 6)$ listed in Appendix B.4 are pairwise non-isometric, rootless, and the sum of $1 slash |"Aut" F|$ over them equals Theorem D.11; so the list is complete. For $(16, 4)$, $p_0 = 5$, the value is attained by a single lattice with $|"Aut" F| = "288,000"$ (the fixed lattice of an element of frame shape $5^4 1^12$, that is of type $5 - (4, 16) - 4$, of $Lambda'_32$, Table A1), which is therefore the only rootless lattice of the genus.

#ptable(
  columns: 2, align: (left, left),
  caption: [Table D5. Checks of (D.10)–(D.12) #calc("C17"). They are consistency checks; the values rest on the proofs above.],
  table.header([check], [result]),
  [$2 m_p$ against [CS88, Table 2] ($p = 3, 5, 7$, $n <= 5$) and the standard mass against [CS88, Tables 3, 6]], [equal],
  [(D.12) against Sage's implementation of the local factors of [CS88], for the 40 parabolic sublattices of $E_8^3$ (including type I and bound components at 2)], [equal],
  [$p_0 = 3$: masses by root system against the complete lists of the genera $(12, 6)$, $(12, 8)$, $(12, 10)$, $(14, 7)$, $(14, 9)$, $(16, 8)$, all classes with roots included], [equal for every root system (e.g. $mu_F^emptyset (16, 8) = 59 slash "1,036,800"$, 79 root systems)],
  [$p_0 = 3$: $(8, 0)$ and $(16, 0)$], [$mu_F^emptyset = 0$ ($E_8$; $E_8^2$, $D_(16)^+$)],
  [$p_0 = 5$: masses by root system against the complete lists of the genera of rank below 24, as for $p_0 = 3$; in rank 24, the rootless mass against a value computed directly], [equal],
)

#rem[
*Remark (vectors of norm $2 slash 3$ in $F^\#$).* For the type $3 - (6, 20) - 6$ the extremal gluings need $F^\#$ free of vectors of norm $2 slash 3$ (Lemma 4.3). For a rootless $F$ and $v in F^\#$ of norm $2 slash 3$, $r = 3 v in F$ is primitive of norm 6 with $(r, F) = 3 ZZ$, and $(F, plus.minus v) <-> Y = r^perp$ is a bijection onto the rootless lattices of the genus $cal(K)_(19)$ of even lattices of rank 19 and determinant $2 dot 3^5$ with discriminant form $chevron.l 1 slash 2 chevron.r perp q'$, where $q'$ is the restriction of $q_F$ to the orthogonal complement of a vector of norm $2 slash 3$ in $D(F) tilde.equiv FF_3^6$ (well defined by Witt's theorem). Indeed $[F : ZZ r perp Y] = 2$, glued along $r slash 2$ with $q = 3 slash 2$, and at 3 the summand $ZZ_3 r$ splits off; the genus is determined by the discriminant form [Nik79, Cor. 1.9.4]. Moreover $|"Aut"(F, plus.minus v)| = 2 |"Aut" Y|$, the factor 2 being the reflection in $r$, which lies in $"Aut" F$ since $(r, F) = 3 ZZ$ and $(r, r) = 6$; the same method gives $mu^emptyset (cal(K)_(19)) = 1 slash 648$ #calc("C17"), so the lattices $F$ with such vectors carry total weight $sum n(F) slash |"Aut" F| = 1 slash "1,296"$, where $n(F)$ is the number of pairs $plus.minus v$ of vectors of norm $2 slash 3$ in $F^\#$. In the list exactly two classes have such vectors, each with one pair $plus.minus v$ and $|"Aut" F| = "2,592"$, which accounts for $1 slash "1,296"$. This is a check; the list is complete by Theorem D.11.
]


= References <refs>
#refs-list(
  [[Aux] roiban1344, _An element-assignment computation of the number of extremal even unimodular lattices of rank 32_, auxiliary note, 2026, distributed with this paper.],
  [[BB92] C. Bachoc, C. Batut, _Étude algorithmique de réseaux construits avec la forme trace_, Experiment. Math. 1 (1992) 183–190.],
  [[BCDYY25] P. Boyvalenkov, D. Cherkashin, P. Dragnev, D. Yorgov, V. Yorgov, _On $T$-avoiding spherical codes and designs in $RR^32$_, arXiv:2512.23092.],
  [[BN97] C. Bachoc, G. Nebe, _Classification of two genera of 32-dimensional lattices of rank 8 over the Hurwitz order_, Experiment. Math. 6 (1997) 151–162.],
  [[BV01] C. Bachoc, B. Venkov, _Modular forms, lattices and spherical designs_, in: J. Martinet (ed.), Réseaux euclidiens, designs sphériques et formes modulaires, Monogr. Enseign. Math. 37, L'Enseignement Mathématique, Genève, 2001, 87–111.],
  [[CE03] H. Cohn, N. Elkies, _New upper bounds on sphere packings I_, Ann. of Math. 157 (2003) 689–714.],
  [[ChT26] G. Chenevier, O. Taïbi, _Unimodular lattices of rank 29 and related even genera of small determinant_, arXiv:2601.19780.],
  [[CS88] J. H. Conway, N. J. A. Sloane, _Low-dimensional lattices IV: the mass formula_, Proc. R. Soc. Lond. A 419 (1988) 259–286.],
  [[D] roiban1344, _The extremal even unimodular lattices of rank 32 with $|"Aut" L| > 2$_, database, version 2 (2026), distributed with this paper: the release asset `database-v2.tar.gz`, one record per lattice of $S$ with its label of §6.8, a Gram matrix and generators of $"Aut" L$, together with the file `data/classified-lattices.jsonl` of the repository, which gives for each label the order, the part of Theorem B and the group.],
  [[GAP] The GAP Group, _GAP – Groups, Algorithms, and Programming_, Version 4.15.1, 2025, https://www.gap-system.org.],
  [[GHY01] W. T. Gan, J. P. Hanke, J.-K. Yu, _On an exact mass formula of Shimura_, Duke Math. J. 107 (2001) 101–131.],
  [[Gor80] D. Gorenstein, _Finite Groups_, 2nd ed., AMS Chelsea Publishing (first edition 1968; second edition 1980).],
  [[Hen13] A. Henn, M. Hentschel, A. Krieg, G. Nebe, _On the classification of lattices over $QQ(sqrt(-3))$ which are even unimodular $ZZ$-lattices of rank 32_, Int. J. Math. Math. Sci. 2013, Art. 837080.],
  [[HK89] K. Hashimoto, H. Koseki, _Class numbers of definite unimodular hermitian forms over the rings of imaginary quadratic fields_, Tohoku Math. J. (2) 41 (1989) 1–30.],
  [[HKN09] M. Hentschel, A. Krieg, G. Nebe, _On the classification of lattices over $QQ(sqrt(-3))$ which are even unimodular $ZZ$-lattices_, Abh. Math. Semin. Univ. Hambg. 80 (2010) 183–192; arXiv:0903.4334.],
  [[HKN12] M. Hentschel, A. Krieg, G. Nebe, _On the classification of even unimodular lattices with a complex structure_, Int. J. Number Theory 8 (2012) 983–992.],
  [[JK07] T. Junttila, P. Kaski, _Engineering an efficient canonical labeling tool for large and sparse graphs_, in: Proc. ALENEX 2007, SIAM, 135–149.],
  [[JK11] T. Junttila, P. Kaski, _Conflict propagation and component recursion for canonical labeling_, in: Proc. TAPAS 2011, Lecture Notes in Comput. Sci. 6595, Springer, 151–162.],
  [[Jür15] M. Jürgens, _Nicht-Existenz und Konstruktion extremaler Gitter_, Dissertation, TU Dortmund, 2015.],
  [[Kin03] O. King, _A mass formula for unimodular lattices with no roots_, Math. Comp. 72 (2003) 839–863; arXiv:math/0012231.],
  [[KN93] H. Koch, G. Nebe, _Extremal even unimodular lattices of rank 32 and related codes_, Math. Nachr. 161 (1993) 309–319.],
  [[KV89] H. Koch, B. B. Venkov, _Über ganzzahlige unimodulare euklidische Gitter_, J. reine angew. Math. 398 (1989) 144–168.],
  [[KV91] H. Koch, B. B. Venkov, _Über gerade unimodulare Gitter der Dimension 32, III_, Math. Nachr. 152 (1991) 191–213.],
  [[MM76] J. M. Masley, H. L. Montgomery, _Cyclotomic fields with unique factorization_, J. Reine Angew. Math. 286/287 (1976) 248–256 (cited in [Aux]).],
  [[MOS75] C. L. Mallows, A. M. Odlyzko, N. J. A. Sloane, _Upper bounds for modular forms, lattices, and codes_, J. Algebra 36 (1975) 68–76.],
  [[Neb13b] G. Nebe, _On automorphisms of extremal even unimodular lattices_, Int. J. Number Theory 9 (2013) 1933–1959; arXiv:1212.0865.],
  [[Neb98q] G. Nebe, _Finite quaternionic matrix groups_, Represent. Theory 2 (1998) 106–223 (cited from the author's preprint; section and theorem numbers as there).],
  [[Nik79] V. V. Nikulin, _Integral symmetric bilinear forms and some of their applications_, Izv. Akad. Nauk SSSR Ser. Mat. 43 (1979) 111–177; English transl. Math. USSR-Izv. 14 (1980) 103–167.],
  [[NS] G. Nebe, N. J. A. Sloane, _A Catalogue of Lattices_, https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/.],
  [[PARI] The PARI Group, _PARI/GP_, version 2.17.3, Univ. Bordeaux, 2025, http://pari.math.u-bordeaux.fr/.],
  [[Rec] roiban1344, _Research record for “The extremal even unimodular lattices of rank 32”_: programs, runs and verification notes, 2026, release assets `research-record.tar.gz` (programs and notes) and `research-runs.tar.gz` (runs) of the repository of §7.4.],
  [[Sage] The Sage Developers, _SageMath, the Sage Mathematics Software System_, Version 10.9, 2026, https://www.sagemath.org.],
  [[Sch98] R. Schoof, _Minus class groups of the fields of the $ell$-th roots of unity_, Math. Comp. 67 (1998) 1225–1245.],
  [[Shi08] G. Shimura, _Arithmetic of Hermitian forms_, Doc. Math. 13 (2008) 739–774.],
  [[SPLAG] J. H. Conway, N. J. A. Sloane, _Sphere Packings, Lattices and Groups_, 3rd ed., Springer, 1999.],
  [[Wei65] A. Weil, _Sur la formule de Siegel dans la théorie des groupes classiques_, Acta Math. 113 (1965) 1–87.],
)
For [HKN09], [HKN12], [Kin03] and [Neb13b] we read the arXiv or author's versions; theorem numbers refer to them.
