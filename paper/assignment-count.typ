#import "style.typ": *
#import "computations.typ": *
#show: paper.with(
  lang: "en", title: [An element-assignment computation of the number of extremal even unimodular lattices of rank 32],
  subtitle: [Auxiliary note to “The extremal even unimodular lattices of rank 32”],
  authors: [roiban1344], date: [Version 2026-11-14.2],
)

#show figure.where(kind: table): set block(breakable: true)

#policy-box[
This note records an organization of the count by assigned group elements: every automorphism other than $plus.minus 1$ is assigned to one element of a short list, and the count becomes a sum of thirteen terms involving the numbers $e(h)$ and the masses $M_2, dots, M_5$. It was used as a check of the classification-based proof of the paper and is not a premise of that proof. It keeps the arrangement of an earlier version of the paper: Sections 3–5, Tables 5–7 and Appendix C.8 of this note are its own, and its references to the paper use that arrangement. There, Lemma 2.4, Propositions 3.4 and 3.5, Lemmas 4.3, 4.4 and 4.5, Theorem 4.6 and Proposition 4.7 are Lemma 2.2, Propositions 3.2 and 3.3, Lemmas 4.2, 4.3 and 4.4, Theorem 4.5 and Proposition 4.6 of the paper, and Tables 4 and 4a are Tables 5 and 6. Its sections keep their numbers in that arrangement, so the numbering has gaps where sections are now in the paper (§4.1.1 follows §3.2), and a reference to a section that this note does not contain says "of the paper". The table of its Appendix C.8 is called Table C8 here so as not to be confused with Table C2 of the paper. Definitions and references are those of the paper.
]

= 3. The counting formula <s3>

#policy-box[
*What this section does.* We rewrite $N$ as a sum over automorphisms (§3.1), assign every automorphism $g != plus.minus 1$ to exactly one element of a small list (Proposition 3.2), and obtain formula (3.3), which involves no individual $|"Aut" L|$. We then list the possible types of the assigned elements (§3.3 of the paper and its Table 4). Everything in this section is elementary, except that Proposition 3.5 uses the exact linear programs and lists of #calc("C2")\; its proof also uses Lemma 4.5 from §4.2 of the paper, which does not depend on §3.
]

== 3.1 Summing over automorphisms
A _pair_ $(L, h)$ consists of $L in cal(L)$ and $h in "Aut" L$; pairs $(L, h)$ and $(L', h')$ are isometric if there is an isometry $u : L -> L'$ with $u h u^(-1) = h'$. For fixed $L$, the isometry classes of pairs $(L, h)$ are the conjugacy classes of $"Aut" L$, and the class of $h$ has $|"Aut" L| slash |C(h)|$ elements. Since $g mapsto 1$ summed over $"Aut" L$ gives $|"Aut" L|$, and $g = 1$, $g = -1$ contribute $2 slash |"Aut" L|$ to each $L$,
#tagged("3.1", $ N = sum_(L in cal(L)) 1 / (|"Aut" L|) sum_(g in "Aut" L) 1 = 2 m + sum_(L in cal(L)) (\#{g in "Aut" L : g != plus.minus 1}) / (|"Aut" L|) $)
and, since $|"Aut" L| = 2$ exactly when $L in.not S$,
#tagged("3.2", $ N = 2 m + sum_(L in S) (1 - 2 / (|"Aut" L|)). $)
As $|"Aut" L|$ is even, each term of (3.2) lies in $[1 slash 2, 1)$, so $2 m + |S| slash 2 <= N < 2 m + |S|$. Formula (3.2) needs the individual $|"Aut" L|$; formula (3.3) below does not.

== 3.2 Assignment of automorphisms
#thm[
*Proposition 3.2.* Let $g in "Aut" L$, $g != plus.minus 1$, of order $n$. Exactly one of the following holds.
+ $n$ is divisible by an odd prime; let $p$ be the smallest one; then $h = g^(n slash p)$ has order $p$ and $g in C(h)$.
+ $n = 2$: $g$ is an involution $sigma != plus.minus 1$.
+ $n = 2^a$ with $a >= 2$ and $sigma = g^(2^(a-1)) != -1$: $sigma$ is an involution $!= plus.minus 1$ and $g in C(sigma)$.
+ $n = 2^a$ with $a >= 2$ and $g^(2^(a-1)) = -1$: the minimal polynomial of $g$ is $Phi_(2^a)$, and $L$ is a Hermitian $ZZ[zeta_(2^a)]$-lattice of rank $32 slash 2^(a-1)$; moreover $a <= 5$.
]
#proofbox[
Since $g != 1$, $n > 1$. If $n$ has an odd prime factor we are in case 1, and $h = g^(n slash p)$ has order $p$ and commutes with $g$. Otherwise $n = 2^a$ with $a >= 1$. If $a = 1$, $g$ is an involution, different from $plus.minus 1$ by assumption: case 2. If $a >= 2$, $sigma = g^(2^(a-1))$ has order 2. Either $sigma = -1$ (case 4) or $sigma != -1$; in the latter case $sigma != 1$ as well, so $sigma$ is an involution $!= plus.minus 1$ commuting with $g$ (case 3). The cases are distinguished by $n$ and by whether $g^(n slash 2) = -1$, so exactly one holds. In case 4, $g$ is a root of $X^(2^(a-1)) + 1 = Phi_(2^a)$, which is irreducible over $QQ$, so it is the minimal polynomial of $g$; then $a <= 5$ by Lemma 2.4 and the remark after it. #qed
]

We say that $g$ is _assigned_ to $h$ in case 1, to $sigma = g$ in case 2 and to $sigma$ in case 3. An element assigned to $h$ lies in $C(h)$, and assignment commutes with conjugation: if $g$ is assigned to $h$, then $u g u^(-1)$ is assigned to $u h u^(-1)$ for every $u in "Aut" L$. For $h$ of odd prime order, or an involution $h != plus.minus 1$, let $e(h)$ be the number of elements of $"Aut" L$ assigned to $h$; it depends only on the isometry class of $(L, h)$, and $e(h) >= 1$ since $h$ is assigned to itself. For $a = 2, dots, 5$ let
#tagged("3.3a", $ M_a = sum_((L, g)) 1 / (|C(g)|), $)
the sum over isometry classes of pairs with $L in cal(L)$ and $g^(2^(a-1)) = -1$. Each $M_a$ is a mass of rootless Hermitian lattices over $ZZ[zeta_(2^a)]$: for $a = 2$ by §4.3 of the paper, for $a = 3, 4, 5$ by Appendix C.8 of this note.

#thm[
*Corollary 3.3.*
#tagged("3.3", $ N = 2 m + sum_((L, h)) e(h) / (|C(h)|) + sum_(a = 2)^5 M_a, $)
where $(L, h)$ runs over isometry classes of pairs with $L in cal(L)$ and $h$ of odd prime order or an involution $!= plus.minus 1$.
]
#proofbox[
Fix $L$ and put $G = "Aut" L$. By Proposition 3.2, $G without {plus.minus 1}$ is the disjoint union of the sets $A(h)$ of elements assigned to $h$, over the elements $h$ of odd prime order and the involutions $h != plus.minus 1$, and of the sets $B_a = {g : g^(2^(a-1)) = -1}$ for $2 <= a <= 5$. Grouping the $h$ into conjugacy classes, a class of size $|G| slash |C(h)|$ contributes $(|G| slash |C(h)|) dot e(h)$ elements; grouping $B_a$ into conjugacy classes gives $|B_a| = sum_([g]) |G| slash |C(g)|$. Divide by $|G|$, sum over $L in cal(L)$ and insert in (3.1). #qed
]

No individual $|"Aut" L|$ occurs in (3.3): only the centralizers of the elements $h$ and $g$, and the numbers $e(h)$, which are computed inside $C(h)$.



=== 4.1.1 Orientation

For an involution, $sigma$ and $-sigma$ are different elements, and the pair $(L, -sigma)$ has fixed lattice $K$. If $F tilde.not K$, the pair $(L, -sigma)$ appears in the list of type $(32 - d, d)$ as the class of $phi^(-1)$, not in the list of type $(d, 32 - d)$; we compute both contributions on the list of $(d, 32 - d)$, and the contribution of the type is $E = E(sigma) + E(-sigma)$ with $E(plus.minus sigma) = sum_([phi]) e(plus.minus sigma) slash |"Stab"(phi)|$. If $F tilde.equiv K$ (type $(16, 16)$ with equal components), the classes $[phi]$ and $[phi^(-1)]$ in the same list correspond to $(L, sigma)$ and $(L, -sigma)$, and the contribution is $E = E(sigma)$ summed over all classes of the list; the two coincide exactly when $sigma tilde -sigma$. In both cases $E$ is the sum of $e(sigma) slash |"Stab"(phi)|$ over the extremal classes $[phi]$ of the list, together with that of $e(-sigma) slash |"Stab"(phi)|$ when $F tilde.not K$.

#thm[
*Lemma 4.2.* Let $sigma != plus.minus 1$ be an involution and $g$ an element of order $2^a$, $a >= 2$, with $g^(2^(a-1)) = sigma$. Then $2^(a-1)$ divides $"rank" K$. In particular $e(sigma) = 1$ if $"rank" K$ is odd.
]
#proofbox[
$g$ commutes with $sigma$, so it preserves $K times.o QQ$, on which $g^(2^(a-1)) = -1$. As in Proposition 3.2, the minimal polynomial of $g|_K$ is the irreducible polynomial $Phi_(2^a)$ of degree $2^(a-1)$, so $2^(a-1)$ divides $"rank" K$. Elements assigned to $sigma$ other than $sigma$ itself are of this form (case 3 of Proposition 3.2). If $"rank" K$ is odd, no such $g$ exists, since $2^(a-1) >= 2$; then $sigma$ is the only element assigned to $sigma$, and $e(sigma) = 1$. #qed
]

For types $(15, 17)$ and $(17, 15)$ both ranks are odd, so $e = 1$ and the contribution is twice the extremal mass of the list of type $(15, 17)$.



*The numbers $e(sigma)$.* If $g$ is assigned to $sigma$ in Proposition 3.2, then $sigma$ is a power of $g$, so $g$ commutes with $sigma$; hence $e(sigma)$ is the number of $g in C(sigma)$ with $g = sigma$ or $g^(2^(a-1)) = sigma$ for some $a >= 2$. For each extremal class, generators of $"Stab"(phi) tilde.equiv C(sigma)$ are lifted to automorphisms of $L_phi$ (checked in integers), and $e(sigma)$ and $e(-sigma)$ are counted in this group with GAP #calc("C34"). By Lemma 4.2 only $a$ with $2^(a-1) | 32 - d$ occur for $sigma$, and only $a$ with $2^(a-1) | d$ for $-sigma$.



*The numbers $e(h)$.* An element $g$ is assigned to $h$ when the smallest odd prime dividing its order $n$ is $p$ and $g^(n slash p) = h$; then $g in C(h)$. If $p^2 divides.not |C(h)|$, the Sylow $p$-subgroup $chevron.l h chevron.r$ is central in $C(h)$, so $C(h) = chevron.l h chevron.r times H$ by the Schur–Zassenhaus theorem, and writing $g = h^a t$ with $t in H$ of order $m$ one finds that $g$ is assigned to $h$ if and only if $m$ has no odd prime factor smaller than $p$ and $a m equiv 1 med (mod p)$. Hence
$ e(h) = \#{t in H : "the order of" t "has no odd prime factor" < p}, $
and, under the same hypothesis $p^2 divides.not |C(h)|$, $e(h) = |C(h)| slash 3$ for $p = 3$. In general $e(h)$ was counted in $C(h)$, built from generators of $"Stab"(phi)$ lifted to $L_phi$ and checked in integers, by conjugacy classes in GAP #calc("C8", "C9")\; in all $1234 + 204$ cases where the group was built, its order equals $|"Stab"(phi)|$.



#ptable(
  columns: 3, align: (left, right, right),
  caption: [Table 5. Elements of odd prime order: the contribution $E = sum e(h) slash |C(h)|$ of each type (exact) #calc("C8", "C9").],
  table.header([type], [$E$], [approx.]),
  [$3 - (9, 14) - 7$], [$3302 slash 6561$], [$0.503$],
  [$3 - (10, 12) - 6$], [$10 slash 81$], [$0.123$],
  [$3 - (10, 12) - 8$], [$7 slash 3$], [$2.333$],
  [$3 - (8, 16) - 8$], [$1438 slash 81$], [$17.753$],
  [$3 - (12, 8) - 8$], [$4187858 slash 2187$], [$1914.887$],
  [$3 - (9, 14) - 9$], [$28961 slash 729$], [$39.727$],
  [$3 - (10, 12) - 10$], [$9640 slash 81$], [$119.012$],
  [$3 - (13, 6) - 5$], [$3676 slash 243$], [$15.128$],
  [$3 - (6, 20) - 6$], [$10048 slash 2187$], [$4.594$],
  [$3 - (16, 0) - 0$], [$413 slash 243$], [$1.700$],
  [*order 3*], [*$13881509 slash 6561$*], [*$2115.761$*],
  [$5 - (4, 16) - 4$], [$38 slash 1875$], [$0.020$],
  [$5 - (5, 12) - 5$], [$214 slash 625$], [$0.342$],
  [$5 - (6, 8) - 4$], [$296 slash 625$], [$0.474$],
  [$5 - (6, 8) - 6$], [$442 slash 15$], [$29.467$],
  [$5 - (7, 4) - 3$], [$836 slash 625$], [$1.338$],
  [$5 - (8, 0) - 0$], [$8192 slash 16875$], [$0.485$],
  [$7 - (4, 8) - 4$], [$134 slash 63$], [$2.127$],
  [$17 - (2, 0) - 0$], [$13 slash 51$], [$0.255$],
  [$31 - (1, 2) - 1$], [$8 slash 31$], [$0.258$],
  [*orders 5, 7, 17, 31*], [*$2164244689 slash 62251875$*], [*$34.766$*],
)



The enumeration also confirms Theorem 4.6. Its sum is a lower bound for $M_2$ whatever the formula says, and equality could hide an error only if an error in the formula were exactly compensated by missing lattices of the same total weight. For $a = 4$ an independent list of 514 classes, of which 238 are rootless, has total weight $M_4$. For $a = 3$ (rank 8) no independent list exists; the value rests on the proof of the formula and on its agreement at ranks 2, 4 and 6 (Table 6).

#ptable(
  columns: 4, align: (left, right, right, left),
  caption: [Table 6. Rootless Hermitian masses $M_a$.],
  table.header([$a$ (rank)], [$M_a$], [approx.], [confirmed by]),
  [2 (16)], [$4128748922231417 slash 120324096000$], [$34313.567$], [enumeration of $76742$ classes, exact #calc("C12")],
  [3 (8)], [$1373459063 slash 7741440$], [$177.416$], [no independent check at rank 8; the same formula agrees at ranks 2, 4, 6, and at rank 6 the rootless value $1 slash 48384$ agrees with a centralizer order in $2."Co"_1$ #calc("C35")],
  [4 (4)], [$6983 slash 640$], [$10.911$], [list of 514 classes #calc("C35")],
  [5 (2)], [$31 slash 24$], [$1.292$], [explicit count of the 16 types of rank 2 (Table C8)],
)



== 4.4 The thirteen terms

#ptable(
  columns: 4, align: (left, right, right, left),
  caption: [Table 7. The thirteen terms of (3.3) #calc("C36"). The term $2 m$ is twice (2.1).],
  table.header([term], [exact value], [approx.], [§]),
  [$2 m$], [$1310037331282023326658917 slash 119431715880960000$], [$10968923.302$], [2.2],
  [$M_2$], [$4128748922231417 slash 120324096000$], [$34313.567$], [4.3],
  [$(16, 16)$, $s = 16$, $delta = 1$], [$2880305926221 slash 293601280$], [$9810.264$], [4.1],
  [$(15, 17)$], [$20874607 slash 2160$], [$9664.170$], [4.1],
  [$(14, 18)$], [$491776351403 slash 70778880$], [$6948.066$], [4.1],
  [$(8, 24)$], [$4903343345383327 slash 792723456000$], [$6185.440$], [4.1],
  [$(12, 20)$], [$24698881091413 slash 4246732800$], [$5815.972$], [4.1],
  [order 3], [$13881509 slash 6561$], [$2115.761$], [4.2],
  [$(16, 16)$, other 10 pairs], [$1524348533472859 slash 5707608883200$], [$267.073$], [4.1],
  [$M_3$], [$1373459063 slash 7741440$], [$177.416$], [C.8],
  [orders 5, 7, 17, 31], [$2164244689 slash 62251875$], [$34.766$], [4.2],
  [$M_4$], [$6983 slash 640$], [$10.911$], [C.8],
  [$M_5$], [$31 slash 24$], [$1.292$], [C.8],
  [*sum*], [*$11044268$*], [], [5],
)

== 4.5 One term in full: $31 - (1, 2) - 1$
Every term of Table 7 for an element of odd prime order or an involution is made in the same way: list the components, list their gluings, keep the extremal ones, and count the elements assigned to $h$. For most types the lists are long and the steps are done by programs. For the type $31 - (1, 2) - 1$ every step can be done by hand, and we do so here; the result is the entry $8 slash 31$ of Table 5.

Let $h$ be an automorphism of order 31 of $L in cal(L)$ of type $31 - (1, 2) - 1$. Then $F = L^h$ has rank 2 and $F^\# slash F tilde.equiv FF_31$, and $Z = F^perp inter L$ has rank 30 and is a Hermitian $ZZ[zeta_31]$-lattice of rank 1 with $Z^\# slash Z tilde.equiv FF_31$ (§3.3 of the paper).

*The lattice $F$.* $F$ is a positive definite even binary lattice of determinant 31, with Gram matrix $mat(2 a, b; b, 2 c)$ and $b^2 - 4 a c = -31$. By the reduction of binary forms, the reduced forms $a x^2 + b x y + c y^2$ with $|b| <= a <= c$ are $(a, b, c) = (1, 1, 8)$ and $(2, plus.minus 1, 4)$, and the last two give isometric lattices. Hence the genus consists of two lattices, $mat(2, 1; 1, 16)$ of minimum 2 and
$ F = mat(4, 1; 1, 8), $
which is the only one of minimum at least 4. Since $0 < b < a < c$, $"Aut" F = {plus.minus 1}$. The dual $F^\#$ has Gram matrix $frac(1, 31) mat(8, -1; -1, 4)$, so $min F^\# = 4 slash 31$; in particular $F^\#$ has no vector of norm $2 slash 31$. On $F^\# slash F tilde.equiv FF_31$ the discriminant form is non-degenerate of dimension 1, so $O(q_F) = {plus.minus 1}$, and $"Aut" F -> O(q_F)$ is onto.

*The lattices $Z$.* Here $s = z = 1$, so $Z$ is unimodular for the form of Lemma 4.5, and $Z$ runs over the genus $(1, 1)$ of Proposition D.1, the rank-1 Hermitian lattices, that is, ideal classes with their forms. This genus has $h^-(QQ(zeta_31)) = 9$ classes: the minus class group of $QQ(zeta_31)$ is cyclic of order 9, a result of Kummer [Sch98, p. 1226]. Independently, the mass of the genus given by the formula of Appendix D is $9 slash 62$ #calc("C2"), which with $|"U"(Z)| = 62$ for every class (below) again gives nine classes. A unitary automorphism of a lattice of rank 1 is multiplication by a unit $u$ with $u overline(u) = 1$; all conjugates of such a $u$ have absolute value 1, so $u$ is a root of unity, and $"U"(Z) = chevron.l -zeta_31 chevron.r$ has order 62 for every class. By (4.2), a root $x$ of $Z$ spans an orthogonal summand $ZZ[zeta_31] x tilde.equiv A_30$, which for rank 1 is all of $Z$; so exactly one class, $A_30$, has roots. The other eight classes are rootless, with total weight $8 slash 62 = 4 slash 31$, the entry of Table D1; being even and rootless, they have minimum at least 4.

*The gluings.* $"Iso"(q_F, -q_Z)$ has $|O(q_F)| = 2$ elements. The image of $"Aut" F$ in $O(q_F) = {plus.minus 1}$ is everything, and so is the image of $"U"(Z)$, because $zeta_31$ acts trivially on $Z^\# slash Z$ (proof of Proposition 4.1(i)) and $-1$ acts as $-1$. So for each $Z$ there is exactly one double coset, with
$ "Stab"(phi) = {(epsilon, epsilon zeta_31^k) : epsilon = plus.minus 1, thin 0 <= k < 31} $
of order 62, in agreement with (4.1): $1 slash 62 = 2 slash (2 dot 62)$. By Lemma 4.4, since $F^\#$ has no vector of norm $2 slash 31$, every gluing is extremal. This gives eight pairs $(L, h)$, one for each rootless $Z$.

*The centralizer and $e(h)$.* By Proposition 4.1(iii), $C(h) tilde.equiv "Stab"(phi) = chevron.l -1 chevron.r times chevron.l h chevron.r$, cyclic of order 62. Since $31^2 divides.not 62$, the formula of the paragraph _The numbers $e(h)$_ above gives $e(h) = \#{t in {plus.minus 1}} = 2$: the elements assigned to $h$ are $h$ itself and $g = -h^(16)$, of order 62 with $g^2 = h^(32) = h$.

*The contribution.* Each of the eight pairs contributes $e(h) slash |C(h)| = 2 slash 62$, so the type contributes
$ E = 8 dot 2 / 62 = 8 / 31 approx 0.258, $
the entry of Table 5. Integrality of the whole sum (§5.3) would detect an error here: a ninth rootless class would change $E$ by $2 slash 62 = 1 slash 31$, and an error of one in one $e(h)$ would change it by $1 slash 62$.

The underlying lattices of these eight pairs need not be pairwise non-isometric and are not needed for this count. They can nevertheless be constructed. Taking for $Z$ the ideal $P slash overline(P)$ for a prime $P$ above a prime $ell equiv 1 med (mod 31)$, the fourteen primes $ell < 3200$ give lattices $Z$ of minimum 4, and the lattices $L$ obtained have the invariants of $"BW"_32$ or of the code lattice of the extended quadratic residue code of length 32 (Appendix A); this sample is not used in the proof #calc("C9").



== 5.3 What integrality checks
That the sum is an integer is not used in the proof, but it detects many errors. A pair $(L, h)$ contributes $e(h) slash |C(h)|$ with $1 <= e(h) <= |C(h)| - 1$, since $1 in C(h)$ is never assigned to $h$. Hence an error in a single $e(h)$, or a single pair wrongly kept or discarded as extremal, changes the sum by a non-zero number of absolute value less than 1, which makes it non-integral. Errors that are not detected are: two errors that cancel modulo 1; an error by an integer amount in one term, which is possible only for the masses $2 m$, $M_a$ and the totals of large lists; and a uniform error of a whole term by an integer factor. Against these the proof relies on §4 and Appendices B–D, not on integrality.

== 5.4 How each term is checked
Every value in Table 7 is protected in three ways, and Table 7a collects them term by term.

+ _Completeness within a term._ No list in §4 is trusted because of how it was found. Each one is complete by an exact equality of masses (Lemma 7.1): (4.1) for the double cosets of each of the 83 pairs of components of involutions and of each pair of odd prime order; (B.1) for the lattices $K$ of type $(8, 24)$; the mass of the genus, the rootless mass and (B.2) for the components of odd prime order (Table B4). A missing class would leave the left side short of the right side, unless the value on the right were wrong by exactly the missing weight. The components of involutions are complete instead by the exhaustive enumeration of Lemma B.1.
+ _Independent confirmations._ A value is confirmed when a second, independent method gives the same fraction: King's table for $m$, the enumeration of Proposition 4.7 for $M_2$, published classifications for two types, and the checks of Appendices B.2, C.7 and D.3. For the pairs with $e(h) > 1$ of the types $(16, 16)$ with $delta = 1$, $(14, 18)$, $(12, 20)$ and $(8, 24)$, $e(h)$ was also recounted inside $"Aut" L$ computed with bliss, for all pairs with $|C(h)| <= 64$ and for all 152 pairs of type $(16, 16)$, $delta = 1$ with $|C(h)| > 64$, with no disagreement; these recounts use the same counting program and check the lifting of the groups, not the counting; the counting programs were checked in rank 24 against $"Co"_0$ (§5.6), which again checks the programs, not the values in rank 32.
+ _Global identities._ The integer auxiliary sum and the historical element-by-element comparisons test the assignment computation. They share component lists with the classification and do not replace the completeness or full-group arguments. For the main proof, (5.1) is aggregated directly from the classified rows #calc("C13").

#[
#set text(size: 0.86em)
#set par(justify: false)
#ptable(
  columns: (7.2em, 1fr, 1fr, 1.25fr), align: (left, left, left, left),
  caption: [Table 7a. How each term of Table 7 is computed and checked. "Complete by" is the identity that proves its lists complete; "confirmed by" lists independent checks of its value. The global identities are in the third item above the table; what remains assumed is in Table 10 of the paper.],
  table.header([term], [computed by], [complete by], [confirmed by]),
  [$2 m$], [King's inclusion–exclusion #calc("C14")], [— (a formula)], [King's table [Kin03], exactly],
  [$M_2$], [formula of App. C #calc("C11")], [— (a formula)], [enumeration of $76742$ classes (Prop. 4.7) #calc("C12")],
  [$(16, 16)$, $s = 16$, $delta = 1$], [gluings #calc("C4")\; $e$ #calc("C34")], [(4.1)\; components by Lemma B.1], [components by a classification of the genera (§B.2)\; $e$ recounted in $"Aut" L$ (lifting only, §5.4)\; the program for $e$ in rank 24 (§5.6)],
  [$(15, 17)$], [gluings #calc("C4")\; $e = 1$ (Lemma 4.2)], [(4.1)\; components by Lemma B.1], [rootless masses of the genera of the components (§B.2)],
  [$(14, 18)$], [gluings #calc("C4")\; $e$ #calc("C34")], [(4.1)\; components by Lemma B.1], [rootless masses (§B.2)\; $e$ recounted for $|C| <= 64$\; the program for $e$ in rank 24 (§5.6)],
  [$(8, 24)$], [Theorem B.3, Niemeier lattices #calc("C6")\; $e$ #calc("C34")], [(B.1) for Leech, $A_1^(24)$, $A_2^(12)$], [$e$ recounted for $|C| <= 64$\; the program for $e$ in rank 24 (§5.6)],
  [$(12, 20)$], [gluings #calc("C4"), components #calc("C5")\; $e$ #calc("C34")], [(4.1)\; components by Lemma B.1], [the 66 components also by the mass formula #calc("C5")\; rootless mass (§B.2)\; $e$ recounted for $|C| <= 64$\; the program for $e$ in rank 24 (§5.6)],
  [order 3], [gluings and $e$ #calc("C8")], [(4.1)\; components by mass and (B.2)], [$3 - (16, 0) - 0$ by [Hen13]\; the program for $e$ in rank 24 (§5.6)],
  [$(16, 16)$, other 10 pairs], [gluings #calc("C4")\; $e$ #calc("C34")], [(4.1)\; components by Lemma B.1], [components by a classification of the genera (§B.2)\; the program for $e$ in rank 24 (§5.6)],
  [$M_3$], [formula of App. C #calc("C35")], [— (a formula)], [none at rank 8\; the same formula agrees at ranks 2, 4, 6, and with $2."Co"_1$ at rank 6],
  [orders 5, 7, 17, 31], [gluings and $e$ #calc("C9")], [(4.1)\; components by Table D1 and (B.2)], [$17 - (2, 0) - 0$ by [Jür15]\; App. D against the Leech lattice (§D.3)\; the program for $e$ in rank 24 (§5.6)],
  [$M_4$], [formula of App. C #calc("C35")], [— (a formula)], [a list of 514 classes],
  [$M_5$], [formula of App. C #calc("C35")], [— (a formula)], [an explicit count of the 16 types of rank 2 (Table C8)],
)
]

The historical checks involving higher cyclotomic masses concern the auxiliary assignment calculation. The classification-based total uses only the inputs of §5.2 of the paper. Agreement of the two arithmetic organizations cannot prove a structured-pair list complete if a class is absent from both.

== 5.5 Checks on known cases
The tools of §4 and of Appendices B–D were run on cases whose answers are known, and Table 7b lists these checks. They are not part of the proofs; they show that the programs reproduce published values and classical classifications wherever such values exist. The checks of one value by a second route inside this paper are in Table 7a and are not repeated here. The tests and scripts are #calc("C37").

#[
#set text(size: 0.86em)
#set par(justify: false)
#show figure: set block(breakable: true)
#ptable(
  columns: (9em, 1.4fr, 0.9fr, 0.9fr), align: (left, left, left, left),
  caption: [Table 7b. Checks of the tools on known cases #calc("C37"). All agreements are exact unless a standard error is given.],
  table.header([tool], [compared with], [range], [result]),
  [mass formula of Minkowski–Siegel], [[SPLAG, Table 16.4]; the classes of ranks 8, 16, 24 with their $|"Aut"|$ ([SPLAG, Tables 16.1, 16.6, 16.7])], [ranks 8, 16, 24, 32], [agree; the 24 Niemeier contributions agree entry by entry],
  [King's table [Kin03]], [the Minkowski–Siegel mass of rank 32 and Siegel's average number of roots $16320 slash 3617$], [all $13218$ rows], [agree; this checks the table for internal consistency only],
  [rootless mass $m$ #calc("C14")], [King's value [Kin03, Table 1]; the known values in ranks 8, 16, 24], [all root systems of rank at most 32], [agree; the method is King's, so this checks our implementation],
  [rootless masses of 2-elementary genera (§B.2)], [the Niemeier lattices and the Leech lattice [SPLAG, Table 16.1]; for $(24, 2, delta)$ a direct count of sublattices of index 2 of the Niemeier lattices], [rank 24], [agree],
  [unimodular lattices without vectors of norm 1 (start of Lemma B.1)], [numbers [SPLAG, Table 2.2], numbers of roots [SPLAG, Table 16.7], masses [SPLAG, Table 16.2], $|"Aut"| slash |W(R)|$ [SPLAG, Table 17.1]], [ranks at most 24], [agree],
  [existence of 2-elementary genera], [[Nik79, Thm 3.6.2]], [ranks at most 32, all $s$, $delta$], [no disagreement],
  [random gluings of involution components], [exhaustive enumeration in ranks at most 5; the average number of roots $4080 slash 691$], [samples of 2000 per case], [agree within 1.91 standard errors],
  [lists of involution gluings], [estimates of $M$ from random gluings], [types $(16, 16)$, $(8, 24)$], [within 0.44 standard errors for $(16, 16)$; within the error for $(8, 24)$],
  [type $(8, 24)$ (Theorem B.3)], [the counts of the 7-design equations $N_0, dots, N_3 = 12240, 122880, 11520, 0$], [all $14482$ lattices $K$], [agree; each $K$ has average 80 roots over its 270 Niemeier lattices],
  [images of the 120 root lines of $E_8$ in $E_8 slash 3 E_8$ (type $3 - (12, 8) - 8$)], [their number $|"GO"^+_8 (3)| slash |W(E_8)| = 56862$], [all images under $"GO"^+_8 (3)$], [agree],
  [rootless masses over $ZZ[zeta_p]$ (Appendix D)], [the Leech lattice and the character table of $2."Co"_1$], [$p = 5, 7, 13$], [agree (§D.3)],
  [elements of order 3 without fixed vectors], [[Hen13, Thm 1, Table 4]], [the 80 classes of the genus], [the table, as transcribed, sums to the mass of [Hen13, Thm 1]; our eleven lattices give its eleven orders $|"U"|$],
  [elements of order 17], [[Jür15, Satz 3.6.9]], [the four lattices], [agree (§4.2 of the paper)],
  [rootless Hermitian mass over $ZZ[i]$ (Appendix C)], [the classes of ranks 4 and 8; a complete enumeration of 41 small types], [ranks 4, 8], [agree],
  [bliss and `qfauto` in rank 32], [$|"Aut" "BW"_32| = 48126558103142400$; the orders predicted by [KN93] for its twelve lattices; isometries between lattices from different chains], [16 lattices], [agree; generators and isometries checked in integers],
  [the codes of [KN93]], [[KN93, Thm 1] and the twelve codes of [KN93, Table 1]; $sum 1 slash |"Aut" C_0|$ against the value forced by King's mass], [12 codes], [one-to-one; for the column $|"Aut" C|$ see the notes to Table A1],
)
]

*Not checked against known values.* King's table was checked for internal consistency, not row by row against an independent computation of the other $13217$ root systems. The lists of the genera $(20, 6)$ over $ZZ$ and $(13, 5)$ over $ZZ[omega]$ for elements of order 3 are complete by local formulas (§D.4) with no second route. The mass $M_3$ at rank 8 has no value to compare with (§5.4). Together with the last paragraph of §5.4, these are the places where the values rest on our programs alone.

== 5.6 The method in rank 24
In rank 24 the answer is known. The Leech lattice $Lambda$ is the only even unimodular lattice of rank 24 without vectors of norm 2 [SPLAG, Table 16.1], so the analogue of $N$ is 1, and $"Aut" Lambda = "Co"_0$. Every term of the analogue of (3.3) is therefore a sum over conjugacy classes of $"Co"_0$, and can be read from the character table of $2."Co"_1$ in GAP. Table 7c compares these values with the values given by the tools of §4 and of Appendices C and D, run in rank 24 #calc("C38").

This is a check of the tools, not an independent proof. We did not redo §3 in rank 24: the types of elements were taken from the classes of $"Co"_0$, and the comparison uses the uniqueness of $Lambda$. For the involutions and for nine types of odd prime order, the contributions $E$ including the numbers $e(h)$ are compared. The programs that count $e(h)$ #calc("C34", "C8", "C9") were used unchanged; the scripts for rank 24 replace only the places of the outputs and the test of extremality of each gluing, which is not automatic in rank 24. Not included are the contribution $E$ of the type $3 - (12, 0) - 0$ (where $C(h) tilde.equiv 6."Suz"$; only its rootless mass is compared), the types $11 - (2, 4) - 2$ and $23 - (1, 2) - 1$, and a version in rank 24 of the reduction of types in §3. Measured by the values from $"Co"_0$, the terms compared in Table 7c add up to $0.693$ of the total 1.

#[
#set text(size: 0.86em)
#set par(justify: false)
#ptable(
  columns: (1.3fr, 1.4fr, auto), align: (left, left, right),
  caption: [Table 7c. The tools in rank 24 against the classes of $"Co"_0 = "Aut" Lambda$ #calc("C38"). In every row the value computed by the tool equals, as a fraction, the value read from $"Co"_0$: $1 slash |"Co"_0|$, the sum of $1 slash |C(g)|$ over the classes of the given kind, the order of the centralizer of an involution, or the contribution $E = sum e(h) slash |C(h)|$ of a type.],
  table.header([term in rank 24], [computed by], [value]),
  [rootless mass $m$ of rank 24], [King's inclusion–exclusion, as for $m$ in rank 32 #calc("C14")], [$1 slash 8315553613086720000 = 1 slash |"Co"_0|$],
  [$M_2$ ($ZZ[i]$, rank 12)], [the formula of Appendix C], [$1 slash 2012774400$],
  [$M_3$ ($ZZ[zeta_8]$, rank 6)], [the formula of Appendix C.8], [$1 slash 48384$],
  [$M_4$ ($ZZ[zeta_16]$, rank 3)], [Lemma 2.4: the rank is odd, so there is no such lattice], [$0$],
  [order 3, no fixed vectors ($ZZ[omega]$, rank 12)], [rootless mass of Appendix D (§D.4)], [$1 slash 2690072985600$],
  [order 5, no fixed vectors ($ZZ[zeta_5]$, rank 6)], [rootless mass of Appendix D; $e(h)$ by the programs of #calc("C9"), unchanged], [$1 slash 6048000$; $E = 191 slash 3375$],
  [order 7, no fixed vectors ($ZZ[zeta_7]$, rank 4)], [as above], [$1 slash 35280$; $E = 92 slash 2205$],
  [order 13, no fixed vectors ($ZZ[zeta_13]$, rank 2)], [as above], [$1 slash 312$; $E = 1 slash 39$],
  [involutions of type $(16, 8)$], [gluings as in §4.1 of the paper: one double coset, extremal, giving $Lambda$; $e(sigma)$ by the program of #calc("C34"), unchanged], [$|"Stab"(phi)| = 178362777600$; $E = 465333121 slash 178362777600$],
  [involutions of type $(8, 16)$], [the same double coset seen from $-sigma$; $e(-sigma)$ as above], [$E = 5721075841 slash 178362777600$],
  [involutions of type $(12, 12)$], [gluings as in §4.1 of the paper: six double cosets, one extremal, giving $Lambda$ with $sigma tilde -sigma$; $e$ as above], [$|"Stab"(phi)| = 389283840$; $E = 5069 slash 77856768$],
  [type $5 - (4, 8) - 4$], [gluings and $e(h)$ as in §4.2 of the paper, by the programs of #calc("C9"), unchanged; one pair, giving $Lambda$], [$E = 38 slash 1125$],
  [type $5 - (5, 4) - 3$], [as above], [$E = 97 slash 1875$],
  [type $7 - (3, 6) - 3$], [as above], [$E = 8 slash 147$],
  [type $3 - (6, 12) - 6$], [as above, with the programs for order 3 #calc("C8")], [$E = 641 slash 6561$],
  [type $3 - (8, 8) - 8$], [as above], [$E = 26 slash 243$],
  [type $3 - (9, 6) - 5$], [as above], [$E = 3736 slash 19683$],
)
]

For $p = 13$ the walk by neighbours in the genus $(2, 2)$ did not close; the list of rootless lattices in the genus $(2, 0)$ was obtained from it by the step of §B.4 from $(z, s + 2)$ to $(z, s)$, and it is complete because its mass equals the value $1 slash 312$ of the mass formula; identity (B.2), which would need the complete list in $(2, 2)$, is not used.

The row of $M_4$ is 0 by Lemma 2.4, not by formula (C.2). Our implementation of (C.2), applied outside its range to the odd rank 3, returns a negative number: it does not impose the condition that the rank $t$ of a complement and the number $d(C)$ of generators of the 2-part of its discriminant group have the same parity. Appendix C.8 uses (C.2) only for even ranks, and for the three cases used in this paper, $(a, r) = (3, 8), (4, 4), (5, 2)$, no term of the computation violates that condition #calc("C38").



== C.8 The levels $a = 3, 4, 5$

Let $3 <= a <= 5$, $zeta = zeta_(2^a)$, $m = 2^(a-1)$, $cal(O) = ZZ[zeta]$, $KK^+ = QQ(zeta + zeta^(-1))$ with $d_(KK^+) = [KK^+ : QQ] = 2^(a-2)$, $pi.alt_a = 1 - zeta$ (as in §C.3) and $pi.alt_(KK^+) = pi.alt_a overline(pi.alt_a)$, a uniformizer of $KK^+$ at 2. For $g in "Aut" L$ with $g^m = -1$ put $H(x, y) = sum_(k = 0)^(m - 1) (x, y g^k) zeta^k$ (for $a = 2$ this is the form $h$ of §C.1). It is $cal(O)$-valued and Hermitian with $(x, y) = "Tr"_(KK slash QQ) H(x, y) slash m$; since the different of $KK slash QQ$ is $m cal(O)$, $L$ is unimodular if and only if $H$ is, and $L$ is even if and only if $H(x, x) in "Tr"_(KK slash KK^+)(cal(O)) = pi.alt_(KK^+) cal(O)_(KK^+)$ for all $x$ (by $(x, x g^(m-k)) = -(x, x g^k)$, $H(x, x) = (x, x) + sum_(0 < k < m slash 2)(x, x g^k)(zeta^k + zeta^(-k))$, and $"Tr"_(KK slash KK^+)(cal(O)) = 2 cal(O)_(KK^+) + (zeta + zeta^(-1)) cal(O)_(KK^+) = pi.alt_(KK^+) cal(O)_(KK^+)$). These lattices, of rank $r = 32 slash m = 8, 4, 2$, play the role of theta lattices; a root is a vector with $(x, x) = 2$. Everything in §C.1–C.7 goes through with the following changes.

_The genus._ For odd $r$ there are no such lattices (Lemma 2.4). For even $r$ these lattices form one genus, whose mass is $mu_r = 2^(1 - r d_(KK^+)) product_(j = 1)^r |L_(KK^+) (1 - j, chi^j)|$, $chi$ the character of $KK slash KK^+$: this genus is the genus $G_r^KK$ of [Jür15, §3.3.1] (locally hyperbolic at the prime above 2 and $chevron.l 1, dots, 1 chevron.r$ elsewhere), and the mass is [Jür15, Satz 3.3.15].

_Root lattices (Proposition C.1)._ Let $C$ be an irreducible component of the root system and $2^(a - b)$ the length of its orbit under $g$; the stabilizer acts on $C$ through an element $g'$ with $g'^(2^(b-1)) = -1$ (for $b >= 1$). For $b = 1$ the orbit gives $C times.o cal(O)$. For $b >= 2$, as in the proof of Proposition C.1, $C = D_n$ with $g'$ a product of negative $2^(b-1)$-cycles, giving the lattice $G_k^((b))$ of §C.3 induced to $cal(O)$, or $C = D_4$, $b = 2$, or $C = E_8$ with $b = 2$ ($E_8^g$) or $b = 3$: the elements of $W(E_8)$ with characteristic polynomial $Phi_8^2$ form one conjugacy class with centralizer of order 192 #calc("C11"). $b >= 4$ is impossible for $E_8$: for $b = 4$ it would make $E_8$ an even unimodular lattice of rank 1 over $ZZ[zeta_16]$, contradicting Lemma 2.4, and for $b >= 5$ the degree $2^(b-1)$ exceeds 8. Proposition C.3 and the base-change paragraph after it were stated for $kappa <= 2$, but hold for every $kappa$ with the same proofs once $G_k^((kappa))$ carries the normalized trace form $(x, y) = 2^(1 - kappa) "Tr" sum_p x_p overline(y_p)$, which is the restriction of $"Tr"_(KK slash QQ) H slash m$: a root of unity $eta != plus.minus 1$ of 2-power order has trace 0, so the half-edges $(1 - eta) e_p$ and the edges have norm 2 while $2 e_p$ has norm 4; every root of unity is $equiv 1$ modulo $pi.alt_kappa$; $v_(pi.alt_kappa)(1 - eta) = 2^(kappa - j)$ for $eta$ of order $2^j$; and the copies $zeta^j C$ in $C times.o_(cal(O)_b) cal(O)_kappa$ are orthogonal because $"Tr"_(KK_kappa slash KK_(max(b, 1))) zeta^j = 0$ for $0 < j < 2^(kappa - max(b, 1))$. So the coefficients $c(T)$ are obtained as before; for $G_1^((3))$ the coefficient is 0. For $C = E_8$ with $b = 3$, the root sublattices over $cal(O)_3$ were enumerated in the model of $E_8$ on which $g$ acts by two negative 4-cycles of the coordinates: besides 0 and $E_8$ they are of the types $A_1$ (18 of them), $A_2$ (8), $G_2$ (4), $G_1^((3))$ (6), $G_2^((3))$ (3), $A_1^2$ (9), $A_1 perp G_1^((3))$ (6) and $G_1^((3)) perp G_1^((3))$ (3), where $A_1$, $A_2$ and $G_2$ stand for their base changes to $cal(O)_3$ (so $G_1^((2)) times.o cal(O)_3 = A_1 times.o cal(O)_3$ is among the 18); with $c(A_1) = -1$, $c(A_2) = c(G_2) = 2$, $c(A_1^2) = 1$ and $c(G_1^((3))) = c(G_2^((3))) = 0$ (the last by Proposition C.3 for $kappa = 3$, $k = 2$) this gives $c = -16$ #calc("C35").

_The Siegel-type formula (Theorem C.4, Lemma C.5)._ The proof is the same with $QQ$ replaced by $KK^+$: the measure at infinity is the product over the $d_(KK^+)$ real places, and Weil's comparison of Tamagawa numbers holds over any number field [Wei65, n°53–54]. Proposition C.7 holds with the same proof at the prime above 2, whose residue field is $FF_2$, and at an odd prime $frak(p)$ of $KK^+$ with $p^(-j)$ replaced by $N(frak(p))^(-j)$ and $chi(p)$ by $chi(frak(p))$. In Theorem C.8,
$ S'_t = 2^(1 - t d_(KK^+) - t [t "odd"]) product_(j = 1)^t |L_(KK^+) (1 - j, chi^j)|, $
and $f_frak(p)$ is defined by the same sums. The constant is obtained as in the proof of (C.6). At infinity the measure is the product over the real places $v$ of the densities of Lemma C.5, divided by the covolumes of the integral lattices: $cal(O)$ has covolume $2^(-d_(KK^+)) sqrt(|D_KK|) = 2^(1 - d_(KK^+)) D_(KK^+)$ in $KK times.o RR = CC^(d_(KK^+))$ (the relative discriminant of $KK slash KK^+$ has norm 4, so $|D_KK| = 4 D_(KK^+)^2$), and $Lambda_n$, with diagonal in $pi.alt_(KK^+) cal(O)_(KK^+)$, has covolume $(2 sqrt(D_(KK^+)))^n (2^(1 - d_(KK^+)) D_(KK^+))^(n (n - 1) slash 2)$. Hence (C.3) holds with
$ alpha_infinity (T) = c^(d_(KK^+)) N_(KK^+ slash QQ)(det T)^(r - n) (2 sqrt(D_(KK^+)))^n (2^(1 - d_(KK^+)) D_(KK^+))^(n (n - 1) slash 2 - r n). $
The functional equation over $KK^+$ gives $|L_(KK^+) (1 - j, chi^j)| = (2 (j - 1)! (2 pi)^(-j))^(d_(KK^+)) D_(KK^+)^(j - 1 slash 2) 4^((j - 1 slash 2) [j "odd"]) L_(KK^+) (j, chi^j)$. Inserting this and $mu_r = S'_r$ into $mu_r alpha_infinity (T) product_(j = t + 1)^r L_(KK^+) (j, chi^j)^(-1)$, the factors $pi$ and $(j - 1)!$ cancel as before, the powers of $D_(KK^+)$ cancel exactly ($n (2 r - n) slash 2 + n^2 slash 2 - r n = 0$), and the terms of the exponent of 2 proportional to $d_(KK^+)$ add up to $d_(KK^+) (n + t - r) = 0$ after comparison with $2^(-t d_(KK^+))$ in $S'_t$. The remaining exponent is the number $E$ of the proof of (C.6), which is 0. For $t = r$ the formula gives $mu_r$, in agreement with [Jür15, Satz 3.3.15].

_2-adic Jordan decomposition (Lemma C.9)._ Over the completion $cal(O)_(pi.alt_a)$ of $cal(O)$ at $pi.alt_a$, all statements and proofs of Lemma C.9 hold with $2$ replaced by $pi.alt_(KK^+)$, $ZZ_2$ by $cal(O)_(KK^+, pi.alt_(KK^+))$, $1 + 4 ZZ_2$ by $1 + pi.alt_(KK^+)^2 cal(O)_(KK^+)$ and the unit $3$ by $1 + pi.alt_(KK^+)$. The facts used are $cal(O)_(KK^+) slash pi.alt_(KK^+) = FF_2$, $2 in pi.alt_(KK^+) cal(O)_(KK^+)$, $overline(pi.alt_a) equiv pi.alt_a med (mod pi.alt_(KK^+))$, $"Tr"_(KK slash KK^+)(pi.alt_a^k cal(O)) = pi.alt_(KK^+)^(floor((k + 2) slash 2)) cal(O)_(KK^+)$ (the relative different is $pi.alt_a^2 cal(O)$), and $N(cal(O)_(pi.alt_a)^times) = 1 + pi.alt_(KK^+)^2 cal(O)_(KK^+)$: the inclusion holds because $N(w + s pi.alt_a) = w^2 + pi.alt_(KK^+)(w s + s^2) equiv 1 med (mod pi.alt_(KK^+)^2)$ for a unit $w$ and $s in cal(O)_(KK^+)$, and equality because the norm group has index 2 in the units of the ramified extension. In the proof of (iv), the norm of $e_1 + lambda e_2$, $lambda = w + s pi.alt_a$, divided by $pi.alt_(KK^+)^(e + 1)$ is $a + b w^2 + w$ modulo $pi.alt_(KK^+)$, again with derivative 1 in $w$. The dynamic programme of §C.6 then runs unchanged.

#ptable(
  columns: (auto, auto, auto, auto, auto, 1fr), align: (left, right, right, right, right, left),
  caption: [Table C8. Rootless masses for $a = 3, 4, 5$ #calc("C35").],
  table.header([$a$, $r$], [root types], [non-zero], [$M_a = mu_r^emptyset$], [$mu_r$], [independent check]),
  [3, 8], [1270], [174], [$1373459063 slash 7741440$], [$34773704473453 slash 112742891520$], [none at $r = 8$; at $r = 6$ the rootless mass $1 slash 48384$ is $1 slash |C_("Co"_0)(g)|$ for the unique class with $g^4 = -1$],
  [4, 4], [122], [29], [$6983 slash 640$], [$11979143 slash 786432$], [list of 514 classes (238 rootless), both masses exact],
  [5, 2], [21], [7], [$31 slash 24$], [$1455 slash 1024$], [direct count of the 16 types of rank 2, using $h(QQ(zeta_32)) = 1$ [MM76, Main Theorem]],
)

*Checks.* For $a = 3$ and $r = 2, 4, 6$, every root system with non-zero mass gives $1 slash m(R)$ equal to the order of the unitary group of the corresponding lattice ($192$ for $r = 2$; $192^2 dot 2$, $2 dot 46080$ and $8^4 dot 4! slash 2$ for $r = 4$; at $r = 6$ the 17 root systems with non-zero mass map to Niemeier root systems or the Leech lattice), and the rootless mass $1 slash 48384$ at $r = 6$ is the value from the character table of $2."Co"_1$. For $a = 4$, the $pi.alt_a$-neighbours of the 514 classes close up, $sum 1 slash |"U"| = mu_4$ and the rootless part is $6983 slash 640$ exactly, and the 29 root types with non-zero mass match the 29 non-zero values of (C.2) one by one. For $a = 5$, for each root lattice $R$ of rank 2 the lattices with root system exactly $R$ correspond to the $g$-stable Lagrangian submodules of $D(R)$ creating no new root, so $m(R) = \#{C} slash |"U"(R)|$; an explicit count over all submodules gives $sum_R m(R) = 397 slash 3072$ over the 16 types of rank 2, and $mu_2 - 397 slash 3072 = 31 slash 24$, the value of (C.2), the formula giving mass 0 to the four types of rank 1 #calc("C35").



== Historical Table 4 with assignment weights

#ptable(
  columns: 7, align: (left, left, left, right, right, right, right),
  caption: [Table 4. Involution types. $F$ and $K$ are the components (number of classes, $s$, $delta$); "classes" and "extremal" count double cosets $[phi]$ and those giving an extremal $L$, that is, isometry classes of pairs $(L, sigma)$. $M = sum_("extremal") 1 slash |"Stab"(phi)|$ counts both orientations; $E$ is the contribution to (3.3); $E = M$ exactly when $e = 1$ throughout. The pairs of type $(16, 16)$ are in Table 6 of the paper #calc("C4", "C5", "C6", "C34", "C7").],
  table.header([type], [$F$], [$K$], [classes], [extremal], [$M$ (approx.)], [$E$ (approx.)]),
  [$(15, 17)$], [$sqrt(2) A_15^+$; 1, $s = 15$, $delta = 1$], [1 class, $s = 15$, $delta = 1$], [$82086$], [$20502$], [$9664.170$], [$9664.170$],
  [$(14, 18)$], [$sqrt(2) (E_7^2)^+$; 1, $s = 14$, $delta = 1$], [5 classes, $s = 14$, $delta = 1$], [$54649$], [$16590$], [$6876.012$], [$6948.066$],
  [$(12, 20)$], [$sqrt(2) D_12^+$; 1, $s = 12$, $delta = 1$], [66 classes, $s = 12$, $delta = 1$], [$27434$], [$13798$], [$5507.962$], [$5815.972$],
  [$(8, 24)$], [$sqrt(2) E_8$; 1, $s = 8$, $delta = 0$], [$14482$ lattices, $s = 8$], [$14482$], [$14482$], [$6017.813$], [$6185.440$],
  [$(16, 16)$], [10 classes, $8 <= s <= 16$], [same], [$233544$], [$46805$], [$9883.635$], [$10077.337$],
)



= Computation references
#calc-table()
