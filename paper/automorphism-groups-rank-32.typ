#import "style.typ": *
#import "computations-p2.typ": *
#let _version = json("release-files.json").at("version", default: "v(research draft)")
#show: paper.with(
  lang: "en",
  title: [Automorphism groups of extremal even unimodular lattices of rank 32],
  subtitle: [The 99 lattices with the largest groups, and why their groups and their minima are what they are],
  authors: [roiban1344],
  date: [Version #_version.slice(1)],
)


#let _a1 = json("p2-lattices-99.json").lattices
#let _grouped(n) = { let s = str(n); let k = s.len(); s.clusters().enumerate().map(((i, c)) => if i > 0 and calc.rem(k - i, 3) == 0 { "," + c } else { c }).join() }
#let _dlab = _a1.map(r => (r.id, r.label)).to-dict()
#let dl(id) = {
  assert(id in _dlab, message: "unknown lattice id " + id)
  raw(_dlab.at(id))
}

#policy-box[
*Version #_version.slice(1).* A versioned preprint, published with [P1] in the same repository and version; later versions may correct or extend both.
]

#abstract-box[
Among the extremal even unimodular lattices of rank 32, exactly $"143,329"$ have an automorphism group larger than ${plus.minus 1}$ [P1, Thm B]. We describe these groups. For the 99 lattices whose groups have order larger than 2,000 or equal to 512, 1,024 or 1,536 (the groups of the other lattices are identified by `IdGroup` in the Small Groups Library), we determine the generalized Fitting subgroup, a chief series with the action on each factor, a presentation whose order is proved, the conjugacy classes with their frame shapes and fixed lattices, and the orbits on the minimal vectors; 27 of them are in the literature, and we have not found the other 72 there. The orders of the presentations are proved without coset enumeration in the large groups themselves (the largest enumeration is in a quotient of order $"1,036,800"$), by presentations along a chief series and, for the simple sections $Omega_8^+(2)$ and $Omega_10^+(2)$, by a chain of sub-presentations. We then explain four phenomena: the four gluings of two copies of one lattice of rank 16, which correspond to the orbits of lengths $1, 56, 63, 120$ of $"Sp"(6,2)$ on the 240 cosets of the image of $"Aut" F_16$ in the orthogonal group of its discriminant form; the two lattices with an element of order 3 of type $3 - (9, 14) - 7$, which are 3-neighbours, and in the larger of which the type of each element of the centre of $O_3$ is determined by a quadratic form over $FF_3$; the six lattices with a fixed-point-free automorphism of order 5, which are modules over $ZZ[(1 + sqrt(5)) slash 2]$; and the origin of extremality: for 92 of the 99 lattices the absence of roots follows from the two components of a decomposition by an involution or an element of order 3 alone, while for the other 7 no such decomposition has this property (Theorem 8.2).
]

= 1. Introduction <s1>

#policy-box[
*What this section does.* §1.1 recalls the classification of [P1] that this paper starts from. §1.2 states the results in six parts, which are the six later sections. §1.3 explains how lattices are referred to. §1.4 describes the organisation.
]

== 1.1 The lattices with automorphisms other than $plus.minus 1$

Let $cal(L)$ be the set of isometry classes of extremal even unimodular lattices of rank 32, that is, of even unimodular lattices of rank 32 with minimum 4, and let $S = \{L in cal(L) : |"Aut" L| > 2\}$. By [P1, Thm B], $|S| = "143,329"$, and $S = S_"gen" union.sq S_(C_2^2) union.sq S_"CQ"$ is the disjoint union of three subsets.
- $S_(C_2^2)$: the $"67,284"$ lattices with $"Aut" L = \{plus.minus 1, plus.minus sigma\} tilde.equiv C_2 times C_2$.
- $S_"CQ"$: the $"63,499"$ lattices for which $"Aut" L$ is a 2-group with $-1$ as its only involution, that is, cyclic of order 4, 8 or 16 or generalized quaternion of order 8, 16, 32 or 64 [P1, Table 8].
- $S_"gen"$: the other $"12,546"$ lattices. These are the lattices that have an automorphism of odd prime order, or an involution $sigma != plus.minus 1$ with $|C(sigma)| > 4$ or with $sigma$ conjugate to $-sigma$.

The groups are known in $S_(C_2^2)$ and $S_"CQ"$. In $S_"gen"$ they are large and varied, and $S_"gen"$ is the subject of this paper. The 27 previously known lattices (Table 5, [P1, App. A]) all lie in $S_"gen"$ and have $|"Aut" L| >= "4,896"$. [BCDYY25, §3.3] mentions 19 extremal lattices of rank 32 in the Magma database and gives the orders of the automorphism groups of 7 of them.

== 1.2 Results

We prove six results, one in each of §§3–8.

+ *The distribution of groups over $S$* (§3). For the $"143,230"$ lattices of $S$ outside the 99 below, $"Aut" L$ is identified in the Small Groups Library. Tables 2–4 give their numbers by order, by isomorphism type and by coarse structure. Only 15 of these lattices have a non-solvable group (13 groups up to isomorphism), and their only non-abelian simple sections are $A_5$ and $A_6$.
+ *The 99 lattices with large groups* (§4). For each of them a profile (§2.2): the identity, a construction explaining the group, $F^*(G)$ and a chief series, a presentation of proved order, the conjugacy classes and their fixed lattices, and the orbits on minimal vectors. All 99 presentations are proved (Theorem 4.1).
+ *Four gluings of one lattice of rank 16* (§5). The four extremal lattices obtained by gluing two copies of the rank-16 lattice $F_16$ with $"Aut" F_16 tilde.equiv 2^(1+8)_+ . "Sp"(6,2)$ correspond to the orbits of lengths $ell = 1, 56, 63, 120$ of $"Sp"(6,2)$ on the 240 cosets of $overline(G) = "Aut" F_16 slash {plus.minus 1}$ in the orthogonal group $O(q)$ of the discriminant form, and $|C(sigma)| = 2 |"Aut" F_16| slash ell$ (Theorem 5.1).
+ *Elements of type $3 - (9, 14) - 7$* (§6). Exactly two extremal lattices have such an element, #dl("9-14-7:0-0/0") and #dl("12-8-8:0-50/2")\; they are 3-neighbours (Theorem 6.1). In the second, the centre of $O_3(G)$ is a 3-dimensional quadratic space over $FF_3$ whose norm determines the type of each element (Proposition 6.2).
+ *The golden ratio* (§7). Exactly six of the 99 lattices have a fixed-point-free automorphism of order 5, and they are lattices of rank 16 over $ZZ[phi]$, $phi = (1 + sqrt(5)) slash 2$ (Lemma 7.1, Table 9).
+ *Origins of extremality* (§8). For 92 of the 99 lattices there is a decomposition $L supset F plus.o K$ in which the minimal norms of the non-zero discriminant classes satisfy $m_F + m_K > 2$, so that every gluing of $F$ and $K$ is extremal; for the other 7 no decomposition by an involution or an element of order 3 has this property (Theorem 8.2).

== 1.3 Labels

#policy-box[
*Convention.* New lattices are not given names. Every lattice of $S$ is referred to by its label `32.`$|"Aut" L|$`.`$i$ in the database [D] (fixed with version 1 and kept since), where $i$ numbers the lattices with the same group order, for example #dl("8-24:83"). Within one order the lattices are numbered by the isomorphism type of $"Aut" L$, then by the subset $S_"gen"$, $S_(C_2^2)$ or $S_"CQ"$ of [P1, Thm B] in this order, then by a key computed from $L$; this order is a convention of [D] and has no mathematical meaning [P1, §6.8].
]

Table 5 gives the labels of the 27 previously known lattices. For instance, the Barnes–Wall lattice $"BW"_32$ is #dl("known:code:r32"), the last label of $S$, since its group is the largest. The lattices $Lambda_32$ and $Lambda'_32$ of [BN97] are #dl("known:hurwitz:s1-22") and #dl("known:hurwitz:s1-21").

#rem[
*Points that are easy to misread.* These remarks add nothing to the proofs.
+ _Which lattices are "the 99"._ They are chosen by the order of the group alone: $|"Aut" L| > "2,000"$, or $|"Aut" L| in \{512, "1,024", "1,536"\}$, the orders at most 2,000 for which the Small Groups Library has no `IdGroup` (§2.1). No structural property of the lattices enters the choice.
+ _"New"._ A lattice is called new if it is not one of the 27 previously known lattices of [P1, App. A] and we have not found it in the literature we read (§4.1). This is not a claim about databases whose lattices are not published.
+ _The orders $|"Aut" L|$ are those of [P1]._ This paper takes them from [P1] (Table 1a). A presentation "of proved order" (§2.3) is a presentation of $"Aut" L$: its relations hold for integer matrices that generate $"Aut" L$, so the presented group maps onto $"Aut" L$, and Lemmas 2.1–2.4 bound its order by $|"Aut" L|$.
+ _A centralizer order is not $|"Aut" L|$._ Tables 8 and 8a give the order $|C(sigma)|$ or $|C(h)|$ of the centralizer of one element. The group order is this order times the size of the class of the element; for instance the new lattice of §5 has $|"Aut" L| = 2 |C(sigma)|$.
]

== 1.4 Organisation

§2 describes the methods: the profile of each lattice, the four lemmas that prove the orders of presentations, and the computation of conjugacy classes and fixed lattices. §3 gives the distribution of groups over $S$. §4 describes the 99 lattices with large groups: their identity, their presentations, a first classification, their simple sections, and the groups explained by one component. §§5–8 treat four phenomena: the four gluings of one lattice of rank 16, the elements of type $3 - (9, 14) - 7$, the golden ratio, and the origins of extremality. §9 lists the computations with their code, and Appendix A gives one row for each of the 99 lattices. Figures 1 and 2 illustrate Theorems 5.1 and 6.1.

= 2. Methods <s2>

#policy-box[
*What this section does.* §2.1 says which lattices are studied and why. §2.2 defines the profile, the record we compute for each lattice. §2.3 proves the four lemmas that bound the order of a finitely presented group without coset enumeration in the whole group; they are what makes the presentations of §4 proved and not only consistent. §2.4 describes the conjugacy classes, frame shapes and fixed lattices. §2.5 says what is assumed.
]

== 2.1 The 99 lattices

The subset $S_"gen"$ of [P1, Thm B] has $"12,546"$ lattices. For each of them $"Aut" L$ was computed by the full-group algorithm of [P1, §7.1] (bliss for the automorphisms of a local Gram graph, GAP for the stabilizers, and a test of every candidate line) [P1, C1], and its order agrees with the order of the permutation group it induces on the minimal vectors (GAP). If $|"Aut" L| <= "2,000"$ and $|"Aut" L| in.not \{512, "1,024", "1,536"\}$, the group has an `IdGroup` in the Small Groups Library and is identified in §3. The remaining 99 lattices are the subject of §4 #acalc("A1"). Of these, 10 have $|"Aut" L| > 10^7$, 14 have $10^5 < |"Aut" L| <= 10^7$, 50 have $"2,000" < |"Aut" L| <= 10^5$, and 25 have order 512, 1,024 or 1,536 (12, 7 and 6 lattices). The 99 are pairwise non-isometric: they are distinguished by $|"Aut" L|$, by their entries in the lists of [P1, §6.1], and by the number of conjugacy classes of $"Aut" L$ #acalc("A1").

== 2.2 The profile

#thm[
*Definition (profile).* The _profile_ of a lattice $L$ with $G = "Aut" L$ consists of six items.
+ _Identity_: whether $L$ is in the literature, and on what evidence (an isometry checked in integers, or a classification together with invariants).
+ _Construction_: elements $sigma$ with small conjugacy classes, their fixed lattices $F$ and the complements $K$, and the images of $C(sigma)$ in $"Aut" F$ and $"Aut" K$, which explain the size of $G$ when they can.
+ _Structure_: the generalized Fitting subgroup $F^*(G)$, a chief series of $G$ with the action of $G$ on each factor (the order of the image, absolute irreducibility), and whether each factor splits.
+ _Presentation_: generators and relators of a presentation whose order is proved (§2.3).
+ _Action on $L$_: for each conjugacy class its size, its centralizer, its frame shape and its fixed lattice (§2.4), and the orbits of $G$ on the $"146,880"$ minimal vectors.
+ _Relations_ to other lattices: common components, neighbours.
]

#rem[
*Group-theoretic terms.* For a reader coming from lattice theory; these are standard.
- $O_p (G)$ is the largest normal $p$-subgroup of $G$, and the _Fitting subgroup_ $F(G)$ is the product of the $O_p (G)$ over all primes $p$.
- A _component_ of $G$ is a subnormal quasisimple subgroup (perfect, and simple modulo its centre). The _layer_ $E(G)$ is the product of the components.
- The _generalized Fitting subgroup_ is $F^*(G) = F(G) E(G)$. Its centralizer in $G$ lies in it, so much of the structure of $G$ is seen in $F^*(G)$ and in the action of $G$ on it.
- $A compose B$ is a central product, and $Phi(P)$ the Frattini subgroup of a $p$-group $P$.
- A _chief series_ is a chain $1 = N_0 < N_1 < dots < N_k = G$ of normal subgroups of $G$ that cannot be refined by normal subgroups of $G$. Its _chief factors_ $N_(i+1) slash N_i$ are minimal normal subgroups of $G slash N_i$, hence direct products of isomorphic simple groups, and they are determined by $G$ up to order.
]

== 2.3 Presentations of proved order

A presentation $chevron.l X | cal(R) chevron.r$ defines the group $Phi = F(X) slash chevron.l chevron.l cal(R) chevron.r chevron.r$. If a matrix group $G$ is generated by elements satisfying the relators, then $Phi -> G$ is onto. To prove that the presentation is a presentation of $G$ we must show $|Phi| <= |G|$; in general this is done by coset enumeration, which did not fit in 4 GB for some groups of order above $10^5$. The lemmas below give upper bounds for $|Phi|$ from presentations of a fixed shape.

#thm[
*Lemma 2.1 (polycyclic presentations).* Let $Phi$ be generated by $x_1, dots, x_m$ subject to $x_i^(r_i) = w_(i i)$ and $x_i^(-1) x_j x_i = w_(i j)$ ($i < j$), where $r_i >= 1$ and $w_(i i), w_(i j)$ are words in $x_(i+1), dots, x_m$. Then $|Phi| <= r_1 r_2 dots.c r_m$.
]
#proofbox[
Put $Phi_i = chevron.l x_i, dots, x_m chevron.r$ and $Phi_(m+1) = 1$; we show $|Phi_i| <= r_i dots.c r_m$ by descending induction. By induction $Phi_(i+1)$ is finite, and the conjugation relations give $x_i^(-1) Phi_(i+1) x_i subset.eq Phi_(i+1)$; both sides have the same finite order, so equality holds and $x_i$ normalizes $Phi_(i+1)$. Hence $Phi_(i+1) lt.tri.eq Phi_i$, and $Phi_i slash Phi_(i+1)$ is cyclic, generated by the image of $x_i$, of order at most $r_i$ by the power relation. #qed
]

#thm[
*Corollary 2.2.* If elements $s_1, dots, s_m$ of a matrix group $G$ satisfy the relations of Lemma 2.1, generate $G$, and $r_1 dots.c r_m = |G|$, then $x_i |-> s_i$ is an isomorphism $Phi tilde.equiv G$.
]

#thm[
*Lemma 2.3 (presentations along a chief series).* Let the generators of $Phi$ be divided into stages $1, dots, s$, each of one of two kinds.
- A _cyclic stage_ $t$ has one generator $x_t$ and a positive integer $r_t$, with the relation $x_t^(r_t) = $ (a word in the stages after $t$).
- A _block stage_ $t$ has $k_t$ groups of generators $Y_(t, 1), dots, Y_(t, k_t)$. For each group $j$ there is a presentation $chevron.l Y_(t, j) | cal(R)_(t, j) chevron.r$ of a group of known order $|T_t|$ (the presentations may differ between groups). The relations are: each relator of $cal(R)_(t, j)$ equals a word in the stages after $t$, and for generators $a$, $b$ in different groups, $a^(-1) b^(-1) a b$ equals a word in the stages after $t$.
Suppose moreover that for every stage $t$, every generator $x$ of stage $t$ and every generator $y$ of a later stage, $x^(-1) y x$ equals a word in the stages after $t$. Then $|Phi| <= product_t m_t$, where $m_t = r_t$ for a cyclic stage and $m_t = |T_t|^(k_t)$ for a block stage.
]
#proofbox[
Let $Phi_t$ be the subgroup generated by the stages from $t$ on, and $Phi_(s+1) = 1$. We show $|Phi_t| <= m_t dots.c m_s$ by descending induction. By induction $Phi_(t+1)$ is finite. By the last relation, $x^(-1) Phi_(t+1) x subset.eq Phi_(t+1)$ for every generator $x$ of stage $t$, with equality since both sides are finite of the same order. So $Phi_(t+1)$ is normal in $Phi_t$, and $Phi_t slash Phi_(t+1)$ is generated by the images of the generators of stage $t$. For a cyclic stage it is cyclic of order at most $r_t$. For a block stage, the images of group $j$ generate a quotient $H_j$ of $chevron.l Y_(t, j) | cal(R)_(t, j) chevron.r$, so $|H_j| <= |T_t|$. The $H_j$ commute pairwise, and a group generated by pairwise commuting subgroups is their product set $H_1 H_2 dots.c H_(k_t)$, of order at most $|T_t|^(k_t)$. #qed
]

#thm[
*Lemma 2.4 (chain of sub-presentations).* Let $Phi = chevron.l X | cal(R) chevron.r$, let $X_0 subset.eq X$, and let $cal(R)_0 subset.eq cal(R)$ be the relators that are words in $X_0$ alone. Then the subgroup of $Phi$ generated by $X_0$ is a quotient of $chevron.l X_0 | cal(R)_0 chevron.r$, and
$ |Phi| <= [Phi : chevron.l X_0 chevron.r] dot |chevron.l X_0 | cal(R)_0 chevron.r| . $
]
#proofbox[
Sending each letter of $X_0$ to the element of $Phi$ of the same name maps every relator of $cal(R)_0 subset.eq cal(R)$ to $1$, so it defines a homomorphism $chevron.l X_0 | cal(R)_0 chevron.r -> Phi$ with image $chevron.l X_0 chevron.r_Phi$. Hence $|chevron.l X_0 chevron.r_Phi| <= |chevron.l X_0 | cal(R)_0 chevron.r|$, and $|Phi| = [Phi : chevron.l X_0 chevron.r_Phi] dot |chevron.l X_0 chevron.r_Phi|$. #qed
]

Iterating Lemma 2.4 along $X supset X_0 supset X_1 supset dots.c$ bounds $|Phi|$ by the product of the indices, each counted by coset enumeration in the group of the previous stage, and the order of the last stage. If the bound equals $|T|$, the relators hold in $T$ and the images of $X$ generate $T$, then $Phi -> T$ is onto and $|Phi| = |T|$. No Reidemeister–Schreier rewriting is needed. Every finitely presented group in these arguments is rebuilt from the free group and the stored relators, so no order attached by the software to the image of an isomorphism onto a finitely presented group is ever used.

Lemma 2.1 is the case of cyclic stages only. The radical-plus-quotient case of Lemma 2.3 has one block stage with $k = 1$ for $Q = G slash R$, where $R$ is the solvable radical, followed by cyclic stages for $R$; its order $|Q|$ is counted by coset enumeration in $Q$, which is much smaller than $G$. A chief factor $p^d$ gives $d$ cyclic stages with $r = p$. A non-abelian chief factor $T^k$ gives one block stage with $k$ groups, one for each component (minimal normal subgroup of the factor). Its $|T|$ is counted from the stored relators of a presentation of $T$: by coset enumeration, or by Lemma 2.4 for $Omega_8^+(2)$ and $Omega_10^+(2)$.

In each case the generators are taken along a chief series (GAP's `ModuloPcgs` for the abelian factors), and the words on the right of the relations are found by sifting through the stages. All relations are then checked with integer matrices, so an error in the exponents returned by GAP would be detected. The generators of $G$ returned by bliss are written as words in the new generators, and these words are also checked with integer matrices. If the product $product m_t$ equals $|G|$, then $|G| <= |Phi| <= product m_t = |G|$ gives $Phi tilde.equiv G$ as in Corollary 2.2. Here $|G|$ is the order of $"Aut" L$ of §2.1, which agrees with the order of the permutation group on the minimal vectors. In Lemma 2.4 the subset $X_(i+1) subset X_i$ is chosen by dropping up to four generators so that the image in $T$ becomes a proper subgroup with as many relators of $cal(R)$ in the remaining letters as possible. A candidate is rejected if its enumerated index differs from the index of its image in $T$. The choice of $X_(i+1)$ does not affect the correctness of the bound.

#ptable(
  columns: 4, align: (left, left, right, left),
  caption: [Table 1. The chains of sub-presentations for the two largest simple sections #acalc("A4").],
  table.header([lattice], [section $T$], [generators, relators], [chain of indices]),
  [#dl("known:code:walk908")], [$Omega_8^+(2)$, $|T| = "174,182,400"$], [15, 88], [$120 dot 63 dot "23,040"$],
  [#dl("known:code:r32") ($"BW"_32$)], [$Omega_10^+(2)$, $|T| = "23,499,295,948,800"$], [21, 250], [$496 dot 255 dot 126 dot 60 dot "24,576"$],
)

== 2.4 Conjugacy classes and fixed lattices

The conjugacy classes of $G$ are computed with GAP on the permutation representation of $G$ on the minimal vectors, restricted to a union of orbits that spans $L$, where it is faithful. For each class we record its size, the order of its elements and of their centralizer, the frame shape (the characteristic polynomial written as $product_k (x^k - 1)^(a_k)$, abbreviated $product k^(a_k)$), and the fixed lattice: its rank, determinant, discriminant group and minimum. In every case the class sizes add up to $|G|$. The classes of the two largest groups, of #dl("known:code:walk908") and #dl("known:code:r32"), were computed in the same way, on the faithful permutation representations on minimal vectors of degrees $"8,640"$ and $"146,880"$ (GAP's search for a smaller faithful representation returned the same degrees), without using a presentation #acalc("A4"). [P1, §6.3] obtained the classes of these two groups from the structure of the group instead; the classes of involutions and of elements of odd prime order found there, with their centralizer orders, are the same as here.

== 2.5 What is assumed

The orders $|"Aut" L|$ are those of [P1, C1], computed by the full-group algorithm of [P1, §7.1], and each agrees with the order of the permutation group on the minimal vectors computed by GAP. Every automorphism, relation and isometry that we use is checked in integers. What is assumed is listed in Table 1a, in the form of [P1, §7.3, Table 10]. The non-abelian simple sections are named with GAP's `IsomorphismTypeInfoFiniteSimpleGroup`, which identifies a finite simple group from its order and a few invariants. By the classification of finite simple groups these determine the group, so we write these names without further proof. Names of other groups written with $tilde.equiv$ are confirmed by explicit isomorphisms #acalc("A6"), those of $"Aut" F_16 slash O_2$ and of the stabilizers of Theorem 5.1 by the computations of §5 #acalc("A5"). Orthogonal groups are written as in the ATLAS except that the simple ones are $Omega_8^+(2)$ and $Omega_10^+(2)$, and $"GO"$ denotes a full orthogonal group. We use $"U"_4 (2) tilde.equiv S_4(3)$ and $"U"_3 (3) : 2 = G_2(2)$ [ATLAS].

#ptable(
  columns: 2, align: (left, left),
  caption: [Table 1a. What the results of this paper assume, beyond [P1, Table 10]. Outputs checked in integers (automorphisms, relations, words, isometries) are not assumptions.],
  table.header([assumption], [used for]),
  [the orders $|"Aut" L|$ of [P1, C1] (full-group algorithm of [P1, §7.1]) and the orders of permutation groups computed by GAP], [every order; the bounds of the presentations (§2.3, Theorem 4.1)],
  [GAP: coset enumeration (orders of finitely presented quotients and simple sections, the indices of Lemma 2.4)], [Theorem 4.1 #acalc("A3", "A4", "A14")],
  [GAP: conjugacy classes and centralizers], [§2.4; Tables 6–9; the action on $L$ in the profiles #acalc("A2", "A4")],
  [GAP: $F^*(G)$, a chief series, the action on each chief factor, absolute irreducibility and splitting], [Table A1; Tables 4 and 7; the structure in the profiles #acalc("A2")],
  [GAP: `IsomorphismTypeInfoFiniteSimpleGroup` and `StructureDescription` (a description, not a proof of isomorphism); the maps of `IsomorphismGroups` are checked to be bijective homomorphisms], [the names of simple sections and of groups (Tables 3, 7, A1) #acalc("A6")],
  [bliss and GAP: the automorphism groups of the components and the images of $C(sigma)$, in particular that the image contains $O_2("Aut" F)$], [§4.5, §5, Theorem 5.1 #acalc("A5")],
  [PARI `qfminim`: the enumeration of the short vectors of the duals of the components is complete], [Lemma 8.1 applied in Theorem 8.2, Table 10 #acalc("A11")],
  [our code: the commutant algebras and the invariants of the fixed lattices (rank, determinant, discriminant group, minimum)], [§2.4, §4.1, §4.6, §6 #acalc("A2", "A7", "A8")],
)

= 3. The distribution of groups over $S$ <s3>

#policy-box[
*What this section does.* We give the orders and isomorphism types of $"Aut" L$ over all of $S$. The data come from the identification of the groups of $S_"gen"$ [P1, C30], recorded in `data/classified-lattices.jsonl` of [P1]; this section only describes them and points to that file. The groups in $S_(C_2^2)$ and $S_"CQ"$ are determined in [P1, §6.4–6.5]; the variety is in $S_"gen"$.
]

#ptable(
  columns: 8, align: (right, right, right, right, right, right, right, right),
  caption: [Table 2. Number of lattices of $S$ by $|"Aut" L|$ and by the subsets of [P1, Thm B], read down the columns, for the orders at most 2,000 other than $512, "1,024", "1,536"$, for which the Small Groups Library has an `IdGroup`. The last row collects the 99 lattices of §4: $|"Aut" L| > "2,000"$, and the orders $512, "1,024", "1,536"$, for which it has none. Order by order the three subsets add up to [P1, Table 9]. In $S_(C_2^2)$ the group is $C_2 times C_2$ throughout; in $S_"CQ"$ it is $C_4$ at order 4 and $C_(2^n)$ or $Q_(2^n)$ above [P1, Table 8].],
  table.header([$|"Aut" L|$], [$S_"gen"$], [$S_(C_2^2)$], [$S_"CQ"$], [$|"Aut" L|$], [$S_"gen"$], [$S_(C_2^2)$], [$S_"CQ"$]),
  [$4$], [], [$"67,284"$], [$"62,468"$], [$192$], [$23$], [], [],
  [$6$], [$"2,308"$], [], [], [$216$], [$4$], [], [],
  [$8$], [$"7,055"$], [], [$968$], [$224$], [$1$], [], [],
  [$10$], [$6$], [], [], [$240$], [$3$], [], [],
  [$12$], [$674$], [], [], [$256$], [$27$], [], [],
  [$16$], [$"1,233"$], [], [$55$], [$288$], [$19$], [], [],
  [$18$], [$4$], [], [], [$320$], [$2$], [], [],
  [$20$], [$13$], [], [], [$324$], [$1$], [], [],
  [$24$], [$272$], [], [], [$384$], [$21$], [], [],
  [$32$], [$284$], [], [$6$], [$432$], [$2$], [], [],
  [$36$], [$17$], [], [], [$480$], [$5$], [], [],
  [$40$], [$23$], [], [], [$576$], [$3$], [], [],
  [$48$], [$163$], [], [], [$672$], [$1$], [], [],
  [$54$], [$2$], [], [], [$768$], [$10$], [], [],
  [$64$], [$110$], [], [$2$], [$800$], [$1$], [], [],
  [$72$], [$17$], [], [], [$864$], [$1$], [], [],
  [$80$], [$5$], [], [], [$896$], [$1$], [], [],
  [$96$], [$39$], [], [], [$960$], [$2$], [], [],
  [$108$], [$2$], [], [], [$"1,152"$], [$5$], [], [],
  [$120$], [$8$], [], [], [$"1,440"$], [$2$], [], [],
  [$128$], [$61$], [], [], [$"1,728"$], [$2$], [], [],
  [$144$], [$11$], [], [], [$"1,920"$], [$1$], [], [],
  [$160$], [$3$], [], [], [], [], [], [],
  [$> "2,000"$ or $512, "1,024", "1,536"$], [99], [], [], [], [], [], [],
)

#ptable(
  columns: 4, align: (right, right, right, left),
  caption: [Table 3. Isomorphism types of $"Aut" L$ in $S_"gen"$ with an `IdGroup` ($|"Aut" L| <= "2,000"$, not $512, "1,024", "1,536"$): the number of lattices, the number of distinct groups, and the most frequent group with its number of lattices. Names are GAP's `StructureDescription` of the `IdGroup`, a description and not a proof of isomorphism.],
  table.header([$|"Aut" L|$], [lattices], [groups], [most frequent `IdGroup`]),
  [$6$], [$"2,308"$], [$1$], [`[6,2]` `C6`],
  [$8$], [$"7,055"$], [$3$], [`[8,2]` `C4 x C2` ($"3,737"$)],
  [$10$], [$6$], [$1$], [`[10,2]` `C10`],
  [$12$], [$674$], [$4$], [`[12,4]` `D12` ($278$)],
  [$16$], [$"1,233"$], [$12$], [`[16,3]` `(C4 x C2) : C2` ($283$)],
  [$18$], [$4$], [$2$], [2 groups with 2 each],
  [$20$], [$13$], [$3$], [`[20,4]` `D20` ($6$)],
  [$24$], [$272$], [$12$], [`[24,5]` `C4 x S3` ($61$)],
  [$32$], [$284$], [$30$], [`[32,22]` `C2 x ((C4 x C2) : C2)` ($72$)],
  [$36$], [$17$], [$4$], [`[36,13]` `C2 x ((C3 x C3) : C2)` ($8$)],
  [$40$], [$23$], [$5$], [`[40,5]` `C4 x D10` ($8$)],
  [$48$], [$163$], [$29$], [`[48,30]` `A4 : C4` ($25$)],
  [$54$], [$2$], [$2$], [2 groups with 1 each],
  [$64$], [$110$], [$48$], [`[64,60]` `(C2 x ((C4 x C2) : C2)) : C2` ($10$)],
  [$72$], [$17$], [$9$], [`[72,46]` `C2 x S3 x S3` ($5$)],
  [$80$], [$5$], [$4$], [`[80,14]` `(C20 x C2) : C2` ($2$)],
  [$96$], [$39$], [$24$], [`[96,87]` `((C4 x C2) : C2) x S3` ($5$)],
  [$108$], [$2$], [$2$], [2 groups with 1 each],
  [$120$], [$8$], [$6$], [2 groups with 2 each],
  [$128$], [$61$], [$37$], [`[128,2163]` `C2 x C2 x ((C2 x C2 x C2 x C2) : C2)` ($8$)],
  [$144$], [$11$], [$8$], [`[144,128]` `S3 x SL(2,3)` ($3$)],
  [$160$], [$3$], [$2$], [`[160,81]` `C5 : ((C4 x C2) : C4)` ($2$)],
  [$192$], [$23$], [$18$], [`[192,991]` `(C2 x (A4 : C4)) : C2` ($5$)],
  [$216$], [$4$], [$4$], [4 groups with 1 each],
  [$224$], [$1$], [$1$], [`[224,75]` `((C4 x C2) : C2) x D14`],
  [$240$], [$3$], [$3$], [3 groups with 1 each],
  [$256$], [$27$], [$22$], [`[256,26531]` `C2 x ((((C2 x C2 x C2) : (C2 x C2)) : C2) : C2)` ($3$)],
  [$288$], [$19$], [$11$], [`[288,851]` `S3 x GL(2,3)` ($5$)],
  [$320$], [$2$], [$2$], [2 groups with 1 each],
  [$324$], [$1$], [$1$], [`[324,63]` `C2 x (((C9 x C3) : C3) : C2)`],
  [$384$], [$21$], [$21$], [21 groups with 1 each],
  [$432$], [$2$], [$2$], [2 groups with 1 each],
  [$480$], [$5$], [$4$], [`[480,943]` `C4 x S5` ($2$)],
  [$576$], [$3$], [$3$], [3 groups with 1 each],
  [$672$], [$1$], [$1$], [`[672,1122]` `C7 : (C2 x (((C4 x C2) : C2) : C3))`],
  [$768$], [$10$], [$10$], [10 groups with 1 each],
  [$800$], [$1$], [$1$], [`[800,957]` `(C5 x C5) : ((C4 x C2) : C4)`],
  [$864$], [$1$], [$1$], [`[864,4687]` `C2 x (((C3 x C3) : Q8) x S3)`],
  [$896$], [$1$], [$1$], [`[896,12517]` `(((C2 x C2 x C2 x C2) : C2) : C2) x D14`],
  [$960$], [$2$], [$2$], [2 groups with 1 each],
  [$"1,152"$], [$5$], [$4$], [`[1152,157851]` `C2 x S4 x S4` ($2$)],
  [$"1,440"$], [$2$], [$2$], [2 groups with 1 each],
  [$"1,728"$], [$2$], [$2$], [2 groups with 1 each],
  [$"1,920"$], [$1$], [$1$], [`[1920,241002]` `(C2 x C2 x C2 x C2) : SL(2,5)`],
)

#ptable(
  columns: 3, align: (left, right, right),
  caption: [Table 4. Coarse structure of $"Aut" L$ over $S$. The middle column counts the $"143,230"$ lattices of $S$ outside the 99, all of whose groups have an `IdGroup` ($S_"gen"$, $S_(C_2^2)$, $S_"CQ"$); the class is read from the Small Groups Library. The last column is from the profiles of the 99.],
  table.header([class of groups], [$L in S$ outside the 99], [the 99]),
  [2-groups], [$"139,553"$], [29],
  [solvable, not 2-groups], [$"3,662"$], [40],
  [non-solvable], [$15$], [30],
)

Among the lattices outside the 99 (those whose group has an `IdGroup`), a non-solvable group occurs for $15$ lattices, all in $S_"gen"$, with $13$ groups. This set is disjoint from the $30$ lattices with a non-solvable group among the 99 (Table 4, last column; §4). Of order $120$: $"SL"(2,5)$ (two lattices) and $C_2 times A_5$. Of order $240$: $A_5 : C_4$ and $C_4 times A_5$. Of order $480$: $C_4 times S_5$ (two lattices), $C_2 times ("SL"(2,5) : C_2)$, $A_5 times D_8$ and $C_2^2 times S_5$. Of order $960$: $A_5 : (C_4 : C_4)$ and $"SL"(2,5) : (C_4 times C_2)$. Of order $"1,440"$: $S_3 times (A_5 : C_4)$ and $C_2 times (A_6 : C_2)$. Of order $"1,920"$: $C_2^4 : "SL"(2,5)$. (Names are GAP's `StructureDescription`.) Outside the 99, the only non-abelian simple sections are $A_5$ and $A_6$. The tables are made from the identification lists by a script that also checks them against [P1, Table 9] #acalc("A13").

Outside the 99, $-1$ is the only central involution of $"Aut" L$ for $"67,844"$ of the $"143,230"$ lattices. These are all $"63,499"$ lattices of $S_"CQ"$, and $"4,345"$ of the $"12,447"$ lattices of $S_"gen"$ outside the 99. In $S_(C_2^2)$ the group $C_2 times C_2$ has three central involutions. A normal subgroup $O_3 != 1$ occurs for $"3,394"$ lattices outside the 99, all in $S_"gen"$. Non-abelian simple sections occur only in the 15 lattices listed above #acalc("A13").

= 4. The 99 lattices with large groups <s4>

#policy-box[
*What this section does.* §4.1 identifies the 27 previously known lattices among the 99. §4.2 states that all 99 presentations are proved. §4.3 gives a first classification of the 99 by central and index-2 elements, and §4.4 the simple sections. §4.5 explains the groups of the lattices with a central involution by one of the two components. The four phenomena of §§5–8 are taken up separately. The full table is Appendix A; the profiles themselves are in the data file `data/largest-groups.jsonl` (Table 12).
]

== 4.1 Identity

Of the 99 lattices, 27 are in the literature; they are the previously known lattices of [P1, App. A], and Table 5 lists them under the constructions and classifications in which they occur. For 19 of them there is an explicit isometry to a Gram matrix of the catalogue [NS]; that a page of [NS] shows the lattice of a given source, and the grade of each identification, are those of [P1, App. A, Table A1]. The catalogue [NS] has 79 pages of rank 32; 25 of them are extremal even unimodular lattices, and each of these 25 is isometric to one of the 99 lattices. We found the isometry $T$ with bliss and checked $T^top G T = G'$ and $|det T| = 1$ in integers #acalc("A12"). The 25 pages give 19 distinct lattices, since several lattices have more than one page; for instance $"BW"_32$ = #dl("known:code:r32") has five. For the other 8 previously known lattices [NS] has no Gram matrix. Their identity rests on the classifications cited in Table 5 together with the agreement of the invariants that distinguish the lattices within each classification ($|"Aut" L|$, and the order $|"U"|$ of the unitary group in [Hen13]), as in [P1, App. A]; we do not claim an isometry for them.

One of the 19 is #dl("12-20-c63:236"), which [P1, Table A1] identifies only through the page `dim32f5` of [NS] and, via that page, with the lattice $[4.L_3(4).2^2]_(32,2)$ of [Neb98q, §21, Thm 21.1]. Here we also compare two invariants given there. The decomposition $"8,064" + "20,160" + 2 dot "23,040" + "32,256" + "40,320"$ of its minimal vectors into orbits agrees with ours. [Neb98q] gives the commutant of the group as $QQ[sqrt(-5)]$, and we confirmed this. The commutant $"End"_G (QQ^32)$ has dimension 2. It is spanned by $1$ and an integral matrix $Z$ with $Z^2 = -5$, $(Z x, Z y) = 5 (x, y)$ and $(Z x, y) = -(x, Z y)$ #acalc("A15"). So $L$ is a module of rank 16 over $ZZ[sqrt(-5)]$, and $"Aut" L$ acts $ZZ[sqrt(-5)]$-linearly. The element $g$ of order 4 that generates $F(G) tilde.equiv C_4$ is not central, because $G$ contains elements that invert it. The commutant has dimension 4 for $C_G (g)$ and for the derived group $G'$, which have index 2 and 4 in $G$.

The other 72 lattices are not isometric to any of the 27, since the 99 are pairwise distinguished by their invariants (§2.1). We have not found any of the 72 in the literature we read; in this sense we call them _new_ below. [BCDYY25, §3.3] lists $|"Aut" L|$ for 7 of the 19 extremal lattices of rank 32 in the Magma database; each of these 7 orders is the order of exactly one previously known lattice among the 99 and of none of the 72. For the other 12 lattices of that database neither orders nor Gram matrices are given, so we cannot exclude that some of the 72 are among them.

#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 5. The 27 previously known lattices among the 99, by the construction or classification in which they occur; a lattice occurring in several rows is listed in each. "$T$" means an explicit isometry to a Gram matrix of [NS], checked in integers #acalc("A12")\; the link from a page of [NS] to the source of the row, and its grade, are those of [P1, Table A1]\; "class." means the classification cited in the first column together with the agreement of the invariants that distinguish its lattices.],
  table.header([layer], [labels], [isometry to [NS] or classification]),
  [code lattices, neighbour defect 0 [KV89]], [#dl("known:code:r32"), #dl("known:code:walk908"), #dl("known:code:cyc0"), #dl("known:kmax24:9"), #dl("known:kmax24:10")], [$T$ (all 5)],
  [the twelve lattices of [KN93] (neighbour defect 8 for ten; $S_1$ = #dl("known:kmax24:10") and $S_2$ = #dl("known:kmax24:9") have neighbour defect 0)], [#dl("known:kmax24:0"), #dl("known:kmax24:1"), #dl("known:kmax24:2"), #dl("known:kmax24:3"), #dl("known:kmax24:4"), #dl("known:kmax24:5"), #dl("known:kmax24:6"), #dl("known:kmax24:7"), #dl("known:kmax24:8"), #dl("known:kmax24:9"), #dl("known:kmax24:10"), #dl("known:kmax24:11")], [$T$ (all 12)],
  [Eisenstein, the 11 classes without roots, nos. 70–80 [Hen13]], [#dl("known:code:r32"), #dl("known:hurwitz:s1-21"), #dl("known:hurwitz:s1-22"), #dl("known:kmax24:8"), #dl("known:seven:s1-4"), #dl("known:zeta17:s1-14"), #dl("known:eis:s1-880"), #dl("known:eis:s1-100"), #dl("known:eis:s1-440"), #dl("known:eis:s1-140"), #dl("known:eis:s1-60") (nos. 70–80 in this order)], [$T$ for nos. 70–73, 75; class. and $|"U"|$ for nos. 74, 76–80],
  [$QQ(zeta_17)$ [Jür15, Satz 3.6.9]], [#dl("known:zeta17:s1-1"), #dl("known:zeta17:s1-10"), #dl("known:zeta17:s1-14"), #dl("known:code:r32")], [$T$ for the last two; class. for the first two],
  [r.i.m.f. groups [Neb98q, §21]], [#dl("known:code:r32"), #dl("known:hurwitz:s1-21"), #dl("12-20-c63:236"), #dl("known:zeta17:s1-14")], [$T$ (all 4)],
  [Hurwitz [BN97]], [#dl("known:code:r32"), #dl("known:hurwitz:s1-21"), #dl("known:hurwitz:s1-22")], [$T$ (all 3)],
)

== 4.2 Presentations

#thm[
*Theorem 4.1.* Each of the 99 groups has a presentation whose order is proved by Lemmas 2.1–2.4 and equals $|"Aut" L|$: by Lemma 2.1 for the 69 solvable groups, by the radical-plus-quotient case of Lemma 2.3 for 27 groups (for #dl("12-20-c63:236") with the quotient bounded by Lemma 2.4), and by Lemma 2.3 with block stages for #dl("known:hurwitz:s1-21"), #dl("known:code:walk908") and #dl("known:code:r32"), the simple sections of the last two by Lemma 2.4 #acalc("A3", "A4").
]

In the radical-plus-quotient case the order of the presentation of $Q = G slash R$ is counted by coset enumeration in $Q$, which is much smaller than $G$. One quotient did not fit: for #dl("12-20-c63:236") we have $|Q| = "80,640"$, and the enumeration on the trivial subgroup exceeded $3 dot 10^7$ cosets. Here one generator $y$ of the presentation $Q_"fp"$ satisfies the relator $y^4$. Lemma 2.4 with $X_0 = {y}$ and $cal(R)_0 = {y^4}$ gives $|Q_"fp"| <= [Q_"fp" : chevron.l y chevron.r] dot 4$. The index, counted by coset enumeration on $chevron.l y chevron.r$, is $"20,160"$, so $|Q_"fp"| <= "80,640" = |Q|$ #acalc("A3").

The presentations of Table 12 are not the only presentations of these groups that occur in our computations. The group identification run of §3 also produced, for each of the $"12,536"$ lattices of $S_"gen"$ with $|"Aut" L| <= 10^7$, generators and relators that hold in $"Aut" L$. For nine of them, with $|"Aut" L|$ between $"196,608"$ and $"393,216"$, the order of that presentation was not counted in the run. We counted it separately for all nine. Coset enumeration over the stabilizer of a point in a permutation representation of small degree gives the index, Reidemeister–Schreier rewriting gives a presentation of the stabilizer, and the last group of the chain is enumerated on the trivial subgroup. The product equals $|"Aut" L|$ in all nine cases #acalc("A14").

== 4.3 A first classification

For an element $sigma$ with a small conjugacy class the group is close to $C(sigma)$. If $sigma$ is central, then $G = C(sigma)$. If the class has two elements, $\{sigma, -sigma\}$ or $\{g, g^(-1)\}$, then $G = C(sigma) . 2$. In both cases $C(sigma)$ maps to the automorphism groups of the fixed lattice $F = L^sigma$ and of its orthogonal complement, and its image consists of the pairs that preserve the gluing [P1, Prop. 4.1]. Table 6 sorts the 99 lattices accordingly. The types of involutions are those of [P1, §4], and the type $3 - (9, 14) - 9$ of an element of order 3 is in the notation of [P1, §3.3].

#ptable(
  columns: 3, align: (left, right, left),
  caption: [Table 6. First classification of the 99 lattices by central elements and elements with a class of size 2. Each lattice is counted in the first row, in the order of the table, whose condition it satisfies. A lattice with central involutions of several types is therefore counted under $(12, 20)$ before $(16, 16)$, and under either of these before $(8, 24)$; the row $(8, 24)$ contains the lattices whose only central involutions other than $-1$ are of type $(8, 24)$.],
  table.header([family], [number], [examples]),
  [previously known (Table 5)], [27], [#dl("known:code:walk908"), #dl("known:zeta17:s1-1"), #dl("12-20-c63:236")],
  [central involution of type $(12, 20)$], [7], [#dl("12-20-c62:148")],
  [central involution of type $(16, 16)$], [11], [#dl("16-16d8:212783"), #dl("16-16-s16d0a:678")],
  [central involution of type $(8, 24)$], [15], [#dl("8-24:83"), #dl("8-24:8")],
  [central $g$ with $g^2 = -1$], [4], [#dl("16-16-s12d1:19"), #dl("14-18e:602"), #dl("16-16d8:98086"), #dl("8-16-8:3-0/4")],
  [central involution of type $(14, 18)$; central element of order 3 of type $3 - (9, 14) - 9$], [1 + 1], [#dl("14-18e:601"), #dl("9-14-9:0-0/131")],
  [an element with a class of size 2: an involution $sigma tilde.op -sigma$ of type $(16,16)$, an element $g tilde.op g^(-1)$ with $g^2 = -1$, or an element $h$ of order 3 with $h tilde.op h^(-1)$], [19], [#dl("16-16-s10d1:0"), #dl("9-14-7:0-0/0"), #dl("16-16-s12d1:69")],
  [none of these], [14], [#dl("5-8-0-0:8"), #dl("12-8-8:0-50/2"), #dl("5-5-12-5:1")],
)

The type $(8, 24)$ is placed last because its component of rank 8 is always $sqrt(2) E_8$, and the group is then that of the other component (§4.5). Among the other 72 lattices, 26 have a central involution of type $(8, 24)$; 11 of them also have one of type $(12, 20)$ or $(16, 16)$ and are counted in those rows. For the lattices of the last row no single fixed lattice explains the group; for #dl("5-8-0-0:8") the smallest conjugacy class of an element other than $plus.minus 1$ has 160 elements. §4.6 and §6 treat two of these lattices.

== 4.4 Simple sections

#ptable(
  columns: 2, align: (left, right),
  caption: [Table 7. Non-abelian simple groups $T$ in chief series of the 30 non-solvable groups: the number of lattices with a chief factor $T$ or $T^k$. The chief factors $A_5^4$ (#dl("known:hurwitz:s1-21")) and $A_6^2$ (#dl("5-8-0-0:8")) are counted under $A_5$ and $A_6$; each group has exactly one non-abelian chief factor, so the numbers add up to 30.],
  table.header([$T$], [lattices]),
  [$A_5$], [9], [$L_2(7)$], [6], [$A_6$], [5], [$L_2(17)$], [3], [$"U"_4 (2) tilde.equiv S_4(3)$], [2],
  [$L_3(4)$, $"U"_3 (3)$, $L_2(31)$, $Omega_8^+(2)$, $Omega_10^+(2)$], [1 each],
)

Of the 99 groups, 69 are solvable, and 29 of these are 2-groups (Table 4). In every group $C_G (F^*(G)) <= F^*(G)$, which we checked as a test of the computation of $F^*(G)$. The non-abelian chief factors of the 30 non-solvable groups are simple except for two, $A_6^2$ for #dl("5-8-0-0:8") and $A_5^4$ for #dl("known:hurwitz:s1-21"). The simple groups are named by their order and invariants (§2.5). A quasisimple group can be the whole group: for #dl("known:zeta17:s1-1") we have $G = F^*(G) = E(G) tilde.equiv "SL"(2, 17)$, a non-split extension of $L_2(17)$ by its centre $chevron.l -1 chevron.r$.

Table 7a lists the names that we confirmed by an explicit isomorphism. In each case GAP's `IsomorphismGroups` returned a map onto a group given by name, and we checked that the map is a bijective homomorphism #acalc("A6"). Sixteen names were confirmed. For two groups every candidate we tried was rejected. For #dl("16-16d8:68772") the group $G$ is not isomorphic to $C_2 times "Aut"(A_6)$, $C_2^2 times S_6$, $C_2^2 times "PGL"(2, 9)$ or $C_4 times S_6$. For #dl("12-20-c62:148") it is not isomorphic to $C_2^2 times "Aut"(A_6)$, $D_8 times S_6$, $C_2^3 times S_6$ or $C_2 times C_4 times S_6$. In both, $E(G) = A_6$.

#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 7a. Group names confirmed by explicit isomorphisms #acalc("A6"). $E$ is the layer, $R$ the solvable radical; $"SL"(2,9) compose "SL"(2,9)$ is the central product with the central elements $-1$ identified.],
  table.header([lattice], [group], [isomorphic to]),
  [#dl("9-14-7:0-0/0"), #dl("9-14-9:0-0/131")], [$G slash F^*$], [$"GL"(2, 3)$],
  [#dl("5-5-12-5:1")], [$G$], [$C_2 times "AGL"(2, 5)$],
  [#dl("known:zeta17:s1-1")], [$G$], [$"SL"(2, 17)$],
  [#dl("16-16-s10d1:0")], [$G slash F^*$], [$"GO"_6^-(2) tilde.equiv W(E_6)$],
  [#dl("16-16-s12d1:69")], [$E$; $G slash R$], [$"U"_4 (2)$; $"GO"_6^-(2)$],
  [#dl("8-16-8:5-0/19")\; #dl("12-8-8:0-49/5")], [$E$], [$"SL"(2, 7)$; $"SL"(2, 5)$],
  [#dl("5-8-0-0:8")], [$E$], [$"SL"(2, 9) compose "SL"(2, 9)$],
  [#dl("8-24:1205")], [$E$], [$L_2(7)$],
  [#dl("14-18e:602"), #dl("16-16d8:98086")], [$G slash F^*$], [$S_5$],
  [#dl("16-16d8:222101")], [$G slash F^*$], [$S_5 times S_3$],
  [#dl("12-8-8:0-49/3")], [$G slash F^*$], [$D_8 times S_3$],
  [#dl("known:code:cyc0")], [$G slash F^*$], [$L_2(31)$],
)

== 4.5 Groups determined by one component

For 96 of the 99 lattices, those with $|"Aut" L| <= 2 dot 10^(10)$ other than #dl("known:zeta17:s1-1"), whose only involution is $-1$, we took an involution $sigma != plus.minus 1$ in a smallest conjugacy class. For each such $sigma$ we computed the fixed lattice $F$, its orthogonal complement $K$, their automorphism groups, and the images of $C(sigma)$ in $"Aut" F$ and $"Aut" K$ #acalc("A5").

*Central involutions of type $(8, 24)$.* Here the component of rank 8 is $sqrt(2) E_8$, and $C(sigma) = G$. By [P1, §4.1.3, type $(8, 24)$] every $v in "Aut" K$, for the component $K$ of rank 24, extends to exactly two elements of $C(sigma)$. So $G -> "Aut" K$ is onto with kernel $\{1, -sigma\}$, the two elements acting trivially on $K$ (on which $sigma$ acts as $-1$), and $|G| = 2 |"Aut" K|$. The group of $L$ is the group of its component of rank 24. In all cases the computed image of $G$ in $"Aut" K$ has order $|G| slash 2$. For the six lattices where $|"Aut" K|$ was also computed directly, namely #dl("8-24:609"), #dl("8-24:680"), #dl("8-24:14202"), #dl("8-24:665"), #dl("16-16-s14d0:628") and #dl("16-16-s16d0a:261"), the image is all of $"Aut" K$. The image of $G$ in $W(E_8)$ is small: its index is at least $"80,640"$ in every case.

*Central involutions of type $(16, 16)$ with $F tilde.equiv K$.* Eight new lattices have such an involution as the $sigma$ chosen above. In seven of them the component $F$ is one of two lattices of [P1, Table 6], both with $|"Aut" F| = "53,271,016,243,200"$: the component with $(s, delta) = (14, 0)$, discriminant group $2^(14)$, for #dl("12-20-c37:49"), #dl("16-16-s14d0:233"), #dl("16-16-s14d0:301"), #dl("16-16-s14d0:503") and #dl("16-16-s14d0:840"), and the component with $(s, delta) = (16, 1)$, discriminant group $2^(16)$, for #dl("16-16-s16d0a:678") and #dl("16-16d8:212783"). For these seven, the image of $G$ in $"Aut" F$ has order between $128$ and $"4,096"$, so its index is more than $10^(10)$. The eighth, #dl("16-16-s12d1:44"), has as $F$ the component with $(s, delta) = (12, 1)$ of [P1, Table 6], with $|"Aut" F| = "754,974,720"$; there the image of $G$ has order $"1,024"$ and index $"737,280"$.

*A central involution of type $(12, 20)$.* For #dl("12-20-c62:148") the image of $G$ in the automorphism group of the component of rank 20 (determinant $2^(12)$, $|"Aut"| = "2,880"$) is the whole group.

*A central $g$ with $g^2 = -1$.* For the four lattices of Table 6 with such an element, $L$ is a Hermitian lattice of rank 16 over $ZZ[i]$ with $i$ acting as $g$. Since $g$ is central, every automorphism commutes with $g$, so $"Aut" L$ is the unitary group $"U"(L)$ of this Hermitian lattice. The groups are $2^8 . (S_3 times S_3)$ for #dl("16-16-s12d1:19"), $2^6 . S_5$ for #dl("14-18e:602") and #dl("16-16d8:98086"), and $(2^6 times 3) . D_(12)$ for #dl("8-16-8:3-0/4"), where the dots denote extensions with the given orders and quotients. The two groups $2^6 . S_5$ are not isomorphic: $O_2(G)$ is elementary abelian for the first and has centre of order 4 and Frattini subgroup of order 2 for the second.

== 4.6 Two copies of one lattice of rank 16 in the genus $"II"_16 (3^(+8))$

The lattice $L$ = #dl("5-8-0-0:8"), with $|"Aut" L| = 2^(10) 3^4 5^2$, is in the last row of Table 6: no conjugacy class of an element other than $plus.minus 1$ has fewer than 160 elements. Its layer is $E(G) tilde.equiv "SL"(2, 9) compose "SL"(2, 9)$, of order $"259,200"$, with $E slash Z(E) tilde.equiv A_6 times A_6$, since $"SL"(2, 9) slash {plus.minus 1} tilde.equiv A_6$, and $G slash F^*(G) tilde.equiv 2^3$, the group $[8, 5]$ of the Small Groups Library #acalc("A2"). The structure of $L$ is seen from $E$ #acalc("A9").

- Restricted to either factor $"SL"(2, 9)$, the character of $QQ^32 times.o CC$ is four times the sum of the two irreducible characters of degree 4. Restricted to $E$ it is the sum of two irreducible characters of degree 16. An irreducible representation of a central product is a tensor product of irreducible representations of the factors, so $QQ^32 = V plus.o V'$ with $V$ and $V'$ of dimension 16, each the tensor product of a 4-dimensional representation of one factor with a 4-dimensional representation of the other.
- The commutant $"End"_E (QQ^32)$ has dimension 2, and the eigenspaces of a non-scalar element are $V$ and $V'$. The sublattices $U = L inter V$ and $U' = L inter V'$ have rank 16, determinant $3^8$, discriminant group $3^8$, minimum 4 and $|"Aut"| = "1,036,800"$, and they are isometric: we found an isometry with bliss and checked it in integers. Both are isometric to the class with $|"Aut"| = "1,036,800"$, the only one with this order among the 163 classes of the genus $"II"_16 (3^(+8))$ [P1, App. B.4, Table B4] (number 158, counted from 0, among the classes of the entry of `genera_F` with $p = 3$, $f = 16$, $s = 8$ in the certificate `odd-prime-masses.json` of [P1]).
- Hence $L$ is obtained by gluing two copies of $U$ along their discriminant groups $3^8$, and $U plus.o U'$ has index $3^8$ in $L$. The lattice $A_2 times.o E_8$ is in the same genus, with $|"Aut"| = "4,180,377,600"$, but it is not $U$.

= 5. Four gluings of one lattice of rank 16 <s5>

#policy-box[
*What this section does.* We show that the four ways of gluing two copies of one lattice $F_16$ of rank 16 to an extremal lattice of rank 32 differ only by an orbit of $"Sp"(6,2)$ on the 240 cosets $O(q) slash overline(G)$ (Theorem 5.1), and that this explains their four automorphism groups at once. One of the four is new; the others are two code lattices and a Hurwitz lattice.
]

Let $F_16$ be the component of rank 16 with $(s, delta) = (10, 1)$ [P1, §4.1.3, Table 6]; it is the only class of $cal(M)(16)$ with these invariants. (The letters $F, K$ of [P1, §4] denote the two components of a pair; $F_16$ is the name of one particular lattice.) The lattice $F_16$ has rank 16, discriminant group $2^(10)$ and minimum 4, and $|"Aut" F_16| = "743,178,240" = 2^9 dot "1,451,520"$, in agreement with [P1, Table 6]. The structure of $"Aut" F_16$ is not determined in [P1], and we determine it here #acalc("A5"). $O_2 ("Aut" F_16)$ is extraspecial of order $2^9$ and type $+$, and the quotient is isomorphic to $"Sp"(6, 2)$: GAP's `IsomorphismGroups` returned an isomorphism onto $"Sp"(6,2)$, which we checked to be a bijective homomorphism. Hence $"Aut" F_16 = 2^(1+8)_+ . "Sp"(6, 2)$. Among the pairs of components of type $(16, 16)$ in [P1, §4.1.3] there is exactly one with $F tilde.equiv K tilde.equiv F_16$. Its gluings form exactly four double cosets, all extremal; the list is complete because the masses add up to the right side of [P1, (4.1)] [P1, C4].


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
  block(width: 440pt, height: 168pt, breakable: false, {
    place(dx: 0pt, dy: 0pt, block(width: 440pt, text(size: 8pt)[the 240 points $O(q) slash overline(G)$, where $overline(G) = "Aut" F_16 slash \{plus.minus 1\}$\; the subgroup $2^8$ of $overline(G)$ acts trivially, and $"Sp"(6,2)$ has four orbits (drawn to scale)]))
    let seg = ((0, 1, luma(40)), (1, 57, rgb("#7fa7d6")), (57, 120, rgb("#a9cba0")), (120, 240, rgb("#e3c98f")))
    for (a, b, c) in seg { place(dx: a * 440pt / 240, dy: 24pt, rect(width: (b - a) * 440pt / 240, height: 14pt, fill: c, stroke: 0.4pt + white)) }
    place(dx: -2pt, dy: 39pt, polygon(fill: luma(40), (2pt, 0pt), (-1pt, 5pt), (5pt, 5pt)))
    let mids = (0.5, 29, 88.5, 180)
    for (i, m) in mids.enumerate() {
      place(line(start: (m * 440pt / 240, 38pt), end: ((52 + 112 * i) * 1pt, 70pt), stroke: 0.5pt + luma(90)))
    }
    let boxes = (
      ([$ell = 1$], [$"Sp"(6,2)$], dl("known:code:walk908"), [code lattice], [$"1,486,356,480"$]),
      ([$ell = 56$], [$"U"_4 (2)$], dl("16-16-s10d1:0"), [new], [$"26,542,080"$]),
      ([$ell = 63$], [$2^5 : S_6$], dl("known:kmax24:9"), [code lattice], [$"23,592,960"$]),
      ([$ell = 120$], [$G_2(2)$], dl("known:hurwitz:s1-22"), [Hurwitz lattice], [$"12,386,304"$]),
    )
    for (i, (l, st, lab, prev, c)) in boxes.enumerate() {
      place(dx: (112 * i) * 1pt, dy: 70pt, box(width: 104pt, height: 92pt, radius: 2pt, inset: 3pt, stroke: 0.5pt + luma(120),
        align(center + horizon, text(size: 7.5pt)[*orbit length* #l \ point stabilizer #st \ #lab \ (#prev) \ $|C(sigma)| = 2 |"Aut" F_16| slash ell$ \ $= #c$])))
    }
  }),
  caption: [Figure 1. Theorem 5.1 in a picture. Each orbit of $"Sp"(6,2)$ on the 240 points is one double coset of gluings of $F_16$ with itself, that is, one extremal lattice $L = F_16 plus.o_phi F_16$ with an involution $sigma$ such that $L^sigma tilde.equiv F_16$. A short orbit means a large centralizer: $|C(sigma)| = 2 |"Aut" F_16| slash ell$, and the image of $C(sigma)$ in $"Sp"(6,2)$ is the stabilizer of a point of the orbit. The four values of $1 slash |C(sigma)|$ add up to the right side of [P1, (4.1)] because the orbit lengths add up to 240.],
)

#thm[
*Theorem 5.1.* Let $overline(G) = "Aut" F_16 slash {plus.minus 1} tilde.equiv 2^8 . "Sp"(6,2)$, acting on the discriminant form $q$ of $F_16$. Then $[O(q) : overline(G)] = 240$, the normal subgroup $2^8$ acts trivially on $O(q) slash overline(G)$ (with $overline(G)$ acting by left multiplication), and $"Sp"(6,2)$ has four orbits on these 240 points, of lengths $ell = 1, 56, 63, 120$. The four extremal lattices $L = F_16 plus.o_phi F_16$ correspond to the four orbits, and for the involution $sigma$ with $L^sigma tilde.equiv F_16$,
$ |C(sigma)| = 2 |"Aut" F_16| slash ell . $
The image of $C(sigma)$ in $"Sp"(6,2)$ is the stabilizer of a point of the orbit: $"Sp"(6,2)$, $"U"_4 (2)$, $2^5 : S_6$, $G_2(2)$.
]
_Strategy of the proof._ We turn the double cosets of [P1, Prop. 4.1(ii)] for the pair $(F_16, F_16)$ into the orbits of a group on 240 points, the cosets of $overline(G)$ in $O(q)$. Up to a kernel of order 4, the centralizer $C(sigma)$ is the stabilizer of a point, which gives $|C(sigma)| = 2 |"Aut" F_16| slash ell$. The list of [P1] has four double cosets, and their centralizers give orbit lengths that add up to 240, so these are all the orbits. Finally, the normal subgroup $2^8$ fixes a point of each orbit, hence every point, so the action is one of $"Sp"(6,2)$.
#proofbox[
We use the notation of [P1, Prop. 4.1] with $h = sigma$, $F tilde.equiv K tilde.equiv F_16$, and $G_F = "Aut" F$, $G_K = "Aut" K$. The kernel of $"Aut" F_16 -> O(q)$ is $\{plus.minus 1\}$ [P1, C4] (pair `16-16-s10d1` of the certificate `involution-gluing-masses.json`), so $overline(G)_F tilde.equiv "Aut" F_16 slash \{plus.minus 1\} = overline(G)$, of order $"371,589,120"$, and likewise for $K$. Since $|O(q)| = "89,181,388,800"$ [P1, C4] (recorded for the same pair in that certificate), we get $[O(q) : overline(G)] = 240$.

_Double cosets as orbits._ Fix $phi_0 in "Iso"(q_F, -q_K)$. Then $phi |-> psi = phi_0^(-1) phi$ is a bijection $"Iso"(q_F, -q_K) -> O(q_F)$, and $overline(v) phi overline(u)^(-1)$ corresponds to $(phi_0^(-1) overline(v) phi_0) psi overline(u)^(-1)$. So the double cosets of [P1, Prop. 4.1(ii)] correspond to the double cosets $cal(H) backslash O(q) slash overline(G)$ with $cal(H) = phi_0^(-1) overline(G)_K phi_0$, a subgroup of $O(q)$ of order $|overline(G)|$. These in turn correspond to the orbits of $cal(H)$ on the 240 cosets $O(q) slash overline(G)$. The stabilizer of $psi overline(G)$ in $cal(H)$ is $cal(H) inter psi overline(G) psi^(-1)$, of order $|cal(H)| slash ell$, where $ell$ is the length of its orbit.

_The stabilizer of $phi$._ By [P1, Prop. 4.1(iii)], $C(sigma) tilde.equiv "Stab"(phi)$. The map $"Stab"(phi) -> overline(G)_F$, $(u, v) |-> overline(u)$, has kernel $\{plus.minus 1\} times \{plus.minus 1\}$ of order 4. Indeed, $overline(u) = 1$ forces $overline(v) = 1$, and the kernels of $"Aut" F -> overline(G)_F$ and $"Aut" K -> overline(G)_K$ are $\{plus.minus 1\}$. Its image is $\{overline(u) in overline(G) : phi overline(u) phi^(-1) in overline(G)_K\}$. Since $phi overline(u) phi^(-1) = phi_0 (psi overline(u) psi^(-1)) phi_0^(-1)$, the condition is $psi overline(u) psi^(-1) in cal(H)$, so the image is $overline(G) inter psi^(-1) cal(H) psi = psi^(-1) (cal(H) inter psi overline(G) psi^(-1)) psi$, the stabilizer above conjugated by $psi$. In particular it has order $|cal(H)| slash ell$. Hence
$ |C(sigma)| = 4 |cal(H)| slash ell = 4 |overline(G)| slash ell = 2 |"Aut" F_16| slash ell . $

_The four orbits._ By [P1, C4] the double cosets are exactly the four classes of Table 8. Their values of $|C(sigma)|$ give $ell = 1, 56, 63, 120$, with sum 240, so these are all the orbits. The orbit of length 1 is that of #dl("known:code:walk908")\; if $psi_1$ corresponds to it, then $cal(H) subset.eq psi_1 overline(G) psi_1^(-1)$, and since the two groups have the same order, $cal(H) = psi_1 overline(G) psi_1^(-1)$. Replacing $phi_0$ by $phi_0 psi_1$ we may assume $cal(H) = overline(G)$, so the orbits are those of $overline(G)$ on $O(q) slash overline(G)$, and the orbit of $phi$ is that of $psi overline(G)$. For each $phi$ the image of $C(sigma)$ in $overline(G)$ is then $overline(G) inter psi^(-1) overline(G) psi$, the stabilizer of the point $psi^(-1) overline(G)$ (not of $psi overline(G)$). It has index $ell$ in $overline(G)$, so the orbit of $psi^(-1) overline(G)$ also has length $ell$; since the four lengths $1, 56, 63, 120$ are distinct, it is the orbit of $psi overline(G)$. So the image of $C(sigma)$ in $overline(G)$ is the stabilizer of a point of the orbit of $phi$. Next, $O_2(overline(G)) = 2^8$ is normal in $overline(G)$. It is the image of $O_2("Aut" F) = 2^(1+8)_+$, which contains $-1$, because this image is a normal 2-subgroup of order $2^8$ and the quotient $"Sp"(6,2)$ is simple. It lies in the stabilizer of a point of each of the four orbits: for #dl("known:code:walk908") the stabilizer is all of $overline(G)$, and for the other three lattices the image of $C(sigma)$ in $"Aut" F$ contains $O_2("Aut" F)$ #acalc("A5"), so the image of $C(sigma)$ in $overline(G)$, a point stabilizer of the orbit by the above, contains $2^8$. If a normal subgroup $N$ of $overline(G)$ fixes a point $x$, then $N = g N g^(-1)$ fixes $g x$ for every $g in overline(G)$. So $2^8$ fixes every point of the four orbits, that is, all 240 points; it acts trivially and the action factors through $overline(G) slash 2^8 tilde.equiv "Sp"(6, 2)$. Since $2^8$ lies in every point stabilizer, the stabilizer of a point in $"Sp"(6,2)$ is the image of its stabilizer in $overline(G)$, so the image of $C(sigma)$ in $"Sp"(6,2)$ is the stabilizer of a point of the orbit of $phi$, as stated. These images were identified by explicit isomorphisms #acalc("A5") as $"Sp"(6,2)$, $"U"_4 (2)$, $2^5 : S_6$ and $G_2(2)$, of indices $1, 56, 63, 120$.

As a check, $sum 1 slash |C(sigma)| = sum ell slash (4 |overline(G)|) = 240 slash "1,486,356,480"$, which equals the right side $|O(q)| slash |"Aut" F_16|^2$ of [P1, (4.1)]. #qed
]

#ptable(
  columns: 5, align: (left, left, right, left, right),
  caption: [Table 8. The four gluings of two copies of $F_16$.],
  table.header([lattice], [previously], [$|C(sigma)|$], [image in $"Sp"(6,2)$], [$ell$]),
  [#dl("known:code:walk908")], [code lattice [KV89]], [$"1,486,356,480"$], [$"Sp"(6,2)$], [1],
  [#dl("16-16-s10d1:0")], [new], [$"26,542,080"$], [$"U"_4 (2) tilde.equiv "PSp"(4,3)$], [56],
  [#dl("known:kmax24:9")], [code lattice [KV89], [KN93]], [$"23,592,960"$], [$2^5 : S_6$], [63],
  [#dl("known:hurwitz:s1-22")], [Hurwitz [BN97]], [$"12,386,304"$], [$G_2(2)$], [120],
)

*The new lattice.* For $L$ = #dl("16-16-s10d1:0") the class of $sigma$ is $\{sigma, -sigma\}$, so $|G : C(sigma)| = 2$, and the elements outside $C(sigma)$ interchange $F$ and $K$. The group is $G = O_2(G) . "GO"_6^-(2)$. Here $F^*(G) = O_2(G)$ has order $2^(10)$, centre of order 4 and Frattini subgroup of order 2, and $G slash O_2(G)$ is isomorphic to the full orthogonal group $"GO"_6^-(2) tilde.equiv W(E_6)$ of order $"51,840"$ (explicit isomorphism #acalc("A6")). Its chief factors are $2, "U"_4 (2), 2^8, 2, 2$ from the top, and $G slash O_2(G)$ acts faithfully and absolutely irreducibly on the factor $2^8$. The presentation is proved by the radical-plus-quotient case of Lemma 2.3. The minimal vectors form five orbits, of lengths $"1,728"$, $"2,304"$, $"34,560"$, $"34,560"$ and $"73,728"$.

*Remark.* The stabilizers $2^5 : S_6$ (index 63) and $"U"_3 (3) : 2 = G_2(2)$ (index 120) are maximal subgroups of $"Sp"(6,2) = S_6(2)$. The stabilizer $"U"_4 (2)$ is the derived group, of index 2, of the maximal subgroup $"U"_4 (2) : 2 tilde.equiv W(E_6)$ of index 28 [ATLAS, $S_6(2)$]; we also checked the list of maximal subgroups against the GAP character table library #acalc("A16"). The 240 points are the cosets $O(q) slash overline(G)$. They should not be confused with the action of $"Sp"(6,2)$ on the 240 cosets of $"U"_3 (3)$, which is imprimitive.

= 6. Elements of type $3 - (9, 14) - 7$ <s6>

#policy-box[
*What this section does.* An element of order 3 of type $3 - (9, 14) - 7$ (notation of [P1, §3.3]) has a fixed lattice of rank 14 and determinant $3^7$. We show that exactly two extremal lattices have such an element, that they are 3-neighbours, and that in the larger one the type of an element of $Z(O_3(G))$ is read off from an invariant quadratic form.
]

#thm[
*Theorem 6.1.* Exactly two extremal even unimodular lattices of rank 32 have an automorphism of order 3 of type $3 - (9, 14) - 7$: #dl("9-14-7:0-0/0"), with $F^*("Aut" L) = 2 times 3^(1+4)_-$ and $"Aut" L slash F^* tilde.equiv "GL"(2,3)$, and #dl("12-8-8:0-50/2") with $|"Aut" L| = 2^7 3^8$ and $O_3(G)$ special of order $3^(3+4)$. They are 3-neighbours.
]
_Strategy of the proof._ The proof is identity [P1, (4.1)] for one pair of components. Each of the two genera contains exactly one lattice of minimum at least 4, so there is one pair $(F, Z)$, and the right side of (4.1) is computed from $|O(q_F)|$, $|"Aut" F|$ and $|"U"(Z)|$. Two lattices known to have such an element contribute two terms that already add up to the right side; as every term is positive, there are no other double cosets. The 3-neighbour relation is then exhibited by a reflection of the discriminant group.
#proofbox[
Let $L$ be an extremal even unimodular lattice of rank 32 and $h in "Aut" L$ of type $3 - (9, 14) - 7$. As in [P1, §4.1.1] put $F = L^h$ and $Z = F^perp inter L$, a Hermitian lattice over $ZZ[omega]$ with $omega$ acting as $h$, and let $G_F = "Aut" F$ and $G_Z = "U"(Z)$, the centralizer of $h|_Z$ in $"Aut" Z$. Then $F$ lies in the genus $(14, 7)$ of even 3-elementary lattices of rank 14 with $F^\# slash F tilde.equiv FF_3^7$, $Z$ lies in the genus $(9, 7)$ of Hermitian lattices, and both have minimum at least 4. By [P1, App. B.4, Table B4] each of these two genera has exactly one class of minimum at least 4. So there is exactly one pair of components $(F, Z)$, and by [P1, Prop. 4.1(ii)] the isometry classes of pairs $(L', h)$ with these components, $L'$ even unimodular of rank 32, correspond to the double cosets $overline(G)_Z backslash "Iso"(q_F, -q_Z) slash overline(G)_F$. For these components $|"Aut" F| = "8,491,392" = 2^7 3^6 7 dot 13$ #acalc("A5"), and $|"U"(Z)| = "22,674,816" = 2^7 3^(11)$, computed in [P1, C8] and recorded in its certificate `odd-prime-masses.json` ([P1, §7.4]); it is half of $|"Aut" Z| = "45,349,632"$ #acalc("A5").

_The right side of [P1, (4.1)]._ Write $q_F (a) = 2 Q(a) slash 3$ with $Q(a) in FF_3$ [P1, Table 3]. Then $Q$ is a quadratic form on $F^\# slash F tilde.equiv FF_3^7$ with polar form $3 b_F mod 3$, which is non-degenerate because the discriminant bilinear form $b_F$ is. Hence $O(q_F) = O(Q)$ is the full orthogonal group $"GO"_7(3)$. For odd $q$, the full orthogonal group of a non-degenerate quadratic form on $FF_q^(2m + 1)$ has order $2 q^(m^2) product_(i = 1)^m (q^(2i) - 1)$; the two isometry classes of such forms differ by a non-square scalar and have the same orthogonal group. With $q = 3$ and $m = 3$,
$ |O(q_F)| = 2 dot 3^9 (3^2 - 1)(3^4 - 1)(3^6 - 1) = 2 dot "19,683" dot 8 dot 80 dot 728 = "18,341,406,720" . $
The set $"Iso"(q_F, -q_Z)$ is not empty, since the lattices below exist, so $|"Iso"(q_F, -q_Z)| = |O(q_F)|$ and the right side of [P1, (4.1)] is
$ ("18,341,406,720") / ("8,491,392" dot "22,674,816") = 5 / ("52,488") . $

_Two pairs exhaust it._ In #dl("9-14-7:0-0/0") the centre of $O_3(G)$ has order 3. A generator $h$ is of type $3 - (9, 14) - 7$, its class is $\{h, h^(-1)\}$, and $|C(h)| = "11,664"$. In #dl("12-8-8:0-50/2") the 8 isotropic elements of $Z(O_3(G))$ (Proposition 6.2) are of this type. They form one conjugacy class, which contains $h^(-1)$ with $h$, and $|C(h)| = "104,976"$ (Table 8a). By [P1, Prop. 4.1(iii)] each pair $(L, h)$ contributes $1 slash |C(h)|$ to the left side of [P1, (4.1)], and the two pairs are not isometric, since the two lattices have different $|"Aut" L|$. Since
$ 1 / "11,664" + 1 / "104,976" = 10 / "104,976" = 5 / "52,488" , $
these two terms already make up the right side. Every term of the left side is positive, so there are exactly two double cosets, and every pair with the components $F$ and $Z$ is isometric to one of the two pairs above. In particular both double cosets are extremal. Every extremal lattice with an automorphism of this type has these components, so it is one of the two lattices, and they are distinct.

The certificate of [P1, C8] lists this pair term by term: $|O(q_F)| = "18,341,406,720"$, $|"Aut" F| = "8,491,392"$, $|"U"(Z)| = "22,674,816"$, two double cosets with $|"Stab"(phi)| = "11,664"$ and $"104,976"$, both extremal, and the extremal mass $5 slash "52,488"$. These agree with the values above; of them the argument uses only $|"U"(Z)|$.

_3-neighbours._ The two lattices are $L(phi_0)$ and $L(phi_1)$ for representatives $phi_0, phi_1$ of the two double cosets, glued along the discriminant groups $D_F tilde.equiv FF_3^7$. Let $r$ be the reflection of $D_F$ in a vector $v$ with $q(v) != 0$, so that $r$ fixes $v^perp$ pointwise, and put $phi' = phi_0 compose r$. Then $L(phi_0)$ and $L(phi')$ both contain $F plus.o Z$ and the gluing over $v^perp$, of order $3^6$, so their intersection has index 3 in each. The double coset of $phi'$ is determined by the order of $overline(G)_F inter phi'^(-1) overline(G)_Z phi'$, which is $"1,296"$ for $phi_0$ and $"11,664"$ for $phi_1$. In general the image of $C(h) tilde.equiv "Stab"(phi)$ in $"Aut" F$ is the preimage of $overline(G)_F inter phi^(-1) overline(G)_Z phi$ under $"Aut" F -> O(q_F)$. Here this map is injective (the certificate of [P1, C8] records a trivial kernel), so $"1,296"$ and $"11,664"$ are also the orders of the images of $C(h)$ in $"Aut" F$ for the two lattices, $"8,491,392" slash "6,552"$ and $"8,491,392" slash 728$ (Table 8a). Of the 729 reflections, 72 move $phi_0$ into the double coset of $phi_1$, for instance the reflection in $v = (2, 1, 0, 1, 0, 0, 0)$ with $q(v) = 2$ #acalc("A7"). #qed
]


#figure(kind: image, numbering: none, gap: 0.8em,
  _x-diagram(440, 222, h: 40, (
    f: (0, 0, 200, "plain", [$F$: the only lattice of minimum $>= 4$ in the genus $(14, 7)$ \ rank 14, $F^\# slash F tilde.equiv FF_3^7$, $|"Aut" F| = "8,491,392"$]),
    z: (240, 0, 200, "plain", [$Z$: the only lattice of minimum $>= 4$ in the Hermitian genus $(9, 7)$ \ rank 9 over $ZZ[omega]$, $|"U"(Z)| = "22,674,816"$]),
    id: (90, 62, 260, "plain", [right side of [P1, (4.1)]: $|O(q_F)| slash (|"Aut" F| |"U"(Z)|) = 5 slash "52,488"$ \ $= 1 slash "11,664" + 1 slash "104,976"$: exactly two double cosets]),
    l0: (0, 124, 170, "result", [$L(phi_0)$ = #dl("9-14-7:0-0/0") \ $|C(h)| = "11,664"$]),
    l1: (270, 124, 170, "result", [$L(phi_0 compose r) tilde.equiv L(phi_1)$ \ = #dl("12-8-8:0-50/2") \ $|C(h)| = "104,976"$]),
    nb: (172, 124, 96, "note", [3-neighbours: \ $phi_0 compose r$ for 72 of the \ 729 reflections $r$ of $D_F$]),
    cap: (110, 182, 220, "plain", [$L(phi_0) inter L(phi_0 compose r)$: index 3 in each\; \ it contains $F plus.o Z$ and the gluing over $v^perp$]),
  ), (
    ("f", "b", "id", ("t", 0.25)), ("z", "b", "id", ("t", 0.75)),
    ("id", ("b", 0.2), "l0", "t"), ("id", ("b", 0.8), "l1", "t"),
    ("cap", "l", "l0", "b"), ("cap", "r", "l1", "b"),
  )),
  caption: [Figure 2. Theorem 6.1 in a picture. The two components are unique (Table B4 of [P1]), so the right side of [P1, (4.1)] is one number; the two lattices already fill it. A reflection $r$ of the discriminant group $D_F$ in a vector $v$ with $q(v) != 0$ changes the gluing on the line of $v$ only, so $L(phi_0)$ and $L(phi_0 compose r)$ share a sublattice of index 3; for 72 of the 729 reflections $L(phi_0 compose r) tilde.equiv L(phi_1)$.],
)

#thm[
*Proposition 6.2.* For $L$ = #dl("12-8-8:0-50/2") the group $G$ acts irreducibly on $cal(Z) = Z(O_3(G)) tilde.equiv FF_3^3$ and preserves a non-degenerate quadratic form $Q_3$, unique up to scalars; we normalize it so that $Q_3(w) = 1$ for the six elements whose fixed lattice is $sqrt(3) E_8$. A non-zero $w in cal(Z)$ is of type $3 - (9, 14) - 7$ if $Q_3(w) = 0$ (8 elements), $3 - (12, 8) - 8$ with fixed lattice $sqrt(3) E_8$ if $Q_3(w) = 1$ (6), and $3 - (6, 20) - 6$ if $Q_3(w) = 2$ (12) #acalc("A8").
]

#ptable(
  columns: 3, align: (left, left, left),
  caption: [Table 8a. The two lattices with an element $h$ of type $3 - (9, 14) - 7$ #acalc("A2", "A5"). Chief factors are listed from the top; $*$ marks a non-split factor. $F$ is the fixed lattice of $h$ (rank 14) and $Z$ the sublattice of rank 18 on which $h^2 + h + 1 = 0$.],
  table.header([], [#dl("9-14-7:0-0/0")], [#dl("12-8-8:0-50/2")]),
  [$|"Aut" L|$], [$"23,328" = 2^5 3^6$], [$"839,808" = 2^7 3^8$],
  [$F^*(G)$], [$chevron.l -1 chevron.r times 3^(1+4)_-$, exponent 9], [$chevron.l -1 chevron.r times 3^(3+4)$, exponent 3],
  [$G slash F^*$], [$"GL"(2,3)$, order 48], [order 192, $O_2 = 2^(1+4)_+$, quotient $S_3$],
  [chief factors], [$2, 3, 2^2, 2^*, 2, 3^2, 3, 3, 3^*$], [$2, 3, 2^2, 2^2, 2^*, 2^*, 3^4, 3^(3 *)$],
  [classes], [140], [107],
  [orbits on minimal vectors], [35, of lengths $108$ to $"23,328"$], [8: $216, 864, "1,944", "3,888",$ #linebreak() $"17,496", "17,496", "34,992", "69,984"$],
  [$|C(h)|$], [$"11,664"$], [$"104,976"$],
  [index of the image of $C(h)$ in $"Aut" F$], [$"6,552" = 2^3 3^2 7 dot 13$], [$728 = 2^3 7 dot 13$],
  [index of the image of $C(h)$ in $"U"(Z)$], [$"1,944"$], [$216$],
)

In both lattices $F$ is the same lattice and $Z$ is the same Hermitian lattice, as the uniqueness of the pair $(F, Z)$ requires. The automorphism group of $Z$ as a $ZZ$-lattice has order $"45,349,632" = 2 |"U"(Z)|$, the factor 2 coming from a semilinear element. $C(h)$ acts faithfully on $Z$, and its image lies in $"U"(Z)$ since it commutes with $h$. The kernel of $C(h) -> "Aut" F$ has order 9: it is generated by $h$ and a second element of order 3 acting only on $Z$. The primes 7 and 13 of $|"Aut" F|$ do not divide $|C(h)|$. The image of $C(h)$ in $"Aut" F times "U"(Z)$ consists of the pairs that preserve the gluing [P1, Prop. 4.1(iii)].

How $G slash F^* = 2^(1+4)_+ . S_3$ of #dl("12-8-8:0-50/2") maps to the orthogonal group of $(cal(Z), Q_3)$, and how the three 3-groups $3^(1+4)_-$, $3^(3+4)$ and $3^4$ that occur over $O_3 slash Phi(O_3) = 3^4$ in #dl("9-14-7:0-0/0"), #dl("12-8-8:0-50/2") and #dl("10-12-10:0-1/358") are related, we leave open.

= 7. The golden ratio <s7>

#policy-box[
*What this section does.* A fixed-point-free automorphism of order 5 makes a lattice a module over $ZZ[phi]$. We find the six such lattices among the 99, say for which of them the whole group is $ZZ[phi]$-linear, and describe the icosian part of one of them.
]

#thm[
*Lemma 7.1.* If $tau in "Aut" L$ has order 5 and no non-zero fixed vector, then $s = tau + tau^4 - tau^2 - tau^3$ and $Psi = -(tau^2 + tau^3)$ are integral, $s^2 = 5$, $Psi = (1 + s) slash 2$, $Psi^2 = Psi + 1$, and $s$ is self-adjoint with $(s x, s x) = 5 (x, x)$. Hence $L$ is a $ZZ[phi]$-module of rank 16, and $H(x, y) = ((x, y) + (x, s y) slash sqrt(5)) slash 2$ is a $ZZ[phi]$-bilinear form with values in $frak(d)^(-1)$ and trace $(x, y)$.
]
#proofbox[
Since $tau$ has no non-zero fixed vector, $tau - 1$ is invertible on $L times.o QQ$, and $(tau - 1)(tau^4 + tau^3 + tau^2 + tau + 1) = tau^5 - 1 = 0$ gives $Phi_5(tau) = 0$, where $Phi_5(x) = x^4 + x^3 + x^2 + x + 1$. As $Phi_5$ is irreducible, $ZZ[tau] tilde.equiv ZZ[zeta_5]$. The identities $(zeta + zeta^4 - zeta^2 - zeta^3)^2 = 5$ and $-(zeta^2 + zeta^3) = (1 + sqrt(5)) slash 2$ in $ZZ[zeta_5]$ therefore hold for $tau$, and $s$ and $Psi$ are integral matrices because $tau$ is. Since $tau$ preserves the form, its adjoint is $tau^(-1) = tau^4$. So the adjoint of $s$ is $tau^4 + tau - tau^3 - tau^2 = s$, and $(s x, s x) = (x, s^2 x) = 5 (x, x)$. Put $a = (x, y)$ and $b = (x, s y)$; then $(x, s^2 y) = 5 a$. With $Psi = (1 + s) slash 2$ this gives $H(x, Psi y) = ((a + b) slash 2 + (b + 5 a) slash (2 sqrt(5))) slash 2 = phi (a + b slash sqrt(5)) slash 2 = phi H(x, y)$. Since $Psi$ is self-adjoint, also $H(Psi x, y) = phi H(x, y)$, so $H$ is $ZZ[phi]$-bilinear. Its trace over $QQ$, the sum of $H$ and its conjugate under $sqrt(5) |-> -sqrt(5)$, is $(x, y)$. Finally $(x, s y) = 2 (x, Psi y) - (x, y)$ is an integer congruent to $(x, y)$ modulo 2, so $2 H(x, y) = (x, y) + (x, s y) slash sqrt(5)$ lies in $ZZ + (1 slash sqrt(5)) ZZ$ with both coefficients of the same parity, and $H(x, y) in (1 slash sqrt(5)) ZZ[phi] = frak(d)^(-1)$. #qed
]

#ptable(
  columns: 4, align: (left, right, right, left),
  caption: [Table 9. The six lattices among the 99 with a fixed-point-free automorphism of order 5 #acalc("A10"). The third column is the size of the class of the element $tau$ chosen for Lemma 7.1. Some lattices have several such classes: #dl("known:hurwitz:s1-21") has three, of sizes $"13,824"$, $"1,728"$ and $96$, and #dl("known:eis:s1-440") has two, both of size 2.],
  table.header([lattice], [$|"Aut" L|$], [class of $tau$], [action of $G$ on $sqrt(5)$]),
  [#dl("12-8-8:0-49/5")], [$2^7 3^4 5$], [24], [depends on $tau$],
  [#dl("5-8-0-0:8")], [$2^(10) 3^4 5^2$], [288], [depends on $tau$],
  [#dl("8-24:1215")], [$2^(10) 3 dot 5$], [4], [index-2 subgroup linear, the rest $sqrt(5) |-> -sqrt(5)$],
  [#dl("known:eis:s1-440")], [$2^8 3^2 5$], [2], [all of $G$ is $ZZ[phi]$-linear],
  [#dl("known:hurwitz:s1-21")], [$2^(13) 3^5 5^4$], [13,824], [depends on $tau$],
  [#dl("known:code:r32")], [$2^(31) 3^5 5^2 7 dot 17 dot 31$], [$"2,901,147,648"$], [depends on $tau$],
)

We found these lattices among the conjugacy classes in the profiles. For $"BW"_32$ = #dl("known:code:r32") we used the classes of #acalc("A4"): the class of frame shape $5^8 1^(-8)$ has no fixed vectors, $"2,901,147,648"$ elements and centralizer of order $"16,588,800"$. For each lattice we took one $tau$ and checked in integers that $s$ and $Psi$ are integral, $s^2 = 5$, $s$ is self-adjoint and $Psi^2 = Psi + 1$ #acalc("A10"). The last column of Table 9 compares $g^(-1) s g$ with $plus.minus s$ for the generators $g$ of $G$ returned by bliss.

- For #dl("known:eis:s1-440") the class of $tau$ has 2 elements. So the number of Sylow 5-subgroups is at most 2 and congruent to 1 modulo 5, hence 1, and $O_5(G) = chevron.l tau chevron.r$ is normal. Every generator commutes with $s$, so $G$ maps $tau$ to $tau^(plus.minus 1)$ only (an element mapping $tau$ to $tau^(plus.minus 2)$ would map $s$ to $-s$). Thus all of $G$ is $ZZ[phi]$-linear, and $"Aut" L$ is the group of $ZZ[phi]$-linear isometries of a $ZZ[phi]$-lattice of rank 16.
- For #dl("8-24:1215") the class of $tau$ has 4 elements, and $O_5(G) = chevron.l tau chevron.r$ is normal by the same argument. The signs $g |-> plus.minus 1$ in $g^(-1) s g = plus.minus s$ define a non-trivial homomorphism $G -> \{plus.minus 1\}$. Its kernel, of index 2, is $ZZ[phi]$-linear, and the other elements act on $ZZ[phi]$ by the Galois automorphism $sqrt(5) |-> -sqrt(5)$.
- For the other four lattices the class of $tau$ is large, and $s$ depends on the choice of $tau$; the group permutes the resulting $ZZ[phi]$-structures.

For the three previously known lattices of Table 9, the $ZZ[phi]$-structure follows at once from the existence of $tau$. We did not find it discussed in the sources of Table 5. The section is restricted to fixed-point-free elements of order 5: if $tau$ has fixed vectors, $s$ vanishes on them and is not a square root of 5.

*The icosian part of #dl("12-8-8:0-49/5").* Here $|"Aut" L| = 2^7 3^4 5$, $F^*(G) = 2 times 3^(1+2) compose E$, and the layer is $E tilde.equiv "SL"(2, 5)$ (explicit isomorphism #acalc("A6")). The character of $E$ on $QQ^32 times.o CC$ is $2 chi_2 + 2 chi'_2 + 6 chi_4$, where $chi_2$ and $chi'_2$ are the two Galois-conjugate characters of degree 2 with values in $QQ(sqrt(5))$ and $chi_4$ is the rational character of degree 4 #acalc("A10"). So $QQ^32 = V_8 plus.o V_24$, with $V_8$ the sum of the isotypic components of $chi_2$ and $chi'_2$, and $V_24$ that of $chi_4$. Since there are no non-zero $E$-homomorphisms between different $QQ$-irreducible constituents, $"End"_E (QQ^32) = "End"_E (V_8) times "End"_E (V_24)$. Each factor is a full matrix ring over the endomorphism division algebra of a $QQ$-irreducible representation, whose centre is the field of character values: $QQ(sqrt(5))$ for $chi_2 + chi'_2$ and $QQ$ for $chi_4$. So the centre of $"End"_E (QQ^32)$ is $QQ(sqrt(5)) times QQ$, and no element of the centre acts as $sqrt(5)$ on all of $QQ^32$. We cut out $L inter V_8$ with the central idempotent $e = (chi_2(1) slash |E|) sum_(g in E) (chi_2 + chi'_2)(g^(-1)) g$, the sum of the central idempotents of $chi_2$ and $chi'_2$ (both of degree 2). Its Gram matrix $G_8$ satisfies: $G_8 slash 3$ is integral, even, of determinant 1. Since $E_8$ is the only even unimodular lattice of rank 8, $L inter V_8 tilde.equiv sqrt(3) E_8$, with minimum 6 #acalc("A10"). By determinants, $(L inter V_8) plus.o (L inter V_24)$ has index $3^8$ in $L$, and $L inter V_24$ has determinant $3^8$ and minimum 4.

The group $E tilde.equiv 2.A_5$ is the abstract group of the 120 icosians, whose finite sums form the icosian ring $cal(I)$ [SPLAG, Ch. 8, §2.1, p. 207]. With the Euclidean norm $"EN"(a + b sqrt(5)) = a + b$ applied to the quaternion norm, $cal(I)$ is isometric to $E_8$ scaled to minimum 1 [SPLAG, Ch. 8, (2)–(3), p. 208]. With $2 "EN"$ it is $E_8$ itself, and with $6 "EN"$ it is $sqrt(3) E_8$. So $L inter V_8$ is isometric to $(cal(I), 6 "EN")$ as a lattice. We have not checked whether $E$ acts on it as the icosians act on $cal(I)$ by multiplication.

= 8. Origins of extremality <s8>

#policy-box[
*What this section does.* An extremal lattice of rank 32 is one without roots. For lattices obtained by gluing two components, we separate the case in which the absence of roots follows from the components alone (Lemma 8.1) from the other cases, and classify the 99 lattices by their decompositions by involutions and by elements of order 3 (Theorem 8.2). On the 83 pairs of components of involution type of [P1] the criterion of Lemma 8.1 coincides with "every gluing is extremal".
]

For a lattice $M$ put $m_M = min{(u, u) : u in M^\#, u in.not M}$, the smallest norm of a non-zero class of $M^\# slash M$.

#thm[
*Lemma 8.1.* Let $L$ be an even unimodular lattice with primitive sublattices $F$ and $K = F^perp$, $min F >= 4$, $min K >= 4$. If $m_F + m_K > 2$, then $min L >= 4$, and the same holds for every even unimodular lattice obtained by gluing $F$ and $K$.
]
#proofbox[
A vector of $L$ is $u + w$ with $u in F^\#$, $w in K^\#$ in glued classes. If the class of $u$ is $0$, then $u in F$ and $w in K$ and $(u + w, u + w) = (u, u) + (w, w)$ is $0$ or at least 4. Otherwise $(u, u) + (w, w) >= m_F + m_K > 2$, and the norm is even, hence at least 4. #qed
]

#thm[
*Theorem 8.2.* Among the 99 lattices, exactly 92 have a decomposition by an involution or by an element of order 3 with $m_F + m_K > 2$. For each of the other 7, every decomposition by a conjugacy class of involutions $sigma != plus.minus 1$ or of elements of order 3 with non-zero fixed vectors has $m_F + m_K = 2$ #acalc("A11").
]

#ptable(
  columns: 4, align: (left, left, center, left),
  caption: [Table 10. Origins of extremality: for each of the 92 lattices of the first kind the decomposition that the search of the Method found first, and for the 7 others the decompositions examined #acalc("A11"). The two lattices of type $3 - (9, 14) - 7$ also have a decomposition by an element of order 3 with $4 slash 3 + 4 slash 3$, and #dl("5-8-0-0:8") the decomposition of §4.6 into two lattices of rank 16 of the genus $"II"_16 (3^(+8))$ with $4 slash 3 + 4 slash 3$, which is not given by an element.],
  table.header([decomposition], [lattices], [$m_F + m_K$], [origin]),
  [an involution of type $(8, 24)$: in a smallest class for 43 lattices, in another class for 34 (among them #dl("known:kmax24:9"), one of the four gluings of §5)], [77], [$1 + 2$], [components],
  [an involution of type $(12, 20)$ with $m = 3 slash 2$ on the side of rank 20: in a smallest class for 6 lattices (among them the two of type $3 - (9, 14) - 7$ of §6 and #dl("5-8-0-0:8")), in another class for 3], [9], [$1 + 3 slash 2$], [components],
  [an involution of type $(16, 16)$ gluing two copies of $F_16$ (§5): #dl("16-16-s10d1:0"), #dl("known:hurwitz:s1-22"), #dl("known:code:walk908")], [3], [$3 slash 2 + 3 slash 2$], [components],
  [an element of order 3 with fixed lattice $sqrt(3) E_8$: #dl("16-16s14:5423"), #dl("known:zeta17:s1-14")], [2], [$2 slash 3 + 2$], [components],
  [the involution class of size 1,054, type $(16, 16)$, of #dl("known:code:r32") ($"BW"_32$)], [1], [$2 + 2$], [components],
  [every involution and every element of order 3 with non-zero fixed vectors of #dl("16-16-s16d0c:221"), #dl("16-16d8:222190"), #dl("16-16d8:47369"), #dl("16-16d8:68772"), #dl("16-16d8:98086")\; the components of rank 16 have many classes of norm 1], [5], [$1 + 1$], [gluing],
  [the same for #dl("known:zeta17:s1-1"), #dl("known:zeta17:s1-10")], [2], [$1 + 1$, $2 slash 3 + 4 slash 3$], [gluing],
)

*Method.* For each lattice we took decompositions $L supset F plus.o K$ by an involution $sigma != plus.minus 1$ ($F = L^sigma$, $K = L^(-sigma)$) or by an element $h$ of order 3 with non-zero fixed vectors ($F = L^h$, $K = ker(h^2 + h + 1)$). For each, we computed the minimal norm of every non-zero class of $F^\# slash F$ and $K^\# slash K$ that has a representative of norm at most 2, by short vectors of the dual lattices #acalc("A11"). Lemma 8.1 does not depend on the type of the discriminant groups. For an element of order 3 the norms of the classes are multiples of $2 slash 3$, and $m_F + m_K > 2$ again forces $min L >= 4$. The order of the search is as follows.
- First, an involution in a smallest conjugacy class: this gives the origin for 51 lattices.
- Of the remaining lattices, 34 also have an involution of type $(8, 24)$, and the decomposition by it gives $1 + 2$: the classes of $sqrt(2) E_8$ have minimal norm 1, those of the component of rank 24 at least 2.
- Three more lattices have an involution of type $(12, 20)$ with $m = 3 slash 2$ on the side of rank 20.
- Two more, #dl("16-16s14:5423") and #dl("known:zeta17:s1-14"), have an element of order 3 whose fixed lattice is $sqrt(3) E_8$ ($m = 2 slash 3$) and whose complement has $m = 2$, so $2 slash 3 + 2 = 8 slash 3$.
- For #dl("known:code:walk908") and #dl("known:code:r32") we used all classes. For #dl("known:code:walk908") the involution class of type $(16, 16)$, which glues two copies of $F_16$, gives $3 slash 2 + 3 slash 2$. For #dl("known:code:r32") the involution class of size $"1,054"$, of type $(16, 16)$, gives $2 + 2$: both components of rank 16 have no non-zero discriminant class with a representative of norm less than 2.
- For the remaining 7 lattices we examined every conjugacy class of involutions other than $-1$ and every class of elements of order 3 with non-zero fixed vectors. All give $m_F + m_K = 2$: either $1 + 1$ for involutions, or $2 slash 3 + 4 slash 3$ for elements of order 3. For #dl("known:zeta17:s1-1") the only involution is $-1$. The lattice #dl("16-16d8:222190") has no element of order 3.

The theorem is restricted to decompositions by involutions and by elements of order 3. Decompositions by elements of order 5 or 17 were not tried.

*Comparison with [P1].* For the 83 pairs of components of involution type in [P1, §4.1.3], all double cosets are extremal exactly for the 7 pairs with $m_F + m_K > 2$. The other 76 pairs have $m_F + m_K = 2$, and for them the extremal part of the mass is a proper fraction of the mass, $M_"ext" slash M$ between $3 slash 34 approx 0.09$ (the pair with list identifier `12-20-c35`) and $7 slash 8$ (`12-20-c4` and `12-20-c61`) #acalc("A11"). On these pairs the criterion of Lemma 8.1 therefore coincides with "every gluing is extremal". In five of the 7 lattices of the second kind, the components of every involution have discriminant classes with a representative of norm 1 (between 16 and 240 such classes for the involutions of type $(16, 16)$); the other two are the lattices of [Jür15]. Such a lattice is extremal because its gluing map does not join two classes of norm 1. In the pairs of components of [P1] from which these lattices come, only 6.5 and 20.1 per cent of the double cosets are extremal (22.6 and 49.4 per cent of the mass), for the pair of [P1] with $(s, delta) = (16, 0)$ and components $sqrt(2) D_16^+$, $sqrt(2) D_16^+$ (list identifier `16-16-s16d0c`; $74 slash "1,141"$ of the double cosets) and the pair with $(s, delta) = (16, 1)$ and both components $sqrt(2) (D_8^2)^+$ (`16-16d8`; $"44,686" slash "222,208"$).

= 9. Computations <s9>

#policy-box[
*What this section does.* We list the computations of this paper with their code. Computations of [P1] are cited as [P1, C$n$], with the numbers of the version of [P1] in the references.
]

#calc2-table(caption: [Table 11. Computations and code. References such as \[A1\] in the text point to the rows of this table.])

The programs and their versions are those of [P1, §7.1]: GAP 4.15.1 [GAP], bliss 0.77 [JK07], PARI/GP 2.17.3 [PARI] and SageMath 10.9 [Sage], with our own code in Python 3.12. The data of this paper are the following files of this repository; their SHA-256 digests are in `MANIFEST.json`. The programs and runs of Table 11 are in the research record of [P1] (§7.4 there).

#ptable(
  columns: 2, align: (left, left),
  caption: [Table 12. Data and checks of this paper.],
  table.header([file], [content]),
  [`data/largest-groups.jsonl`], [for each of the 99 lattices: its label, Gram matrix and generators of $"Aut" L$ (for six lattices a list that differs from that of [D] by one added or replaced matrix; both generate $"Aut" L$), its presentation (generators as integer matrices, relations, relative orders), its profile (§2.2: the structure, the conjugacy classes with their frame shapes and fixed lattices, the orbits on the minimal vectors, and the construction by an involution) and the columns of Table A1],
  [`paper/p2-lattices-99.json`], [the columns of Table A1, read by this paper],
  [`verify/verify_p2.py`], [checks in integers that the generators are automorphisms, that the relations of each presentation hold and that the generators of $"Aut" L$ of [P1] are products of the stored words, and that the bound of its lemma is $|"Aut" L|$; the orders obtained by coset enumeration are inputs (§2.5); standard Python only],
)

= Appendix A. The 99 lattices <appA>

#policy-box[
*What this appendix does.* One row per lattice, in the order of the labels. The columns are read from the profiles; the full profiles are in `data/largest-groups.jsonl` (Table 12).
]

#let chief-math(c) = {
  let toks = c.split(" ").filter(t => t != "").map(t => {
    let mark = if t.ends-with("*") { "*" } else if t.ends-with("?") { "?" } else { "" }
    let b = if mark != "" { t.slice(0, t.len() - 1) } else { t }
    let m = if mark == "" { b } else if b.contains("^") and not b.contains("_") { b.replace(regex("\\^(\\d+)$"), m => "^(" + m.captures.at(0) + " " + mark + ")") } else { b + "^" + mark }
    eval(m, mode: "math")
  })
  toks.join([, ])
}
#let fstar-math(segs) = segs.map(((kind, s)) => if kind == "math" { eval(s, mode: "math") } else { [#s] }).join([ ])

#{
  set text(size: 7pt)
  set par(justify: false)
  show figure: set block(breakable: true)
  ptable(
    columns: (auto, auto, 1.6fr, auto, 2fr, auto, auto), align: (left, left, left, left, left, center, left),
    caption: [Table A1. The 99 lattices. $F^*(G)$: for a $p$-group that is not extraspecial, its order with the orders of its centre and of its Frattini subgroup ($Phi$ in the table) and its nilpotency class; $compose E$ denotes the layer. Chief factors are listed from the top; $*$ marks a non-split factor and $?$ a factor whose splitting was not computed; for $p$-groups the chief series is not unique and the marks depend on the series chosen. "pres." is the lemma of §2.3 that proves the presentation; "orig." is the origin of extremality of §8 (comp.: the components alone; gluing: the gluing map).],
    table.header([label], [$|"Aut" L|$], [$F^*(G)$], [$|G slash F^*|$], [chief factors], [pres.], [orig.]),
    .._a1.map(r => (
      raw(r.label),
      [#eval(r.aut_factor.replace(" ", " dot "), mode: "math")],
      [#fstar-math(r.fstar)],
      [#_grouped(r.g_mod_fstar)],
      [#chief-math(r.chief)],
      [#r.presentation],
      [#r.origin],
    )).flatten(),
  )
}

= References <refs>
#refs-list(
  [[ATLAS] R. A. Wilson et al., _ATLAS of Finite Group Representations_, version 3, https://brauer.maths.qmul.ac.uk/Atlas/v3/, the pages of $S_6(2)$, $"U"_4 (2)$, $"U"_3 (3)$, $Omega_8^+(2)$ and $Omega_10^+(2)$; the orders and indices of the maximal subgroups were checked with the character table library of GAP.],
  [[BCDYY25] P. Boyvalenkov, D. Cherkashin, P. Dragnev, D. Yorgov, V. Yorgov, _On $T$-avoiding spherical codes and designs in $RR^32$_, arXiv:2512.23092.],
  [[BN97] C. Bachoc, G. Nebe, _Classification of two genera of 32-dimensional lattices of rank 8 over the Hurwitz order_, Experiment. Math. 6 (1997) 151–162.],
  [[D] The database of extremal even unimodular lattices of rank 32 with $|"Aut" L| > 2$, version 2: the release asset `database-v2.tar.gz` of this repository, version #raw(_version), together with `data/classified-lattices.jsonl` of the repository.],
  [[GAP] The GAP Group, _GAP – Groups, Algorithms, and Programming_, Version 4.15.1, 2025, https://www.gap-system.org.],
  [[Hen13] A. Henn, M. Hentschel, A. Krieg, G. Nebe, _On the classification of lattices over $QQ(sqrt(-3))$ which are even unimodular $ZZ$-lattices of rank 32_, Int. J. Math. Math. Sci. 2013, Art. 837080.],
  [[JK07] T. Junttila, P. Kaski, _Engineering an efficient canonical labeling tool for large and sparse graphs_, in: Proc. ALENEX 2007, SIAM, 135–149.],
  [[Jür15] M. Jürgens, _Nicht-Existenz und Konstruktion extremaler Gitter_, Dissertation, TU Dortmund, 2015.],
  [[KN93] H. Koch, G. Nebe, _Extremal even unimodular lattices of rank 32 and related codes_, Math. Nachr. 161 (1993) 309–319.],
  [[KV89] H. Koch, B. B. Venkov, _Über ganzzahlige unimodulare euklidische Gitter_, J. reine angew. Math. 398 (1989) 144–168.],
  [[Neb98q] G. Nebe, _Finite quaternionic matrix groups_, Represent. Theory 2 (1998) 106–223 (cited from the author's preprint; section and theorem numbers as there).],
  [[NS] G. Nebe, N. J. A. Sloane, _A Catalogue of Lattices_, https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/.],
  [[P1] roiban1344, _The extremal even unimodular lattices of rank 32: their number, and the classification of those with automorphisms other than $plus.minus 1$_, version #raw(_version), the file `paper/extremal-lattices-rank-32.pdf` of this repository.],
  [[PARI] The PARI Group, _PARI/GP_, version 2.17.3, Univ. Bordeaux, 2025, http://pari.math.u-bordeaux.fr/.],
  [[Sage] The Sage Developers, _SageMath, the Sage Mathematics Software System_, Version 10.9, 2026, https://www.sagemath.org.],
  [[SPLAG] J. H. Conway, N. J. A. Sloane, _Sphere Packings, Lattices and Groups_, 3rd ed., Springer, 1999.],
)
