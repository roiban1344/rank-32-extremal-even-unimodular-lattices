#import "style.typ": *
#import "computations.typ": *

*Closed sums at odd primes.* Put $R_p=ZZ_p[i]$ and $t=r-n$. At odd primes the root inputs of Gaussian rank at most 16 have
$ D_(p)(T) tilde.equiv (R_p slash p)^a perp (R_p slash p^2)^b, $
where $b=0,1$, or $(p,a,b)=(3,0,2)$. Only $A_8[i]$ contributes a $p^2$ factor. For a residual module requiring $d$ generators, abbreviate
$ ell_(p)(d)=product_(j=t-d+1)^t (1-chi(p)^j p^(-j)). $
Every submodule contributes $|C|^(-t)ell_(p)(d(C))$ to (C.6); $d>t$ gives the zero factor at $j=0$.

If $p equiv 1 med (mod 4)$, then $b=0$ and the discriminant splits into two mutually dual $a$-dimensional spaces over $FF_p$. An isotropic submodule is a pair $(C_+,C_-)$ with $C_- subset C_+^perp$. For dimensions $k_1,k_2$ its multiplicity is
$ binom(a,k_1)_p binom(a-k_1,k_2)_p, $
where $binom(n,k)_p$ counts $k$-subspaces over $FF_p$. It has $|C|=p^(k_1+k_2)$ and $d(C)=a-k_1-k_2$.

If $p equiv 3 med (mod 4)$, the residue field is $FF_(p^2)$. A non-degenerate Hermitian $a$-space has
$ N(a,k)=(product_(j=0)^(k-1) (p^(a-2j)-(-1)^(a-2j))(p^(a-2j-1)-(-1)^(a-2j-1))) slash (product_(j=1)^k (p^(2j)-1)) $
totally isotropic $k$-subspaces. This is zero when $2k>a$ and one when $k=0$. At step $j$, each chosen nonzero isotropic vector of the orthogonal quotient has $p^(2j)$ lifts. These factors cancel the same powers of $p$ in the number of ordered bases of a $k$-space over $FF_(p^2)$, leaving the displayed denominator. For $b=0$ these subspaces have $|C|=p^(2k)$ and $d(C)=a-2k$.

For $b=1$, apply the graph transition to the elementary summand and the cyclic summand $B=R_p slash p^2$. After choosing $C_0$ of dimension $k$ in the elementary summand, its quotient $A_0$ has dimension $a'=a-2k$. The projection $Z$ to $B$ is $0$, $p B$ or $B$. The last case is impossible: the image of a generator in $A_0$ has norm denominator at most $p$, which cannot cancel the unit denominator $p^2$ on $B$. The remaining choices are

#table(columns: 4,
  table.header([projection and map], [multiplicity], [$log_p|C|$], [$d(C)$]),
  [$Z=0$], [1], [$2k$], [$a'+1$],
  [$Z=p B$, zero map], [1], [$2k+2$], [$a'$],
  [$Z=p B$, nonzero isotropic image], [$N(a',1)(p^2-1)$], [$2k+2$], [$a'-1$],
)

The final quotient in the last row replaces one hyperbolic plane of $A_0$ and $B$ by a cyclic $p^2$ factor; the other $a'-2$ elementary factors remain. Multiplying each row by $N(a,k)$ and its local weight gives the three terms used in the code.

For $(p,a,b)=(3,0,2)$, $D_3 (A_8 [i]^2) = (R_3 slash 9)^2$ with the form $"diag"(8 slash 9, 8 slash 9)$; since 8 is a unit at 3, the implementation uses $"diag"(1 slash 9, 1 slash 9)$, which has the same isotropic submodules and the same orthogonals $C^perp$, hence the same $|C|$ and $d(C)$, which are all that (C.6) uses. It works with the actual finite set $(R_3 slash 9)^2$. Starting with zero, it adds $x,i x$ for isotropic $x in C^perp$ and takes their additive span, discarding duplicate subsets. Every isotropic submodule is reached by adjoining its generators in this way. The subgroup of $C^perp slash C$ killed by 3 has size $9^(d(C))$, which determines the required generator number. This finite enumeration, the closed odd-prime sums and the 2-adic transition below are the local operations of #calc("C11").

