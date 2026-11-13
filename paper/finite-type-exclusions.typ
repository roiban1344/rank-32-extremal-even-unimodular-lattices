#import "style.typ": *
#import "computations.typ": *

=== B.1.1 Exact exclusion of the six remaining types

The six types excluded by linear programming in Table B1 can be checked by finite rational certificates. Here we specify the lattice-to-linear-system step and the bases used in that calculation #calc("C3").

*One coset.* Let $X$ have minimum at least 4. Distinct vectors $v_1,dots,v_r$ of norm $alpha$, $0<alpha<2$, in the same coset of $X^\# slash X$ satisfy $(v_i,v_j)<=alpha-2$. Hence
$ 0<=norm(sum_i v_i)^2<=r alpha+r(r-1)(alpha-2)=r(2-r(2-alpha)). $
There are therefore at most $m(alpha)=floor(2 slash (2-alpha))=1+floor(alpha slash (2-alpha))$ such vectors. For $alpha<1$ the bound is 1. Vectors are counted with both signs.

*Gluing inequality.* Write $N(j)$ for the number of discriminant classes of $F$ of value $2j slash p$ modulo $2ZZ$, and $c_(F)(j),c_(Z)(j)$ for the coefficients of $q^j$ in the theta series of $sqrt(p)F^\#$ and $sqrt(p)Z^\#$. For $1<=j<p$,
#tagged("B.LP1", $ c_(F)(j) slash m(2j slash p)+c_(Z)(p-j) slash m(2-2j slash p)<=N(j). $)
Indeed the classes containing the first kind of vector and the inverse images of the classes containing the second kind are disjoint: otherwise their gluing would give a root. The coset bound applied to each set proves the inequality.

*The spaces and Fricke action.* For an even $p$-elementary lattice $X$ of rank $n=2k$ and determinant $p^s$, put $theta_(X)(tau)=sum_x e^(pi i tau(x,x))$. The lattice $sqrt(p)X^\#$ is even: $p(x,x)$ is integral and $p^2(x,x)$ is even. The theta series lies in $M_(k)(Gamma_(0)(p),chi)$ with $chi(a)=((-1)^k p^s slash a)$ [BV01, Thm 2.1]. Our Fricke convention is
$ (f|W_p)(tau)=p^(k slash 2)(p tau)^(-k) f(-1 slash (p tau)). $
Poisson summation gives $theta_X|W_p=i^(-k)p^((k-s) slash 2)theta_(sqrt(p)X^\#)$. The matrix in [BV01, Prop. 3.1] is the negative of this convention, accounting for its phase $i^k$. Thus a basis element coming from determinant $p^(s_i)$ has normalized dual component $p^((s-s_i) slash 2)theta_(sqrt(p)X_i^\#)$. The exponent is integral since $s_i equiv s$ modulo 2.

The following bases have the stated dimension. Their 64-coefficient matrices have full rank over $QQ$. Since every listed theta belongs to the specified space and the rank equals its dimension, they span that space; no converse inference from a truncated identity is used. Dimensions are computed with SageMath; the theta coefficients use PARI enumeration and the specified code. These are explicit computation assumptions.

#table(columns: (auto,auto,auto,auto,1fr),
  table.header([$p$], [rank], [$s$], [dimension], [basis of component theta series; determinant exponent]),
[3], [8], [6], [2], [$A_2 perp sqrt(3) E_6^\#$ ($s_i=6$); $A_2 perp E_6$ ($s_i=2$)],
[3], [10], [7], [2], [$(A_2)^2 perp sqrt(3) E_6^\#$ ($s_i=7$); $(A_2)^2 perp E_6$ ($s_i=3$)],
[3], [10], [9], [2], [$(A_2)^2 perp sqrt(3) E_6^\#$ ($s_i=7$); $(A_2)^2 perp E_6$ ($s_i=3$)],
[3], [14], [5], [3], [$A_2 perp (sqrt(3) E_6^\#)^2$ ($s_i=11$); $A_2 perp E_6 perp sqrt(3) E_6^\#$ ($s_i=7$); $A_2 perp (E_6)^2$ ($s_i=3$)],
[3], [16], [6], [3], [$(A_2)^2 perp (sqrt(3) E_6^\#)^2$ ($s_i=12$); $(A_2)^2 perp E_6 perp sqrt(3) E_6^\#$ ($s_i=8$); $(A_2)^2 perp (E_6)^2$ ($s_i=4$)],
[3], [18], [5], [4], [$(sqrt(3) E_6^\#)^3$ ($s_i=15$); $E_6 perp (sqrt(3) E_6^\#)^2$ ($s_i=11$); $(E_6)^2 perp sqrt(3) E_6^\#$ ($s_i=7$); $(E_6)^3$ ($s_i=3$)],
[3], [22], [7], [4], [$(A_2)^2 perp (sqrt(3) E_6^\#)^3$ ($s_i=17$); $(A_2)^2 perp E_6 perp (sqrt(3) E_6^\#)^2$ ($s_i=13$); $(A_2)^2 perp (E_6)^2 perp sqrt(3) E_6^\#$ ($s_i=9$); $(A_2)^2 perp (E_6)^3$ ($s_i=5$)],
[3], [22], [9], [4], [$(A_2)^2 perp (sqrt(3) E_6^\#)^3$ ($s_i=17$); $(A_2)^2 perp E_6 perp (sqrt(3) E_6^\#)^2$ ($s_i=13$); $(A_2)^2 perp (E_6)^2 perp sqrt(3) E_6^\#$ ($s_i=9$); $(A_2)^2 perp (E_6)^3$ ($s_i=5$)],
[3], [24], [6], [5], [$(sqrt(3) E_6^\#)^4$ ($s_i=20$); $E_6 perp (sqrt(3) E_6^\#)^3$ ($s_i=16$); $(E_6)^2 perp (sqrt(3) E_6^\#)^2$ ($s_i=12$); $(E_6)^3 perp sqrt(3) E_6^\#$ ($s_i=8$); $(E_6)^4$ ($s_i=4$)],
[7], [14], [3], [5], [$(Q_7)^7$ ($s_i=7$); $sqrt(7) A_6^\# perp (Q_7)^4$ ($s_i=9$); $(sqrt(7) A_6^\#)^2 perp Q_7$ ($s_i=11$); $A_6 perp (Q_7)^4$ ($s_i=5$); $A_6 perp sqrt(7) A_6^\# perp Q_7$ ($s_i=7$)],
[7], [18], [3], [7], [$(Q_7)^9$ ($s_i=9$); $sqrt(7) A_6^\# perp (Q_7)^6$ ($s_i=11$); $(sqrt(7) A_6^\#)^2 perp (Q_7)^3$ ($s_i=13$); $(sqrt(7) A_6^\#)^3$ ($s_i=15$); $A_6 perp (Q_7)^6$ ($s_i=7$); $A_6 perp sqrt(7) A_6^\# perp (Q_7)^3$ ($s_i=9$); $A_6 perp (sqrt(7) A_6^\#)^2$ ($s_i=11$)],
)

Here $sqrt(3)E_6^\#$ has Gram $3G(E_6)^(-1)$, $sqrt(7)A_6^\#$ has Gram $7G(A_6)^(-1)$, and $Q_7$ has Gram $mat(2,1;1,4)$. Powers denote orthogonal sums of identical factors. Positive definiteness, evenness, determinants, evenness of the $p$-scaled duals, and equality of the characters are checked in integers and rationals. The machine names of the two scaled duals and of $Q_7$ are `rt3E6d`, `rt7A6d`, and `Q7`, respectively.

*The discriminant sign.* At the odd prime $p$, write $F_p=U_0 perp p U_1$, $"rank" U_1=s$, and $epsilon=(det U_1 slash p)$. The other odd places have excess zero and the even unimodular 2-adic part has oddity zero. The oddity formula [SPLAG, Ch. 15, §7.7, (30), (35)] gives
$ d+s(p-1)+4[epsilon=-1] equiv 0 med (mod 8), quad epsilon=(-1)^((d+s(p-1)) slash 4). $
On $F^\# slash F$ define $Q(a)=j$ when $q(a)=2j slash p$ modulo $2ZZ$. Its matrix is $2^(-1)U_1^(-1)$ over $FF_p$, so its determinant character is $(2 slash p)^s epsilon$. For $p=3,7$ this equals $(-1)^((d-2s) slash 4)$. Counting values in a diagonal form with this determinant gives $N(j)$; only the zero vector is removed from $N(0)$.

*Finite system and certificates.* The variables are the rational coefficients in the displayed bases for $F$ and $Z$. Their theta and normalized dual series have constant term 1; the original series have zero $q^1$ coefficient and non-negative coefficients. The dual series are non-negative, and Lemma 4.2 makes the dual $Z$ coefficients at $q^1,dots,q^(p-2)$ zero. Add (B.LP1). The implementation initially imposes these conditions through degree 39. Each certificate below uses only the indicated equations and inequalities, of degree at most 6.

Each equation is represented by its two inequalities; “reverse” denotes the negative of that equation. `F:0:j` and `F:1:j` mean the equations for the coefficients of $q^j$ of the original and of the normalized dual series of $F$, respectively, and `F:dual:j>=0` the non-negativity of the coefficient of $q^j$ of the normalized dual series; similarly for $Z$. `gluing:j` means (B.LP1). Multiplying the indicated rows by the non-negative rational multipliers gives zero for every variable and right side $-1$.

#table(columns: (auto,auto,1fr),
  table.header([type $p - z - f - s$], [multiplier], [row]),
[3-8-16-6], [1/2], [F:0:0=1],
[3-8-16-6], [1/12], [F:0:1=0],
[3-8-16-6], [3/2], [F:1:0=1 (reverse)],
[3-8-16-6], [1/4], [F:dual:1>=0],
[3-9-14-5], [1/24], [F:0:1=0],
[3-9-14-5], [14/9], [F:1:0=1 (reverse)],
[3-9-14-5], [1/8], [F:dual:1>=0],
[3-9-14-5], [1/162], [Z:dual:1>=0],
[3-9-14-5], [1/162], [gluing:2],
[3-11-10-7], [1/4], [F:0:0=1 (reverse)],
[3-11-10-7], [1/12], [F:0:1=0],
[3-11-10-7], [3/4], [F:1:0=1 (reverse)],
[3-11-10-9], [1/2], [F:0:0=1 (reverse)],
[3-11-10-9], [1/6], [F:0:1=0],
[3-11-10-9], [1/2], [F:1:0=1 (reverse)],
[3-12-8-6], [1/6], [F:0:1=0],
[3-12-8-6], [40/27], [F:1:0=1 (reverse)],
[3-12-8-6], [1/486], [Z:dual:1>=0],
[3-12-8-6], [1/486], [gluing:2],
[7-3-14-3], [541/464], [F:0:0=1 (reverse)],
[7-3-14-3], [1/58], [F:0:1=0],
[7-3-14-3], [315/464], [F:1:0=1 (reverse)],
[7-3-14-3], [13/58], [F:dual:1>=0],
[7-3-14-3], [2/29], [F:dual:2>=0],
[7-3-14-3], [7/464], [Z:dual:1>=0],
[7-3-14-3], [7/464], [gluing:6],
)

Thus a necessary system for each of the six types implies $0<=-1$, excluding that type. The rational verification needs no linear-programming solver. The full matrices, coefficients, basis Grams and code hashes are fixed with the certificates #calc("C3"). This finite contradiction does not independently revalidate the software computations of the theta coefficients or space dimensions.
