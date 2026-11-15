"""Density bounds that exclude types p-(z, d)-s of automorphisms of prime order.

An automorphism sigma of odd prime order p of an even unimodular lattice of rank n has a type
p-(z, d)-s [Neb13b, Def. 4.1]: d = dim Fix(sigma), z is the rank of the cyclic part over Z[zeta_p]
(n = d + z (p - 1)), and det Fix = det Cyc = p^s. Known conditions: s <= min(z, d), and z - s is even
[Neb13b, Lemma 4.3]. Fix(sigma) (rank d) and Cyc(sigma) (rank z (p - 1)) both have minimum at least
that of the lattice and determinant p^s, so an upper bound on the density gives a necessary condition.
  ok(n, p, s, mn)          the rank-n lattices of minimum >= mn and determinant p^s do not violate
                           the bound: for n <= 36 the centre density of [CE03, Table 3] ("New upper
                           bound", copied into CE below), above that Blichfeldt's bound
                           gamma_n <= (2/pi) Gamma(2 + n/2)^(2/n); n = 0 is never excluded.
  remark32(p, z, d, s, mn) the conditions of [Neb16, Remark 3.2 (b), (d)].
Used by computation C2 of the first paper (the table of types). The callers of ok() in rank 32
(rank32.joint_theta_lp.types_for and the function table of produce/element_types/prime_order_type_masses.py)
pass mn = 4, the minimum of an extremal lattice of rank 32; the default mn = 10 is not used there.
[Neb16] G. Nebe, Automorphisms of extremal unimodular lattices in dimension 72, J. Number Theory 161 (2016);
arXiv:1409.8473.

素数位数の自己同型の型 p-(z, d)-s を、密度の上界で除く。

階数 n の偶ユニモジュラー格子の奇素数位数 p の自己同型 sigma には型 p-(z, d)-s がある
[Neb13b, Def. 4.1]：d = dim Fix(sigma)、z は巡回部分の Z[zeta_p] 上の階数（n = d + z (p - 1)）、
det Fix = det Cyc = p^s。既知の条件は s <= min(z, d) と、z - s が偶数であること [Neb13b, Lemma 4.3]。
Fix(sigma)（階数 d）と Cyc(sigma)（階数 z (p - 1)）は、どちらも最小が格子の最小以上で、行列式が
p^s である。そこで密度の上界から必要条件が出る。
  ok(n, p, s, mn)          最小 >= mn・行列式 p^s の階数 n の格子が上界に反しない。n <= 36 では
                           [CE03, Table 3] の中心密度（「New upper bound」の列。下の CE に転記）、
                           それより上は Blichfeldt の上界 gamma_n <= (2/pi) Gamma(2 + n/2)^(2/n)。
                           n = 0 は除かない。
  remark32(p, z, d, s, mn) [Neb16, Remark 3.2 (b), (d)] の条件。
第 1 の論文の計算 C2（型の表）が使う。階数 32 で ok() を呼ぶ側（rank32.joint_theta_lp.types_for と、
produce/element_types/prime_order_type_masses.py の関数 table）は、階数 32 の極値格子の最小 mn = 4 を渡す。
既定値 mn = 10 はそこでは使わない。
[Neb16] G. Nebe, Automorphisms of extremal unimodular lattices in dimension 72, J. Number Theory 161 (2016);
arXiv:1409.8473.
"""
import math

EXACT = {1: 1.0, 2: (4 / 3) ** (1 / 2), 3: 2 ** (1 / 3), 4: 4 ** (1 / 4), 5: 8 ** (1 / 5), 6: (64 / 3) ** (1 / 6), 7: 64 ** (1 / 7), 8: 2.0, 24: 4.0}


def gamma_upper(n):
    """An upper bound on the Hermite constant gamma_n (exact for the known n, Blichfeldt otherwise).

    エルミート定数 gamma_n の上界（値の知られた n では厳密な値、それ以外は Blichfeldt）。"""
    if n in EXACT:
        return EXACT[n]
    return (2 / math.pi) * math.gamma(2 + n / 2) ** (2 / n)


# Cohn-Elkies, New upper bounds on sphere packings I, Table 3, column "New upper bound" (centre density, n <= 36)
# Cohn-Elkies, New upper bounds on sphere packings I, Table 3 の「New upper bound」の列（中心密度、n <= 36）
CE = [None, 0.5, 0.28868, 0.18616, 0.13126, 0.09975, 0.08084, 0.06933, 0.06251, 0.05900, 0.05804, 0.05932, 0.06279,
      0.06870, 0.07750, 0.08999, 0.10738, 0.13150, 0.16503, 0.21202, 0.27855, 0.37389, 0.51231, 0.71601, 1.01998,
      1.48001, 2.18614, 3.28537, 5.02059, 7.79782, 12.30390, 19.71397, 32.06222, 52.90924, 88.55925, 150.29783, 258.54994]


def ok(n, p, s, mn=10):
    """True unless a rank-n lattice of minimum >= mn and determinant p^s would violate the density bound.

    最小 >= mn・行列式 p^s の階数 n の格子が密度の上界に反しないなら True。"""
    if n == 0:
        return True
    if n <= 36:
        # centre density (mn/4)^(n/2) / sqrt(p^s) <= CE[n] (the table has 5 digits, so the bound is relaxed slightly)
        # 中心密度 (mn/4)^(n/2) / sqrt(p^s) <= CE[n]（表は 5 桁なので僅かに緩める）
        return (n / 2) * math.log(mn / 4) - (s / 2) * math.log(p) <= math.log(CE[n] * (1 + 1e-4))
    return mn <= gamma_upper(n) * (p ** s) ** (1 / n) + 1e-12


def remark32(p, z, d, s, mn=10):
    """The conditions of [Neb16, Remark 3.2 (b), (d)].

    If s = 0, then z is even. If s = d (p odd), then F/sqrt(p) is even unimodular, so d is a multiple
    of 8 and its minimum is at least mn/p; the minimum of an even unimodular lattice is at most
    2 floor(d/24) + 2 (Rains-Sloane), so only the d for which this fits are kept.

    [Neb16, Remark 3.2 (b), (d)] の条件。
    s = 0 なら z は偶数。s = d（p は奇数）なら F/sqrt(p) は偶ユニモジュラーなので d は 8 の倍数で、
    その最小は mn/p 以上。偶ユニモジュラー格子の最小は 2 floor(d/24) + 2 以下（Rains-Sloane）なので、
    これに収まる d だけを残す。"""
    if s == 0 and z % 2:
        return False
    if d > 0 and s == d:
        if d % 8:
            return False
        need = math.ceil(mn / p)
        need += need % 2  # the minimum of an even lattice is even / 偶格子の最小は偶数
        if 2 * (d // 24) + 2 < need:
            return False
    return True
