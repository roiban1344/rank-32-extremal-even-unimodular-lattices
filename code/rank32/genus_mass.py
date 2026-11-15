"""Genera, masses and the mass of pairs (L, sigma) for the types p-(z, d)-s of automorphisms of odd prime order.

Types as in [Neb13b, Def. 4.1]: F = Fix(sigma) (rank d) and the cyclic part Z (rank z over Z[zeta_p]) are
both even, p-elementary and of determinant p^s. The central formula is
    M(p; z, d, s) := sum_{(L, sigma) up to isomorphism} 1/|C_{Aut L}(sigma)| = mass(genus of F) * mass_H(genus of Z) * |O(q_s)|,
where L runs over all even unimodular lattices of rank n (not only the extremal ones) and sigma over the
elements of type (z, d, s);
    mass_H(Z) = mass_U(z) * I(z, m) / prod_{i=1}^{m} (p^i + 1),   m = (z - s)/2,
    mass_U(z) = 2^(-z r) * prod_{j=1}^{z} |L_{K+}(1 - j, chi^j)| * lambda_z,   r = (p - 1)/2,
    lambda_z = 1 (z odd),  p^(z/2) + ((-1)^(z/2) / p) (z even).
mass_U is the formula of Hashimoto-Koseki (cited in the proof of [HKN09, Theorem 1]) written for a general p.
For p = 3 and z = 4, 8, 12, 16 it is compared with the values of [HKN09]; the whole formula is compared
with the actual numbers of elements for E8 and in rank 16 (the mode `test` of
code/produce/element_types/prime_order_type_masses.py).
Used by computation C2 of the first paper and by the programs that build on it.

奇素数位数の自己同型の型 p-(z, d)-s についての、種・質量・対 (L, sigma) の質量。

型は [Neb13b, Def. 4.1] のとおり。F = Fix(sigma)（階数 d）と巡回部分 Z（Z[zeta_p] 上の階数 z）は、
どちらも偶・p-基本・行列式 p^s である。中心の式は
    M(p; z, d, s) := sum_{(L, sigma)/同型} 1/|C_{Aut L}(sigma)| = mass(F の種) * mass_H(Z の種) * |O(q_s)|
である。左辺の L は階数 n の偶ユニモジュラー格子の全部（極値に限らない）、sigma は型 (z, d, s) の元を動く。
    mass_H(Z) = mass_U(z) * I(z, m) / prod_{i=1}^{m} (p^i + 1),   m = (z - s)/2
    mass_U(z) = 2^(-z r) * prod_{j=1}^{z} |L_{K+}(1 - j, chi^j)| * lambda_z,   r = (p - 1)/2
    lambda_z = 1（z が奇数）、p^(z/2) + ((-1)^(z/2) / p)（z が偶数）
mass_U は Hashimoto-Koseki の式（[HKN09, Theorem 1] の証明が引く）を、一般の p について書いたもの。
p = 3、z = 4, 8, 12, 16 では [HKN09] の値と照合する。全体の式は、E8 と階数 16 の実際の元の数と照合する
（code/produce/element_types/prime_order_type_masses.py のモード `test`）。
第 1 の論文の計算 C2 と、その上に立つプログラムが使う。
"""
from sage.all import ZZ, QQ, bernoulli, DirichletGroup, CyclotomicField, prod, kronecker_symbol
from sage.quadratic_forms.genera.genus import genera

from . import density_sieve as sieve


def genus_exists(p, n, s):
    """Whether the genus of positive definite, even, p-elementary lattices of rank n and determinant p^s exists.

    For an odd prime p there is at most one such genus (SPLAG, Chapter 15, Theorem 13). The oddity
    condition of the local symbol 1^{eps n}_{II} 1^{eps1 (n - s)} p^{eps2 s}: n + s (p - 1) + 4k = 0 mod 8,
    where k = 1 if and only if eps2 = -1.

    偶・p-基本・階数 n・行列式 p^s の正定値格子の種が存在するか。
    奇素数 p なら高々 1 個（SPLAG 第 15 章 Theorem 13）。局所記号 1^{eps n}_{II} 1^{eps1 (n - s)} p^{eps2 s}
    の奇数条件（oddity formula）：n + s (p - 1) + 4k = 0 mod 8、k = 1 は eps2 = -1 と同値。"""
    if n % 2:
        return False
    if s < 0 or s > n:
        return False
    if s == 0 or s == n:
        return n % 8 == 0  # one Jordan component only: the sign is + (the determinant is the square class 1) / 片方の Jordan 成分だけ：符号は + に決まる（行列式が平方類 1）
    return (n + s * (p - 1)) % 4 == 0


def eps2(p, n, s):
    """The sign of the p-adic Jordan component p^{eps2 s}.

    p 進の Jordan 成分 p^{eps2 s} の符号。"""
    return 1 if (n + s * (p - 1)) % 8 == 0 else -1


def q_type(p, n, s):
    """The type of the discriminant quadratic form (on F_p^s): +1 split, -1 non-split for s even; 0 for s odd.

    判別式の 2 次形式（F_p^s 上）の型。s が偶数なら分裂で +1、非分裂で -1。s が奇数なら 0。"""
    if s % 2:
        return 0
    m = s // 2
    return kronecker_symbol(-1, p) ** m * eps2(p, n, s)


def order_O(q, s, e):
    """The order of the orthogonal group of a non-degenerate quadratic space of dimension s over F_q (e: its type for s even).

    F_q 上の s 次元の非退化 2 次空間の直交群の位数。e は型（s が偶数のとき +1 か -1）。"""
    if s == 0:
        return ZZ(1)
    if s % 2:
        m = (s - 1) // 2
        return 2 * ZZ(q) ** (m * m) * prod(ZZ(q) ** (2 * i) - 1 for i in range(1, m + 1))
    m = s // 2
    return 2 * ZZ(q) ** (m * (m - 1)) * (ZZ(q) ** m - e) * prod(ZZ(q) ** (2 * i) - 1 for i in range(1, m))


_gen_cache = {}


def genus_F(p, n, s):
    """The genus (a Genus of Sage) of even p-elementary lattices of rank n and determinant p^s, or None.

    偶・p-基本・階数 n・行列式 p^s の種（Sage の Genus）。無ければ None。"""
    if n == 0:
        return None
    key = (p, n, s)
    if key in _gen_cache:
        return _gen_cache[key]
    g = None
    for G in genera((n, 0), ZZ(p) ** s, even=True):
        sym = G.local_symbol(p).symbol_tuple_list()
        scales = {t[0]: t[1] for t in sym}
        if all(k in (0, 1) for k in scales) and scales.get(1, 0) == s:
            g = G
            break
    _gen_cache[key] = g
    return g


def mass_Z_genus(p, n, s):
    """The mass of the genus of even p-elementary Z-lattices (Sage, the mass formula of Conway-Sloane).

    偶・p-基本の Z-格子の種の質量（Sage、Conway-Sloane の質量公式）。"""
    if n == 0:
        return QQ(1)
    g = genus_F(p, n, s)
    return None if g is None else g.mass()


_char_cache = {}


def _odd_even_chars(p):
    if p not in _char_cache:
        D = DirichletGroup(p, CyclotomicField(p - 1))
        odd = [c for c in D if c(-1) == -1]
        even = [c for c in D if c(-1) == 1 and not c.is_trivial()]
        _char_cache[p] = (odd, even)
    return _char_cache[p]


def l_value(p, j):
    """|L_{K+}(1 - j, chi^j)|.

    For j odd: prod_{chi odd} L(1 - j, chi). For j even: zeta_{K+}(1 - j) = zeta(1 - j) prod_{chi even, chi != 1} L(1 - j, chi).
    L(1 - j, chi) = -B_{j, chi}/j.

    |L_{K+}(1 - j, chi^j)|。
    j が奇数なら prod_{chi 奇} L(1 - j, chi)。j が偶数なら zeta_{K+}(1 - j) = zeta(1 - j) prod_{chi 偶, chi != 1} L(1 - j, chi)。
    L(1 - j, chi) = -B_{j, chi}/j。"""
    odd, even = _odd_even_chars(p)
    if j % 2:
        val = prod(-c.bernoulli(j) / j for c in odd)
    else:
        val = (-bernoulli(j) / j) * prod((-c.bernoulli(j) / j for c in even), z=QQ(1))
    return abs(QQ(val))


def lambda_U(p, z):
    """lambda_z: 1 for z odd, p^(z/2) + ((-1)^(z/2) / p) for z even.

    lambda_z：z が奇数なら 1、偶数なら p^(z/2) + ((-1)^(z/2) / p)。"""
    if z % 2:
        return ZZ(1)
    return ZZ(p) ** (z // 2) + kronecker_symbol(-1, p) ** (z // 2)


def mass_U(p, z):
    """The mass of the genus of Hermitian unimodular lattices of rank z over Z[zeta_p] (the dual for the Hermitian form is the lattice itself).

    Z[zeta_p] 上の階数 z のエルミート・ユニモジュラー格子（エルミート形式についての双対が自分）の種の質量。"""
    if z == 0:
        return QQ(1)
    r = (p - 1) // 2
    return QQ(1) / 2 ** (z * r) * prod(l_value(p, j) for j in range(1, z + 1)) * lambda_U(p, z)


def isotropic_count(q, z, m):
    """The number of totally isotropic subspaces of dimension m of U/pi U (dimension z over F_q, discriminant the square class 1).

    U/pi U（F_q 上 z 次元、判別式は平方類 1）の、m 次元の全等方部分空間の数。"""
    if m == 0:
        return ZZ(1)
    if z % 2:
        # odd dimension, by the general formula / 奇数次元：一般の公式で
        a = (z - 1) // 2
        # number of totally isotropic subspaces of dimension m = prod_{i=0}^{m-1} (q^(2(a-i)) - 1)/(q^(i+1) - 1)
        # 全等方 m 次元の数 = prod_{i=0}^{m-1} (q^(2(a-i)) - 1)/(q^(i+1) - 1)
        num = prod(ZZ(q) ** (2 * (a - i)) - 1 for i in range(m))
    else:
        a = z // 2
        e = kronecker_symbol(-1, q) ** a  # +1 if split / 分裂なら +1
        # prod_{i=0}^{m-1} (q^(a-i) - e)(q^(a-i-1) + e)/(q^(i+1) - 1)
        num = prod((ZZ(q) ** (a - i) - e) * (ZZ(q) ** (a - i - 1) + e) for i in range(m))
    den = prod(ZZ(q) ** (i + 1) - 1 for i in range(m))
    r = QQ(num) / den
    assert r.denominator() == 1, (q, z, m, r)
    return ZZ(r)


def mass_H(p, z, s):
    """The mass of the Hermitian genus (z, s) (measured with the orders of the unitary groups); None if there is no such genus.

    エルミートの種 (z, s) の質量（ユニタリ群の位数で測る）。種が無ければ None。"""
    if (z - s) % 2 or s < 0 or s > z:
        return None
    m = (z - s) // 2
    return mass_U(p, z) * isotropic_count(p, z, m) / prod((ZZ(p) ** i + 1 for i in range(1, m + 1)), z=ZZ(1))


def M_pair(p, n, z, s):
    """sum_{(L, sigma)} 1/|C(sigma)| over all even unimodular lattices of rank n; 0 if a genus is missing.

    階数 n の偶ユニモジュラー格子の全体での sum_{(L, sigma)} 1/|C(sigma)|。種が無ければ 0。"""
    d = n - z * (p - 1)
    if d < 0 or s > min(z, d) or (z - s) % 2:
        return QQ(0)
    if not genus_exists(p, z * (p - 1), s) or (d > 0 and not genus_exists(p, d, s)):
        return QQ(0)
    if d == 0 and s != 0:
        return QQ(0)
    mF = mass_Z_genus(p, d, s)
    mH = mass_H(p, z, s)
    if mF is None or mH is None:
        return QQ(0)
    return mF * mH * order_O(p, s, q_type(p, d, s))
