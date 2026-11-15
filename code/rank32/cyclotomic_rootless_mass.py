r"""Rootless masses of the Hermitian genera (z, s) over Z[zeta_p] (p = 5, 7, 17, ...) without enumerating lattices.

The irreducible components of the root system: A = A_{p-1} (rank 1, |U| = 2p, any of its roots generates it)
and, for p = 5 only, E = E8 (rank 2, |U| = 600). With c(A) = -1 and c(E) = 5,
    m(empty) = sum_{a, b} c(A)^a c(E)^b A(A^a E^b) / ((2p)^a a! 600^b b!),
    A(T) = sum_{C isotropic} A_prim(T_C),
    A_prim(T') = sum_{H of each type} emb(H, D(T')) emb(H^-, D(K)) / |Aut H| * mass_H(z - n, s - s_{T'} + 2 dim H).
mass_H is the mass of the Hermitian genus from rank32.genus_mass (0 when the genus does not exist); the
standard mass agrees with [Jür15, Satz 3.3.15 and (45)] (the mode `std` of
code/produce/element_types/rootless_cyclotomic_mass_oddp.py).
  rootless(p, z, s)   m(empty), the mass of the rootless lattices of the genus
  table(p, z, s)      m(A^a E^b) for every root system A^a E^b (their sum is mass_H(p, z, s))
  order_O, sing, frames, gl, emb_nondeg, disc_V   counts in finite quadratic spaces over F_p
Used by computation C9 of the first paper.

Z[zeta_p] 上のエルミートの種 (z, s)（p = 5, 7, 17, ...）のルート無しの質量を、格子を列挙せずに求める。

ルート系の既約成分は A = A_{p-1}（階数 1、|U| = 2p、どのルートも全体を生成する）と、p = 5 のときだけの
E = E8（階数 2、|U| = 600）である。c(A) = -1、c(E) = 5 として
    m(空) = sum_{a, b} c(A)^a c(E)^b A(A^a E^b) / ((2p)^a a! 600^b b!)
    A(T) = sum_{C 等方} A_prim(T_C)
    A_prim(T') = sum_{H の型} emb(H, D(T')) emb(H^-, D(K)) / |Aut H| * mass_H(z - n, s - s_{T'} + 2 dim H)
mass_H は rank32.genus_mass のエルミートの種の質量（種が無ければ 0）。標準質量は
[Jür15, Satz 3.3.15 と式 (45)] と一致する（code/produce/element_types/rootless_cyclotomic_mass_oddp.py の
モード `std`）。
  rootless(p, z, s)   m(空)。種のルート無しの格子の質量
  table(p, z, s)      全てのルート系 A^a E^b についての m(A^a E^b)（和は mass_H(p, z, s)）
  order_O, sing, frames, gl, emb_nondeg, disc_V   F_p 上の有限 2 次空間の数え上げ
第 1 の論文の計算 C9 が使う。
"""
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial

from . import genus_mass as ptmass


def F(x):
    """A rational number of Sage as a Fraction.

    Sage の有理数を Fraction に直す。"""
    return Fr(int(x.numerator()), int(x.denominator()))


def leg(a, p):
    """The Legendre symbol (a / p) for an odd prime p.

    奇素数 p についての Legendre 記号 (a / p)。"""
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


# finite quadratic spaces (dimension n, discriminant class d = +1 / -1)
# 有限の 2 次空間（次元 n、判別式の類 d = +1 / -1）
def order_O(p, n, d):
    """The order of the orthogonal group of the quadratic space over F_p of dimension n and discriminant class d.

    F_p 上の次元 n・判別式の類 d の 2 次空間の直交群の位数。"""
    if n == 0:
        return 1 if d == 1 else 0
    if n % 2:
        k = (n - 1) // 2
        out = 2 * p ** (k * k)
        for i in range(1, k + 1):
            out *= p ** (2 * i) - 1
        return out
    k = n // 2
    eps = leg((-1) ** k, p) * d
    out = 2 * p ** (k * (k - 1)) * (p ** k - eps)
    for i in range(1, k):
        out *= p ** (2 * i) - 1
    return out


def sing(p, n, d):
    """The number of nonzero singular vectors.

    0 でない特異ベクトルの数。"""
    if n <= 0:
        return 0
    if n % 2:
        return p ** (n - 1) - 1
    k = n // 2
    eps = leg((-1) ** k, p) * d
    return (p ** k - eps) * (p ** (k - 1) + eps)


def frames(p, n, d, rho):
    """The number of ordered bases of totally isotropic subspaces of dimension rho, counted step by step.

    rho 次元の全等方部分空間の順序つきの基底の数（1 本ずつ数える）。"""
    out = 1
    for i in range(rho):
        out *= p ** i * sing(p, n - 2 * i, d * leg(-1, p) ** i)
    return out


def gl(p, k):
    """|GL_k(F_p)|.

    |GL_k(F_p)|。"""
    out = 1
    for i in range(k):
        out *= p ** k - p ** i
    return out


def emb_nondeg(p, nu, dN, w, dW):
    """The number of isometric injections of a non-degenerate N (dimension nu, discriminant dN) into W.

    非退化な N（次元 nu、判別式 dN）の W への等長な単射の数。"""
    if nu > w:
        return 0
    if nu == w and dN != dW:
        return 0
    oW = order_O(p, w, dW)
    oC = order_O(p, w - nu, dW * dN)
    if oC == 0 or oW == 0:
        return 0
    return Fr(oW, oC)


def disc_V(p, z, s):
    """The discriminant class of the discriminant form of the genus (z, s): an isotropic quotient of <-1/2>^z.

    種 (z, s) の判別式の形の判別式の類。<-1/2>^z の等方商。"""
    return (leg(-2, p) ** z) * (leg(-1, p) ** ((z - s) // 2))


@lru_cache(maxsize=None)
def mass_H(p, z, s):
    """The mass of the Hermitian genus (z, s) as a Fraction; 0 if the genus does not exist.

    エルミートの種 (z, s) の質量（Fraction）。種が無ければ 0。"""
    if z == 0:
        return Fr(1) if s == 0 else Fr(0)
    if s > z or (z - s) % 2:
        return Fr(0)
    n = z * (p - 1)
    if not ptmass.genus_exists(p, n, s):
        return Fr(0)
    return F(ptmass.mass_H(p, z, s))


def glue_sum(p, s_T, d_T, z_K, s_target):
    """sum_H emb(H, V) emb(H^-, D(K)) / |Aut H| mass_H(z_K, s_K).

    sum_H emb(H, V) emb(H^-, D(K)) / |Aut H| mass_H(z_K, s_K)。"""
    tot = Fr(0)
    for h in range(s_T + 1):
        s_K = s_target - s_T + 2 * h
        if s_K < 0:
            continue
        m = mass_H(p, z_K, s_K)
        if m == 0:
            continue
        dK = disc_V(p, z_K, s_K)
        for rho in range(h + 1):
            nu = h - rho
            for dN in ((1, -1) if nu > 0 else (1,)):
                e1 = emb_nondeg(p, nu, dN, s_T, d_T)
                if e1 == 0:
                    continue
                f1 = frames(p, s_T - nu, d_T * dN, rho)
                dNm = dN * leg(-1, p) ** nu
                e2 = emb_nondeg(p, nu, dNm, s_K, dK)
                if e2 == 0:
                    continue
                f2 = frames(p, s_K - nu, dK * dNm, rho)
                autH = gl(p, rho) * p ** (rho * nu) * (order_O(p, nu, dN) if nu > 0 else 1)
                tot += e1 * f1 * e2 * f2 / autH * m
    return tot


def A_value(p, z, s, a, b):
    """A(A^a E^b) in the genus (z, s).

    種 (z, s) での A(A^a E^b)。"""
    n = a + 2 * b
    if n > z:
        return Fr(0)
    tot = Fr(0)
    d0 = leg(-2, p) ** a
    for k in range(a // 2 + 1):
        cnt = Fr(frames(p, a, d0, k), gl(p, k))
        if cnt == 0:
            continue
        s_T = a - 2 * k
        d_T = d0 * leg(-1, p) ** k
        tot += cnt * glue_sum(p, s_T, d_T, z - n, s)
    return tot


def types(p, z):
    """The root systems (a, b) = A^a E^b that fit in rank z (b = 0 unless p = 5).

    階数 z に収まるルート系 (a, b) = A^a E^b（p = 5 以外は b = 0）。"""
    bmax = z // 2 if p == 5 else 0
    for b in range(bmax + 1):
        for a in range(z - 2 * b + 1):
            yield a, b


def sub_count(p, a, b, a2, b2):
    """The number of root sublattices A^{a2} E^{b2} of A^a E^b.

    A^a E^b の中の、部分ルート格子 A^{a2} E^{b2} の個数。"""
    if b2 > b:
        return 0
    from math import comb

    tot = 0
    # b2 of the E stay E; each of the other b - b2 becomes 0, A or A^2
    # E のうち b2 個は E のまま。残り b - b2 個の E は 0, A, A^2 のどれかに落ちる
    rest = b - b2
    for i in range(rest + 1):  # i of them become A (12 ways) / i 個が A（12 通り）
        for j in range(rest - i + 1):  # j of them become A^2 (6 ways) / j 個が A^2（6 通り）
            need = a2 - i - 2 * j
            if need < 0 or need > a:
                continue
            tot += comb(b, b2) * comb(rest, i) * comb(rest - i, j) * 12 ** i * 6 ** j * comb(a, need)
    return tot


def autT(p, a, b):
    """|U(A^a E^b)| = (2p)^a a! 600^b b!.

    |U(A^a E^b)| = (2p)^a a! 600^b b!。"""
    return (2 * p) ** a * factorial(a) * 600 ** b * factorial(b)


def table(p, z, s):
    """All m(A^a E^b), as a dict keyed by (a, b).

    全ての m(A^a E^b) を、(a, b) をキーとする dict で返す。"""
    ts = sorted(types(p, z), key=lambda t: (-(t[0] + 2 * t[1]), -t[1]))
    res = {t: A_value(p, z, s, *t) for t in ts}
    M = {}
    for R in ts:
        val = res[R] / autT(p, *R)
        M[R] = val
        if val == 0:
            continue
        for T in ts:
            if T != R:
                c = sub_count(p, R[0], R[1], T[0], T[1])
                if c:
                    res[T] -= val * c * autT(p, *T)
    return M


def rootless(p, z, s):
    """m(empty): the mass of the rootless lattices of the Hermitian genus (z, s).

    m(空)：エルミートの種 (z, s) のルート無しの格子の質量。"""
    tot = Fr(0)
    for a, b in types(p, z):
        c = Fr((-1) ** a * 5 ** b)
        tot += c * A_value(p, z, s, a, b) / autT(p, a, b)
    return tot
