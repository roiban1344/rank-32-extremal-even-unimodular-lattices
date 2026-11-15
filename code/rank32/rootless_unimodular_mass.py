"""Masses of the rootless Hermitian unimodular lattices over Z[zeta_p].

A root x (a vector of norm 2 of the trace form) of a Hermitian unimodular lattice U of rank z over
Z[zeta_p] spans Z[zeta_p] x = A_{p-1} (unimodular of rank 1), which splits off as an orthogonal summand.
The decomposition U = A_{p-1}^k + U' (U' rootless) is unique and the unitary group has order
(2p)^k k! |U(U')|, so
    mass_U(z) = sum_{k=0}^{z} mass_noroot(z - k) / ((2p)^k k!),
which this module inverts. mass_noroot(z) = 0 means that there is no rootless U of rank z.
Used by computation C2 of the first paper.

Z[zeta_p] 上のルートを持たないエルミート・ユニモジュラー格子の質量。

Z[zeta_p] 上の階数 z のエルミート・ユニモジュラー格子 U のルート x（トレース形式のノルム 2 の
ベクトル）は Z[zeta_p] x = A_{p-1}（階数 1 のユニモジュラー）を張り、直和因子として分かれる。
分解 U = A_{p-1}^k + U'（U' はルートを持たない）は一意で、ユニタリ群の位数は (2p)^k k! |U(U')| である。
よって
    mass_U(z) = sum_{k=0}^{z} mass_noroot(z - k) / ((2p)^k k!)
が成り立つ。このモジュールはこれを逆に解く。mass_noroot(z) = 0 なら、階数 z のルートを持たない U は無い。
第 1 の論文の計算 C2 が使う。
"""
from sage.all import ZZ, factorial

from . import genus_mass as ptmass


def noroot_masses(p, zmax):
    """The list [mass_noroot(0), ..., mass_noroot(zmax)] for the prime p.

    素数 p についての列 [mass_noroot(0), ..., mass_noroot(zmax)]。"""
    out = []
    for z in range(zmax + 1):
        m = ptmass.mass_U(p, z)
        for k in range(1, z + 1):
            m -= out[z - k] / (ZZ(2 * p) ** k * factorial(k))
        out.append(m)
    return out
