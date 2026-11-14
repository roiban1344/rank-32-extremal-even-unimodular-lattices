"""Exact checks that an integer matrix is an isometry between two lattices.

Convention (columns, as PARI qfauto / qfisom and rank32.bliss_automorphisms): T maps the
lattice with Gram matrix A onto the lattice with Gram matrix B when T^T A T = B and
|det T| = 1. An automorphism of the lattice with Gram matrix G is the case A = B = G.
The checks are exact (Sage integer matrices); they do not search for anything.

整数の行列が二つの格子の間の等長写像であることの、厳密な検算。

約束（列。PARI の qfauto / qfisom、rank32.bliss_automorphisms と同じ）：T^T A T = B かつ
|det T| = 1 のとき、T は Gram 行列 A の格子を Gram 行列 B の格子に移す。Gram 行列 G の格子の
自己同型は、A = B = G の場合である。
検算は厳密（Sage の整数の行列）で、何も探索しない。
"""
from sage.all import ZZ, matrix


def is_isometry(T, A, B):
    """True when T^T A T = B and |det T| = 1 (T, A, B with integer entries).

    T^T A T = B かつ |det T| = 1 なら True（T、A、B は整数の成分）。"""
    T, A, B = matrix(ZZ, T), matrix(ZZ, A), matrix(ZZ, B)
    return T.transpose() * A * T == B and abs(T.det()) == 1


def is_automorphism(g, G):
    """True when g^T G g = G and |det g| = 1.

    g^T G g = G かつ |det g| = 1 なら True。"""
    return is_isometry(g, G, G)


def all_automorphisms(gens, G):
    """True when every matrix of gens is an automorphism of G.

    gens の行列がすべて G の自己同型なら True。"""
    return all(is_automorphism(g, G) for g in gens)
