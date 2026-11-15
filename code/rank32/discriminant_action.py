"""The action of automorphisms of a 2-elementary lattice on its discriminant group, and the components of the involution gluing.

Used by the involution gluing (computation C4 of the first paper).
disc_basis(G) chooses the basis of G^#/G in the same way as rank32.antiisometry_sampler.disc
(columns of G^{-1}, kept when they raise the rank mod 2), and action(G, g, ...) writes the action of
an automorphism g (g^T G g = G, acting from the left on column vectors of coordinates) as the
list "mask of the image of basis vector i".
PARTS names the components F, K of the pairs: name -> (file of the list of 2-elementary lattices of
minimum 4, index of the class in that file).

2-初等的な格子の自己同型の、判別式群への作用と、対合の貼り合わせの成分。

対合の貼り合わせ（第 1 の論文の計算 C4）が使う。
disc_basis(G) は、G^#/G の基底を rank32.antiisometry_sampler.disc と同じ方法で選ぶ（G^{-1} の列を、
mod 2 の階数が増えるときに残す）。action(G, g, ...) は、自己同型 g（g^T G g = G。座標の列ベクトルに
左から作用する）の作用を「基底 i の像の mask」の列で書く。
PARTS は組の成分 F, K の名前である：名前 -> (最小 4 の 2-初等的な格子の一覧のファイル、その中の類の番号)。
"""
from sage.all import GF, matrix, vector

PARTS = {"14": ("min4-twoel-14.json", 0), "15": ("min4-twoel-15.json", 0), "16d8": ("min4-twoel-16.json", 2), "16s14": ("min4-twoel-16.json", 4),
         "17s15": ("min4-twoel-17-s15.json", 1), "18a": ("min4-twoel-18-s14.json", 8), "18b": ("min4-twoel-18-s14.json", 9),
         "18c": ("min4-twoel-18-s14.json", 10), "18d": ("min4-twoel-18-s14.json", 11), "18e": ("min4-twoel-18-s14.json", 12)}
PARTS.update({f"16c{i}": ("min4-twoel-16.json", i) for i in (0, 1, 3, 5, 6, 7, 8, 9)})


def disc_basis(G):
    """Basis of G^#/G (rational column vectors), the pivot columns J and the matrix Sinv that reads coordinates.

    G^#/G の基底（有理の列ベクトル）、ピボットの列 J、座標を読む行列 Sinv。"""
    d = G.nrows(); Xi = G.inverse(); gens = []; span = matrix(GF(2), 0, d)
    for j in range(d):
        c = Xi.column(j); v2 = vector(GF(2), [int(2 * x) % 2 for x in c])
        if span.stack(v2).rank() > span.rank():
            gens.append(c); span = span.stack(v2)
    J = span.pivots(); Sinv = span.matrix_from_columns(J).inverse()
    return gens, J, Sinv


def action(G, g, gens, J, Sinv):
    """Action of g (g^T G g = G, acting from the left on column vectors of coordinates) on the discriminant group: the list of the masks of the images of the basis vectors.

    g（g^T G g = G、座標の列ベクトルに左から作用）の判別式への作用：基底 i の像の mask の列。"""
    out = []
    for c in gens:
        v = vector(GF(2), [int(2 * x) % 2 for x in g * c])
        w = vector(GF(2), [v[j] for j in J]) * Sinv
        out.append(sum(int(w[l]) << l for l in range(len(gens))))
    return out
