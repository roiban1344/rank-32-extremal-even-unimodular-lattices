"""A faithful permutation representation of Aut L on vectors of norm 4, for GAP.

perm_model turns a list of generators of Aut L (integer matrices g with g^T G g = G) into a
permutation group of GAP. The group acts on the vectors x and -x of norm 4 (not on the
lines +-x: on lines the sign changes of an orthogonal frame can act trivially, so the action
on lines need not be faithful). The action is restricted to a union of small orbits that
spans Q^n, which is faithful. Vectors are matched by their rows themselves (the bytes of the
int64 rows), so no two vectors can share a key.
perm_to_row turns a permutation of those points back into the row action S (x -> x S) on
the coordinates of L. With these, GAP computes conjugacy classes, centralizers and the like
inside Aut L.

ノルム 4 のベクトルの上の Aut L の忠実な置換表現（GAP 用）。

perm_model は、Aut L の生成元の列（g^T G g = G となる整数の行列 g）を、GAP の置換群にする。
群はノルム 4 のベクトル x と -x の上に作用する（直線 +-x の上ではない。直線の上では、直交する
枠の符号の反転が自明に働きうるので、忠実とは限らない）。作用は、Q^n を張る小さい軌道の和に
制限する。張っていれば忠実である。ベクトルは行そのもの（int64 の行のバイト列）で照合するので、
二つのベクトルのキーが重なることはない。
perm_to_row は、その点の上の置換を、L の座標の行の作用 S（x -> x S）に戻す。これらを使って、
GAP が Aut L の中の共役類や中心化群などを求める。
"""
import numpy as np
from sage.all import ZZ, QQ, matrix, libgap

from .forms import short_vectors
from .order3_gluing import orbits


def perm_model(G, gens, extra=()):
    """The permutation group of GAP of Aut L on norm-4 vectors (generators in the column convention g^T G g = G).

    ノルム 4 のベクトルの上の Aut L の GAP の置換群（生成元は列の約束 g^T G g = G）。
    Returns (K, the permutations of the row actions in extra); with extra = "points", (K, the points).
    (K, extra の行の作用の置換の列) を返す。extra = "points" なら (K, 点) を返す。
    """
    V = short_vectors(G, 4)
    V = np.concatenate([V, -V])
    # The key is the row itself (the bytes of the int64 row), so keys never collide.
    # キーは行そのもの（int64 の行のバイト列）なので、衝突しない
    def rowkeys(A):
        A = np.ascontiguousarray(A, dtype=np.int64)
        return A.view(np.dtype((np.void, A.dtype.itemsize * A.shape[1]))).ravel()
    keys = rowkeys(V)
    order = np.argsort(keys, kind="stable")
    ks = keys[order]
    assert (ks[1:] != ks[:-1]).all(), "duplicate norm-4 vectors"
    vmax = int(np.abs(V).max())
    def perm(R):
        Rm = np.array([[int(x) for x in r] for r in R.rows()], dtype=np.int64)
        assert vmax * int(np.abs(Rm).max()) * Rm.shape[0] < 2 ** 62        # V @ Rm does not overflow int64 / V @ Rm が int64 であふれない
        k = rowkeys(V @ Rm)
        pos = np.minimum(np.searchsorted(ks, k), len(ks) - 1)
        assert (ks[pos] == k).all()
        return order[pos]
    perms = [perm(g.transpose()) for g in gens]
    label, reps, sizes = orbits(perms, len(V))
    chosen = []
    rank = 0
    for i in np.argsort(sizes, kind="stable"):
        mem = np.nonzero(label == reps[i])[0]
        r = matrix(QQ, V[np.concatenate(chosen + [mem])[:4000]].tolist()).rank()
        if r > rank:
            chosen.append(mem)
            rank = r
        if rank == G.nrows():
            break
    assert rank == G.nrows()
    pts = np.concatenate(chosen)
    ix = np.full(len(V), -1, dtype=np.int64)
    ix[pts] = np.arange(len(pts))
    def restrict(p):
        q = ix[p[pts]]
        assert (q >= 0).all()
        return libgap.PermList([int(x) + 1 for x in q])
    K = libgap.Group([restrict(p) for p in perms])
    if extra == "points":
        return K, V[pts]
    return K, [restrict(perm(S)) for S in extra]


def perm_to_row(x, pts):
    """The row action S (x -> x S) on the coordinates of L of a permutation x of the points pts.

    点 pts の上の置換 x を、L の座標の行の作用 S（x -> x S）に戻す。"""
    n = pts.shape[1]
    M = matrix(ZZ, pts.tolist())
    piv = M.transpose().pivots()
    assert len(piv) == n
    idx = [int(i) for i in M.pivot_rows()]
    im = np.array([pts[int(libgap.OnPoints(i + 1, x)) - 1] for i in idx], dtype=np.int64)
    S = matrix(QQ, pts[idx].tolist()).inverse() * matrix(QQ, im.tolist())   # A S = image for the basis rows A / 基底の行 A について A S = 像
    assert S.denominator() == 1
    return S.change_ring(ZZ)
