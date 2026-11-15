"""Even 2-neighbours of extremal even unimodular lattices of rank 32, and the replay of a step of the seeded walk.

Paper: used by the checks of the tools (computation C37 of the first paper): the check of the bliss
automorphism tool on the lattices of a random walk (the seven repeated fingerprints of the walk are
rebuilt by replay and compared by isometry tests).
Facts behind the step: let L be even unimodular with min L >= 4 and c in L/2L a nonzero class whose
norms are divisible by 4 (an admissible class). For a representative v of c with (v, v) = 0 mod 8,
the even 2-neighbour L'_c = L_v + Z v/2 (L_v = {x : (x, v) even}) has min L'_c >= 4 if and only if c
contains no vector of norm 8. If c contains k pairs of vectors of norm 8, they are mutually orthogonal
and the root system of L'_c is A_1^k. So the test of a class is one closest-vector enumeration in L
with target -v/2 and squared radius 2, and L is reduced once.
A step of the walk: draw admissible classes uniformly until the type of the class is not ("k", k)
(it is "min4" or "none", so the class has no vector of norm 8), lift it mod 8, build the neighbour and
LLL-reduce it. The random generator of step s of walker w is random.Random("pilot-w-s").
Needs: Sage and fpylll.

階数 32 の極値偶ユニモジュラー格子の偶 2 近傍と、種つきのランダムウォークの 1 歩の再現。

論文：道具の検査（第 1 の論文の計算 C37）が使う。bliss の自己同型の道具をランダムウォークの格子で確かめる検算で、
ランダムウォークで指紋が重複した 7 件を replay で作り直し、同型判定で比べる。
1 歩の背景：L を偶ユニモジュラーで min L >= 4、c in L/2L をノルムが 4 で割れる 0 でない類（許容類）とする。
c の代表 v で (v, v) = 0 mod 8 のものから作る偶 2 近傍 L'_c = L_v + Z v/2（L_v = {x : (x, v) は偶数}）について、
min L'_c >= 4 であることと、c にノルム 8 のベクトルが無いことは同値である。c にノルム 8 のベクトルが k 対あれば、
それらは互いに直交し、L'_c のルート系は A_1^k になる。だから類の検査は、L で目標 -v/2、半径の 2 乗 2 の
最近ベクトルの全列挙を 1 回行うことであり、L の簡約は 1 回で済む。
ランダムウォークの 1 歩：類の型が ("k", k) でなくなる（"min4" か "none" になり、ノルム 8 のベクトルが無い）まで許容類を一様に引き、mod 8 に持ち上げ、
近傍を作って LLL 簡約する。歩行者 w の第 s 歩の乱数は random.Random("pilot-w-s") である。
必要なもの：Sage と fpylll。
"""
import random

from fpylll import BKZ, GSO, LLL, Enumeration, EnumerationError, IntegerMatrix
from fpylll.algorithms.bkz2 import BKZReduction
from sage.all import QQ, ZZ, matrix, vector


def reduce(G, bs=20):
    """(GSO M, transformation T, reduced Gram R = T G T^t) after LLL and, for bs > 2, BKZ with block size bs.

    LLL と（bs > 2 なら）ブロックの大きさ bs の BKZ の後の (GSO M, 変換 T, 簡約後の Gram R = T G T^t)。"""
    n = G.nrows()
    A = IntegerMatrix.from_matrix([[int(x) for x in r] for r in G.rows()])
    U = IntegerMatrix.identity(n)
    M = GSO.Mat(A, U=U, gram=True, float_type="ld")
    M.update_gso()
    LLL.Reduction(M)()
    if bs > 2:
        BKZReduction(M)(BKZ.Param(block_size=bs, max_loops=8, strategies=BKZ.DEFAULT_STRATEGY))
    M.update_gso()
    T = matrix(ZZ, n, n, lambda i, j: U[i, j])
    return M, T, T * G * T.transpose()


def random_admissible_class(G, rng):
    """A uniformly random class of L/2L with norm = 0 mod 4 (not 0), as a 0/1 coefficient vector.

    L/2L の一様な乱択類で、ノルム = 0 mod 4 のもの（0 以外）。0/1 の係数ベクトルを返す。"""
    n = G.nrows()
    while True:
        v = vector(ZZ, [rng.randrange(2) for _ in range(n)])
        if v != 0 and (v * G * v) % 4 == 0:
            return v


def mod8_rep(G, v):
    """A representative of the class v + 2L with norm = 0 mod 8.

    類 v + 2L の代表で、ノルム = 0 mod 8 のもの。"""
    if (v * G * v) % 8 == 4:
        a = G * v
        k = next(i for i in range(len(v)) if a[i] % 2)
        v = v + 2 * vector(ZZ, [1 if i == k else 0 for i in range(len(v))])
    assert (v * G * v) % 8 == 0
    return v


class ClassTester:
    """Reduce L once and list the short vectors (norm <= 8) of classes v + 2L.

    L を 1 回簡約し、類 v + 2L の短いベクトル（ノルム <= 8）を全部挙げる。"""

    def __init__(self, G, bs=20):
        self.G = G
        self.M, self.T, self.R = reduce(G, bs)
        self.Tinv = self.T.inverse()
        self.n = G.nrows()

    def short_in_class(self, v, radius2=2.02, maxsol=400):
        """The norms of all w = v + 2l in v + 2L with norm <= 4*radius2 (checked in integers).
        |w|^2 = 4 |v/2 + l|^2: a full closest-vector enumeration (no pruning) in L with target -v/2.

        類 v + 2L のノルム <= 4*radius2 の元 w = v + 2l のノルムをすべて返す（整数で確かめる）。
        |w|^2 = 4 |v/2 + l|^2。L（簡約基底）で目標 -v/2 の最近ベクトルの全列挙（枝刈りなし）。"""
        n = self.n; M = self.M
        t = -(vector(QQ, v) * self.Tinv) / 2                        # target in the reduced basis / 簡約基底の座標での目標
        tf = [float(x) for x in t]
        c = [tf[i] + sum(M.get_mu(j, i) * tf[j] for j in range(i + 1, n)) for i in range(n)]
        E = Enumeration(M, nr_solutions=maxsol)
        try:
            sols = E.enumerate(0, n, radius2, 0, target=c)
        except EnumerationError:
            return []
        assert len(sols) < maxsol
        out = []
        for _, coords in sols:
            z = vector(ZZ, [int(round(x)) for x in coords]) * self.T   # a point of L in the original coordinates / L の点（元の座標）
            # a point z close to the target -v/2 gives w = v + 2z / 目標 -v/2 に近い点 z から w = v + 2z
            w = vector(ZZ, v) + 2 * z
            nw = int(w * self.G * w)
            out.append(nw)
        return sorted(out)

    def class_type(self, v):
        """('min4', 0): the class contains minimal vectors; ('k', k): k pairs of norm 8; ('none', 0): no element of norm <= 8.

        ('min4', 0)：最短ベクトルの類。('k', k)：ノルム 8 が k 対。('none', 0)：ノルム <= 8 の元なし。"""
        norms = self.short_in_class(v)
        assert all(x in (4, 8) for x in norms), norms
        n4 = norms.count(4); n8 = norms.count(8)
        if n4:
            assert n4 == 2 and n8 == 0, norms
            return "min4", 0
        assert n8 % 2 == 0
        return ("k", n8 // 2) if n8 else ("none", 0)


def neighbour_gram(G, v):
    """The Gram matrix of the even 2-neighbour L' = L_v + Z v/2 (v a representative with (v, v) = 0 mod 8).

    偶 2 近傍 L' = L_v + Z v/2 の Gram（v は (v, v) = 0 mod 8 の代表）。"""
    n = G.nrows(); a = [int(x) % 2 for x in G * v]
    k = next(i for i in range(n) if a[i])
    rows = []
    for i in range(n):
        e = [0] * n
        if i == k:
            e[k] = 2
        else:
            e[i] = 1
            if a[i]:
                e[k] = 1
        rows.append(e)
    Bv = matrix(ZZ, rows)
    H = matrix(ZZ, list(2 * Bv) + [v]).hermite_form(include_zero_rows=False)
    Bp = H.change_ring(QQ) / 2
    Gp = Bp * G * Bp.transpose()
    assert all(x in ZZ for x in Gp.list())
    Gp = Gp.change_ring(ZZ)
    assert Gp.det() == 1 and all(Gp[i, i] % 2 == 0 for i in range(n))
    return Gp


def lll_gram(G):
    """The LLL-reduced Gram matrix.

    LLL 簡約した Gram 行列。"""
    M, T, R = reduce(G, 2)
    return R


def replay(G, w, s):
    """Rebuild step s of walker w from the lattice G of step s - 1 (the random generator is seeded by "pilot-w-s").

    第 s - 1 歩の格子 G から、歩行者 w の第 s 歩を作り直す（乱数の種は "pilot-w-s"）。"""
    rng = random.Random(f"pilot-{w}-{s}")
    ct = ClassTester(G)
    while True:
        v = random_admissible_class(G, rng)
        typ, k = ct.class_type(v)
        if typ != "k":
            break
    v = mod8_rep(G, v)
    return lll_gram(neighbour_gram(G, v))
