"""Gluing for elements of order 3: double cosets of the anti-isometries phi of discriminant forms on F_3^s.

Notation: an element of a discriminant group is a row vector a in F_3^s, and a linear map is a
matrix T acting as a -> a T (from the right). B is 3(y_i, y_j) mod 3 (y_i lifts of the basis).
phi: D_F -> D_Z is an anti-isometry if and only if T B_Z T^t = -B_F (mod 3).
A "line" is a pair +-a ((3^s - 1)/2 of them, 0 left out).

The set X = Iso / bar G (bar G the image of the group of one component in its discriminant group) is
represented by the O(q)-orbit of a colouring of a bar G-invariant set of lines. The orbits of the
group of the other component on X are the double cosets, and
|Stab(phi)| = |group| |kernel| / |orbit| (the kernel is that of Aut -> O(q) on the bar G side).

位数 3 の貼り合わせ：判別式 F_3^s の上の反等長 phi の両側剰余類を扱う道具。

記法：判別式の元は行ベクトル a in F_3^s、線形写像は行列 T で a -> a T（右から）。B は
3(y_i, y_j) mod 3（y_i は基底の持ち上げ）。phi: D_F -> D_Z が反等長であることと、
T B_Z T^t = -B_F (mod 3) は同値である。「直線」は +-a の組（0 を除く (3^s - 1)/2 個）。

点の集合 X = Iso / bar G（bar G は片側の成分の群の判別式への像）は、bar G で不変な直線の集合の彩色の
O(q) の軌道で表す。もう片方の群が X に働く軌道が両側剰余類で、|Stab(phi)| = |群| |カーネル| / |軌道|
である（カーネルは bar G の側の Aut -> O(q) のカーネル）。
"""
import numpy as np

P = 3


def vec_table(s):
    """All elements of F_3^s (a 3^s x s int8 array); the digits of the index i are the coordinates, the first coordinate most significant.

    F_3^s の全元（3^s x s の int8）。添字 i の桁が座標で、先頭の座標が最上位。"""
    n = 3 ** s
    V = np.zeros((n, s), dtype=np.int8)
    t = np.arange(n)
    for i in range(s - 1, -1, -1):
        V[:, i] = t % 3
        t //= 3
    return V


def vec_index(a):
    """The index of an array of row vectors (... x s).

    行ベクトルの配列（... x s）の添字。"""
    a = np.asarray(a, dtype=np.int64)
    s = a.shape[-1]
    return a @ (3 ** np.arange(s - 1, -1, -1, dtype=np.int64))


def inv_mod3(M):
    """The inverse matrix over F_3.

    F_3 上の逆行列。"""
    M = np.asarray(M, dtype=np.int64) % 3
    n = len(M)
    A = np.concatenate([M, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        r = c + int(np.argmax(A[c:, c] != 0))
        if A[r, c] == 0:
            raise ValueError("singular matrix")
        if r != c:
            A[[c, r]] = A[[r, c]]
        A[c] = A[c] * A[c, c] % 3   # 1^-1 = 1, 2^-1 = 2
        for i in range(n):
            if i != c and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[c]) % 3
    return A[:, n:]


def nullspace_mod3(M, s):
    """A basis (as rows) of {x : M x^t = 0}; M is k x s.

    {x : M x^t = 0} の基底（行の列）。M は k x s。"""
    if M is None or len(M) == 0:
        return np.eye(s, dtype=np.int64)
    A = np.array(M, dtype=np.int64) % 3
    r = 0
    piv = []
    for c in range(s):
        nz = np.nonzero(A[r:, c])[0] if r < len(A) else []
        if len(nz) == 0:
            continue
        i = r + int(nz[0])
        A[[r, i]] = A[[i, r]]
        A[r] = A[r] * A[r, c] % 3
        for j in range(len(A)):
            if j != r and A[j, c]:
                A[j] = (A[j] - A[j, c] * A[r]) % 3
        piv.append(c)
        r += 1
    free = [c for c in range(s) if c not in piv]
    out = []
    for f in free:
        v = np.zeros(s, dtype=np.int64)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-A[i, f]) % 3
        out.append(v)
    return np.array(out, dtype=np.int64).reshape(-1, s)


def orthogonal_basis(B, rng):
    """An orthogonal basis V and diagonal d (V B V^t = diag(d)) of a symmetric nondegenerate B over F_3.

    対称な B（F_3 上の非退化）の直交基底 V と対角 d（V B V^t = diag(d)）。"""
    s = len(B)
    V = []
    for _ in range(s):
        M = (np.array(V, dtype=np.int64) @ B) % 3 if V else None
        N = nullspace_mod3(M, s)
        while True:
            c = rng.integers(0, 3, size=len(N))
            v = (c @ N) % 3
            if int(v @ B @ v) % 3:
                V.append(v)
                break
    V = np.array(V, dtype=np.int64)
    return V, np.diag(V @ B @ V.T) % 3


def anti_isometry(BF, BZ, rng):
    """One T with T B_Z T^t = -B_F (mod 3).

    T B_Z T^t = -B_F (mod 3) となる T を 1 つ作る。"""
    s = len(BF)
    for _ in range(1000):
        VF, dF = orthogonal_basis(BF, rng)
        VZ, dZ = orthogonal_basis(BZ, rng)
        # order VZ so that the diagonal values match up to the factor -1; if the counts differ, draw other orthogonal bases
        # 対角の値が -1 倍で揃うように、VZ の並びを選ぶ。個数が合わなければ別の直交基底を引く
        want = (-dF) % 3
        order = []
        pool = {1: [i for i in range(s) if dZ[i] == 1], 2: [i for i in range(s) if dZ[i] == 2]}
        ok = True
        for w in want:
            if not pool[int(w)]:
                ok = False
                break
            order.append(pool[int(w)].pop())
        if not ok:
            continue
        T = inv_mod3(VF) @ VZ[order] % 3
        assert ((T @ BZ @ T.T + BF) % 3 == 0).all()
        return T
    raise ArithmeticError("no anti-isometry found (the discriminant forms do not match)")


def reflection(B, v):
    """The matrix of a -> a - 2 (a,v)/(v,v) v (acting on row vectors from the right).

    a -> a - 2 (a,v)/(v,v) v の行列（行ベクトルに右から）。"""
    v = np.asarray(v, dtype=np.int64)
    q = int(v @ B @ v) % 3
    assert q
    s = len(B)
    # -2/q = 1/q = q (q^{-1} = q in F_3) / -2/q = 1/q = q（F_3 で q^{-1} = q）
    return (np.eye(s, dtype=np.int64) + q * np.outer(B @ v % 3, v)) % 3


def line_tables(s):
    """Return (the line number of every vector (-1 for 0), the indices of the representatives of the lines); a representative has first nonzero coordinate 1.

    (各ベクトルの直線の番号（0 は -1）, 直線の代表の添字) を返す。代表は、先頭の 0 でない座標が 1 のもの。"""
    V = vec_table(s)
    n = 3 ** s
    first = np.argmax(V != 0, axis=1)
    lead = V[np.arange(n), first]
    reps = np.nonzero(lead == 1)[0]
    line = np.full(n, -1, dtype=np.int32)
    line[reps] = np.arange(len(reps))
    other = vec_index((2 * V[reps].astype(np.int64)) % 3)
    line[other] = np.arange(len(reps))
    return V, line, reps


class Space:
    """F_3^s with the form B; builds the tables of permutations of the lines.

    F_3^s と形式 B。直線の上の置換の表を作る。"""

    def __init__(self, B):
        self.B = np.asarray(B, dtype=np.int64) % 3
        self.s = len(self.B)
        self.V, self.line, self.reps = line_tables(self.s)
        self.nlines = len(self.reps)
        self.Vr = self.V[self.reps].astype(np.int64)

    def line_perm(self, T):
        """The permutation of the lines induced by the matrix T (a -> a T).

        行列 T（a -> a T）の直線の上の置換。"""
        return self.line[vec_index((self.Vr @ T) % 3)]

    def lines_of(self, vecs):
        """The set (ascending) of line numbers of an array of vectors.

        ベクトルの配列の直線の番号の集合（昇順）。"""
        vecs = np.asarray(vecs, dtype=np.int64).reshape(-1, self.s)
        return np.unique(self.line[vec_index(vecs % 3)])

    def q_lines(self):
        """The value of q on every line (aBa mod 3: 0, 1 or 2).

        各直線の q の値（aBa mod 3、0, 1, 2）。"""
        return np.einsum("ij,jk,ik->i", self.Vr, self.B, self.Vr) % 3


class FrameSet:
    """The O(q)-orbit of a coloured set of lines. A point is the concatenation of the ascending lists of line numbers of each colour.

    sets = [S_1, S_2, ...] (line numbers per colour). Points are looked up by a Zobrist sum
    (64 bits), and after the lookup the contents are checked to agree. The orbit under the
    reflections gens is built breadth first. When its size equals |O(q)|/|bar G|, the stabilizer of
    the colouring is exactly bar G and the orbit of gens is everything.

    彩色された直線の集合の O(q) の軌道。点は「色ごとの直線の番号の昇順の列」を連ねた配列で持つ。

    sets = [S_1, S_2, ...]（色ごとの直線の番号）。点は Zobrist 和（64 ビット）で引き、引いたあとに
    実際の内容も一致することを確かめる。幅優先で反射 gens の軌道を作る。軌道の大きさが
    |O(q)|/|bar G| に等しいとき、彩色の安定化群がちょうど bar G で、gens の軌道が全部である。"""

    def __init__(self, space, sets, gens, expected=None, seed=20261031, limit=None):
        self.space = space
        rng = np.random.default_rng(seed)
        self.ncol = len(sets)
        self.zob = rng.integers(1, 2 ** 63, size=(space.nlines, self.ncol), dtype=np.uint64)
        self.sizes = [len(x) for x in sets]
        self.gens = [np.asarray(g, dtype=np.int64) for g in gens]
        perms = [space.line_perm(g) for g in self.gens]
        dtype = np.uint16 if space.nlines < 65536 else np.int32
        pts = np.concatenate([np.sort(np.asarray(x, dtype=np.int64)) for x in sets]).astype(dtype)[None, :]
        keys = self.key(pts)
        allkeys = keys.copy()     # ascending / 昇順
        allpos = np.zeros(1, dtype=np.int64)
        parents = [np.array([0], dtype=np.int64)]
        pgens = [np.array([-1], dtype=np.int32)]
        blocks = [pts]
        cnt = 1
        front, fstart = pts, 0
        while len(front):
            newk, newpts, newpar, newgen = [], [], [], []
            for gi, perm in enumerate(perms):
                im = self.apply(perm, front.astype(np.int64))
                k = self.key(im)
                pos = np.searchsorted(allkeys, k)
                pos[pos == len(allkeys)] = 0
                fresh = allkeys[pos] != k
                ii = np.nonzero(fresh)[0]
                if len(ii):
                    newk.append(k[ii]); newpts.append(im[ii].astype(dtype)); newpar.append(fstart + ii); newgen.append(np.full(len(ii), gi, dtype=np.int32))
            if not newk:
                break
            k = np.concatenate(newk)
            uk, first = np.unique(k, return_index=True)
            front = np.concatenate(newpts)[first]
            npar = np.concatenate(newpar)[first]; ngen = np.concatenate(newgen)[first]
            blocks.append(front); parents.append(npar); pgens.append(ngen)
            fstart = cnt
            pos2 = np.arange(cnt, cnt + len(uk))
            cnt += len(uk)
            ins = np.searchsorted(allkeys, uk)
            allkeys = np.insert(allkeys, ins, uk)
            allpos = np.insert(allpos, ins, pos2)
            if limit is not None and cnt > limit:
                raise OverflowError(f"the number of points exceeds the limit {limit}")
        self.pts = np.concatenate(blocks)
        self.parent = np.concatenate(parents)
        self.pgen = np.concatenate(pgens)
        self.keys_sorted = allkeys
        self.pos_sorted = allpos
        self.n = cnt
        # that distinct points do not share a key is checked by the comparison of contents in perm_of and by bijectivity
        # 異なる点が同じキーに落ちていないことは、perm_of の内容の照合と、置換が全単射であることで確かめる
        if expected is not None and cnt != expected:
            raise ArithmeticError(f"orbit size {cnt} differs from the expected {expected}")

    def split(self, arr):
        """Split an array of points into the parts of each colour.

        点の配列を色ごとの部分に分ける。"""
        out = []
        i = 0
        for z in self.sizes:
            out.append(arr[:, i:i + z])
            i += z
        return out

    def key(self, arr):
        """The Zobrist sums of an array of points.

        点の配列の Zobrist 和。"""
        k = np.zeros(len(arr), dtype=np.uint64)
        for c, part in enumerate(self.split(arr)):
            k += self.zob[part.astype(np.int64), c].sum(axis=1, dtype=np.uint64)
        return k

    def apply(self, perm, arr):
        """Apply a permutation of the lines to an array of points (re-sorted ascending per colour).

        点の配列に直線の置換を施す（色ごとに昇順に並べ直す）。"""
        parts = [np.sort(perm[p.astype(np.int64)], axis=1) for p in self.split(arr)]
        return np.concatenate(parts, axis=1)

    def perm_of(self, T):
        """The permutation of X induced by the matrix T (i -> number of the image of point i); raises if an image leaves X or the contents disagree.

        行列 T の X の上の置換（i -> 点 i の像の番号）。X の外に出たり、内容が食い違えば例外。"""
        im = self.apply(self.space.line_perm(np.asarray(T, dtype=np.int64)), self.pts.astype(np.int64))
        k = self.key(im)
        pos = np.searchsorted(self.keys_sorted, k)
        pos[pos == self.n] = 0
        if (self.keys_sorted[pos] != k).any():
            raise ArithmeticError("left X (not an element of O(q), or X is not closed under O(q))")
        out = self.pos_sorted[pos]
        if not (self.pts[out].astype(np.int64) == im).all():
            raise ArithmeticError("collision of Zobrist sums")
        if len(np.unique(out)) != self.n:
            raise ArithmeticError("the permutation is not a bijection")
        return out

    def element(self, i):
        """An element h of O(q) (a matrix) with point i = image of point 0: multiplying the lines of point 0 by h from the right gives point i.

        点 i = 点 0 の像とする O(q) の元 h（行列）。点 0 の直線に右から h を掛けると点 i。"""
        word = []
        while self.pgen[i] >= 0:
            word.append(int(self.pgen[i]))
            i = int(self.parent[i])
        h = np.eye(self.space.s, dtype=np.int64)
        for g in reversed(word):
            h = (h @ self.gens[g]) % 3
        return h


def orbits(perms, n):
    """Orbits of a list of permutations: (label (the least point of the orbit of each point), array of representatives, array of orbit sizes).

    置換の列の軌道。(ラベル（各点の軌道の最小の点）, 代表の配列, 軌道の大きさの配列) を返す。"""
    label = np.arange(n)
    while True:
        new = label.copy()
        for t in perms:
            np.minimum.at(new, t, label)
            new = np.minimum(new, new[t])
        new = new[new]
        if (new == label).all():
            break
        label = new
    reps, sizes = np.unique(label, return_counts=True)
    return label, reps, sizes


def assert_form_preserved(T, B):
    """Assert T B T^t = B (mod 3).

    T B T^t = B (mod 3) を確かめる。"""
    T = np.asarray(T, dtype=np.int64)
    assert ((T @ B @ T.T - B) % 3 == 0).all()


def reflection_gens(B, rng, count=24):
    """A list of random reflections meant to generate O(q) (whether they do is checked by the orbit size).

    O(q) を生成するつもりの反射の列（生成するかどうかは軌道の大きさで確かめる）。"""
    s = len(B)
    out = []
    while len(out) < count:
        v = rng.integers(0, 3, size=s)
        if int(v @ B @ v) % 3:
            out.append(reflection(B, v))
    return out


def config_edges(vecs, colors, B):
    """The edges (i, j, c) of the graph that tests the coloured configuration of a set vecs of elements of D (rows, closed under +-) up to O(q) (the input of rank32/_canonical.c).

    Colour 1 means (x, y) = 1, colour 2 means y = -x, and 10 + the vertex colour is a loop.
    (x, y) = 2 is determined by (-x, y) = 1, so an isomorphism of the graph is a map of vertices
    keeping the inner products, and if vecs spans D it extends uniquely to an element of O(q).

    D の元の集合 vecs（行、+- で閉じる）の色つき配置を O(q) の同型で調べるグラフの辺 (i, j, c)（rank32/_canonical.c の入力）。

    色 1 は (x, y) = 1、色 2 は y = -x、10 + 頂点の色は自己ループ。(x, y) = 2 は (-x, y) = 1 で
    決まるので、グラフの同型は内積を保つ頂点の写像で、vecs が D を張れば O(q) の元に一意に延びる。"""
    n = len(vecs)
    M = (vecs @ B @ vecs.T) % 3
    ii, jj = np.nonzero(M == 1)
    pos = {int(k): i for i, k in enumerate(vec_index(vecs))}
    neg = np.array([pos[int(k)] for k in vec_index((-vecs) % 3)], dtype=np.int64)
    loop = np.arange(n, dtype=np.int64)
    return np.vstack([np.column_stack([ii, jj, np.ones_like(ii)]),
                      np.column_stack([loop, neg, np.full(n, 2, dtype=np.int64)]),
                      np.column_stack([loop, loop, 10 + np.asarray(colors, dtype=np.int64)])]).astype(np.int32)
