"""Rootless mass of the Gaussian Hermitian lattices without enumeration (King's method over Z[i]).

Paper: computation C11 of the first paper (the rootless mass M_2 of the Hermitian lattices of
rank 16 over Z[i] whose trace form is even unimodular); the entry point with the commands is
code/produce/gaussian/rootless_hermitian_mass_gaussian.py.
Formula:
  m(emptyset) = sum_T c(T) A_r(T) / |U(T)|,
  T runs over multisets of irreducible Z[i] root lattices of rank n <= r, and c is the Moebius
  function of the partial order of root sublattices (multiplicative in the irreducible components).
  A_r(T) = S'_t det(T)^t prod_p f_p(T, t), t = r - n, S'_t = 2^{-t [t odd]} S_t,
  S_t = 2^{1-t} prod_{j<=t} |L(1-j, chi^j)|, S'_0 = 1.
  f_p(T, t) = sum_{C subset D_p(T) isotropic} |C|^{-t} prod_{j in J} (1 - chi(p)^j p^{-j}),
    p odd: J = {t-d+1, ..., t}; p = 2: the even j in {t-d+2, ..., t} (chi(2) is taken as 1);
    d is the rank of C^perp/C.
table(r) solves the triangular system for the mass m(R) of every root system R.
The 2-adic statistics are cached in OUT/cache-r<r>.json (JSON; a run can be resumed from it); OUT is
results/C11 under the current directory by default.
Needs: Python 3 and numpy (no Sage).

Z[i] 上のエルミート格子のルート無しの質量を、列挙せずに厳密に求める（King の方法のエルミート版）。

論文：第 1 の論文の計算 C11（トレース形式が偶ユニモジュラーになる、Z[i] 上の階数 16 のエルミート格子の
ルート無しの質量 M_2）。コマンドの入口は code/produce/gaussian/rootless_hermitian_mass_gaussian.py である。
式：
  m(emptyset) = sum_T c(T) A_r(T) / |U(T)|
  T は既約な Z[i] ルート格子の多重集合（階数 n <= r）。c は部分ルート格子の半順序の Moebius 関数
  （既約成分について乗法的）。
  A_r(T) = S'_t det(T)^t prod_p f_p(T, t)、t = r - n、S'_t = 2^{-t [t odd]} S_t、
  S_t = 2^{1-t} prod_{j<=t} |L(1-j, chi^j)|、S'_0 = 1。
  f_p(T, t) = sum_{C subset D_p(T) 等方} |C|^{-t} prod_{j in J} (1 - chi(p)^j p^{-j})
    p 奇数：J = {t-d+1, ..., t}。p = 2：J = {t-d+2, ..., t} の偶数（chi(2) は 1 とみなす）。
    d は C^perp/C の階数。
table(r) は三角の連立を解いて、各ルート系 R の質量 m(R) を求める。
2 進の統計は OUT/cache-r<r>.json（JSON）にキャッシュする。途中から再開できる。OUT の既定は今のディレクトリの
下の results/C11 である。
必要なもの：Python 3 と numpy（Sage は要らない）。
"""
import itertools
import json
import os
import time
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

from .io import write_json_atomic_indent1 as _write_json

# Output directory of solve() and table(); the command-line entry sets it / solve() と table() の出力先。コマンドの入口が設定する
OUT = Path("results/C11")

# ---------------------------------------------------------------- L-values / L 値


def bernoulli_list(n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


def euler_abs(n):
    """|E_{2k}| (secant numbers); L(1-j, chi_{-4}) = E_{j-1}/2 for odd j.

    |E_{2k}|（セカント数）。j が奇数のとき L(1-j, chi_{-4}) = E_{j-1}/2。"""
    E = [Fr(0)] * (n + 1)
    E[0] = Fr(1)
    for m in range(2, n + 1, 2):
        E[m] = -sum(comb(m, k) * E[k] for k in range(0, m, 2))
    return [abs(x) for x in E]


_B = bernoulli_list(40)
_E = euler_abs(40)


def L_abs(j):
    if j % 2:
        return _E[j - 1] / 2
    return abs(_B[j]) / j


def S_mass(t):
    """Standard mass S_t = 2^{1-t} prod_{j<=t} |L(1-j, chi^j)| (S_0 = 1 by convention).

    標準質量 S_t = 2^{1-t} prod_{j<=t} |L(1-j, chi^j)|（t = 0 では 1 と約束する）。"""
    if t == 0:
        return Fr(1)
    m = Fr(1, 2 ** (t - 1))
    for j in range(1, t + 1):
        m *= L_abs(j)
    return m


def Sprime(t):
    return S_mass(t) * (Fr(1, 2 ** t) if t % 2 else 1)


# ---------------------------------------------------------------- Gaussian integers (denominators powers of 2 allowed) / ガウス整数（2 冪の分母を許す）
# An element is a + b i with Fraction parts; matrices are computed directly with these pairs.
# 要素は実部と虚部の Fraction の組。行列の計算はこの組で素直に行う。


class G:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fr(a)
        self.b = Fr(b)

    def __add__(self, o):
        return G(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        return G(self.a - o.a, self.b - o.b)

    def __neg__(self):
        return G(-self.a, -self.b)

    def __mul__(self, o):
        return G(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a)

    def conj(self):
        return G(self.a, -self.b)

    def norm(self):
        return self.a * self.a + self.b * self.b

    def inv(self):
        n = self.norm()
        return G(self.a / n, -self.b / n)

    def iszero(self):
        return self.a == 0 and self.b == 0

    def __repr__(self):
        return f"({self.a}+{self.b}i)"

    def __eq__(self, o):
        return self.a == o.a and self.b == o.b


ONE = G(1)
ZERO = G(0)
PI = G(1, 1)
PIB = G(1, -1)


def v2(x: Fr):
    """2-adic valuation of x != 0.

    2 進付値。x != 0。"""
    n, d = x.numerator, x.denominator
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    while d % 2 == 0:
        d //= 2
        v -= 1
    return v


def vpi(z: G):
    """pi-adic valuation (pi = 1 + i) of z != 0: v_pi(z) = v_2(N(z)).

    pi 進付値（pi = 1 + i）。z != 0。v_pi(z) = v_2(N(z))。"""
    return v2(z.norm())


def pi_pow(s):
    """pi^s (s may be negative).

    pi^s（s は負でもよい）。"""
    out = ONE
    if s >= 0:
        for _ in range(s):
            out = out * PI
    else:
        inv = G(Fr(1, 2), Fr(-1, 2))
        for _ in range(-s):
            out = out * inv
    return out


# ---------------------------------------------------------------- 2-adic Jordan blocks / 2 進の Jordan ブロック
# ('d', e, u): <2^e u> (e >= 1, u = 1, 3 mod 4). D = O/2^e, generator y = e_1/2^e, q(y) = u/2^e.
# ('H', e): 2^e [[0,1],[1,0]] (e >= 1). D = (O/2^e)^2.
# ('P', e, t): the plane of type pi^{2e+1}, [[2^{e+1} t, pi^{2e+1}], [conj, 2^{e+1} t]] (t = 0, 1). D = (O/pi^{2e+1})^2.
# 上の 3 種：対角の <2^e u>、双曲平面 2^e H、pi^{2e+1} 型の平面。D はそれぞれの判別式群。


def block_gram(b):
    if b[0] == "d":
        return [[G(2 ** b[1] * b[2])]]
    if b[0] == "H":
        x = G(2 ** b[1])
        return [[ZERO, x], [x, ZERO]]
    _, e, t = b
    s = pi_pow(2 * e + 1)
    a = G(2 ** (e + 1) * t)
    return [[a, s], [s.conj(), a]]


def block_rank(b):
    return 1 if b[0] == "d" else 2


def block_denom_pi(b):
    """The coordinates of D are pi^{-s} c (c in a system of representatives of O mod pi^s); return s.

    D の座標は pi^{-s} c（c は O の代表系）で、s を返す。"""
    if b[0] == "d":
        return 2 * b[1]
    if b[0] == "H":
        return 2 * b[1]
    return 2 * b[1] + 1


def residues_pi(s):
    """A system of representatives of O / pi^s.

    O / pi^s の代表系。"""
    # pi^s O has index 2^s; choose among a + b i, 0 <= a, b < 2^{ceil(s/2)}, those with pairwise differences not divisible by pi^s
    # pi^s O の指数は 2^s。a + b i（0 <= a, b < 2^{ceil(s/2)}）から、差が pi^s で割れないものを選ぶ
    m = 2 ** ((s + 1) // 2)
    reps = []
    pis = pi_pow(s)
    pinv = pis.inv()
    for a in range(m):
        for b in range(m):
            z = G(a, b)
            ok = True
            for w in reps:
                q = (z - w) * pinv
                if q.a.denominator == 1 and q.b.denominator == 1:
                    ok = False
                    break
            if ok:
                reps.append(z)
    assert len(reps) == 2 ** s
    return reps


def ann_of(vec_pi_denoms):
    """The exponent of the order of x (x_i = pi^{-s} c_i): the least c with pi^c x integral.

    x_i = pi^{-s} c_i の x の位数の指数：pi^c x が整になる最小の c。"""
    m = 0
    for z in vec_pi_denoms:
        if not z.iszero():
            m = max(m, -vpi(z))
    return m


def herm(x, A, y):
    """h(x, y) = sum x_i A_ij conj(y_j), linear in the first variable.

    h(x, y) = sum x_i A_ij conj(y_j)。第 1 変数について線形。"""
    out = ZERO
    for i, xi in enumerate(x):
        if xi.iszero():
            continue
        for j, yj in enumerate(y):
            if yj.iszero() or A[i][j].iszero():
                continue
            out = out + xi * A[i][j] * yj.conj()
    return out


def block_autos(b, limit=4000):
    """Automorphisms of D(b) (O-linear, preserving q) as permutations. Not necessarily all of them: a subgroup only
    makes the orbits finer, and the count stays correct. Returns (list of coordinates of the elements, list of permutations).

    D(b) の自己同型（O 線形で q を保つ）を置換として返す。全部とは限らない（部分群では軌道が細かくなるだけで、数え方は正しいまま）。
    返り値：(元の座標のリスト, 置換のリスト)。"""
    A = block_gram(b)
    s = block_denom_pi(b)
    R = _res(s)
    pinv = pi_pow(-s)
    pis = pi_pow(s)

    def red1(z):
        c = z * pis
        for j, w in enumerate(R):
            q = (c - w) * pinv
            if q.a.denominator == 1 and q.b.denominator == 1:
                return j
        raise AssertionError

    rc = {}

    def red(z):
        k = (z.a, z.b)
        if k not in rc:
            rc[k] = red1(z)
        return rc[k]

    if block_rank(b) == 1:
        els = [[c * pinv] for c in R]
        idx = {(red(x[0]),): i for i, x in enumerate(els)}
        q = [herm(x, A, x).a % 2 for x in els]
        g0 = idx[(red(pinv),)]
        perms = []
        for i in range(len(els)):
            if ann_of(els[i]) == s and q[i] == q[g0]:
                c = els[i][0] * pis
                perms.append([idx[(red(c * x[0]),)] for x in els])
        return els, perms
    els = [[c1 * pinv, c2 * pinv] for c1 in R for c2 in R]
    idx = {(red(x[0]), red(x[1])): i for i, x in enumerate(els)}
    q = [herm(x, A, x).a % 2 for x in els]
    ann = [ann_of(x) for x in els]
    g1 = idx[(red(pinv), red(ZERO))]
    g2 = idx[(red(ZERO), red(pinv))]
    b0 = herm(els[g1], A, els[g2])
    gens1 = [i for i in range(len(els)) if ann[i] == s and q[i] == q[g1]]
    gens2 = [i for i in range(len(els)) if ann[i] == s and q[i] == q[g2]]
    perms = []
    tried = 0
    for i in gens1:
        for j in gens2:
            z = herm(els[i], A, els[j]) - b0
            if z.a.denominator != 1 or z.b.denominator != 1:
                continue
            tried += 1
            if tried > limit:
                return els, perms
            perm = [idx[(red(x[0] * pis * els[i][0] + x[1] * pis * els[j][0]), red(x[0] * pis * els[i][1] + x[1] * pis * els[j][1]))] for x in els]
            if len(set(perm)) == len(els):
                perms.append(perm)
    return els, perms


@lru_cache(maxsize=None)
def block_elements(b):
    """The elements of D of a block, split into orbits of automorphisms: {(ann, q, index): (size, coordinates of a representative)}.

    ブロックの D の元を自己同型の軌道で分けた類：{(ann, q, 番号): (個数, 代表の座標)}。"""
    A = block_gram(b)
    els, perms = block_autos(b)
    parent = list(range(len(els)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for perm in perms:
        for a, c in enumerate(perm):
            parent[find(a)] = find(c)
    orb = defaultdict(list)
    for i in range(len(els)):
        orb[find(i)].append(i)
    out = {}
    for k, (root, mem) in enumerate(sorted(orb.items())):
        x = els[mem[0]]
        out[(ann_of(x), herm(x, A, x).a % 2, k)] = (len(mem), x)
    return out


# ---------------------------------------------------------------- Jordan decomposition (2-adic, over Z_2[i]; odd denominators are 2-adic units)
# Jordan 分解（2 進、Z_2[i]。奇数の分母は 2 進の単数）


def jordan(A):
    """Jordan decomposition of an even Hermitian matrix A (entries in Z_(2)[i]); returns the list of blocks. Unimodular planes are dropped.

    偶エルミート行列 A（Z_(2)[i] の成分）を Jordan 分解して、ブロックの列を返す。ユニモジュラーの平面は捨てる。"""
    n = len(A)
    A = [row[:] for row in A]
    idx = list(range(n))
    blocks = []
    while idx:
        # find the least scale / 最小のスケールを探す
        best = None
        for i in idx:
            if not A[i][i].iszero():
                v = vpi(A[i][i])
                if best is None or v < best[0] or (v == best[0] and best[1] == "off"):
                    best = (v, "diag", i, i)
        for i in idx:
            for j in idx:
                if i < j and not A[i][j].iszero():
                    v = vpi(A[i][j])
                    if best is None or v < best[0]:
                        best = (v, "off", i, j)
        if best is None:
            raise ValueError("degenerate")
        v, kind, i, j = best
        if kind == "diag":
            a = A[i][i]
            assert a.b == 0 and v % 2 == 0
            e = v // 2
            assert e >= 1, "not even"
            u = a.a / 2 ** e
            un, ud = u.numerator, u.denominator
            um = (un * pow(ud, -1, 4)) % 4
            blocks.append(("d", e, um))
            ainv = a.inv()
            rest = [k for k in idx if k != i]
            newA = {}
            for k in rest:
                for l in rest:
                    newA[(k, l)] = A[k][l] - A[k][i] * ainv * A[i][l]
            for (k, l), val in newA.items():
                A[k][l] = val
            idx = rest
        else:
            B = A[i][j]
            Bb = A[j][i]
            a, c = A[i][i], A[j][j]
            det = a * c - B * Bb
            dinv = det.inv()
            rest = [k for k in idx if k not in (i, j)]
            if v % 2 == 0:
                e = v // 2
                if e > 0:
                    blocks.append(("H", e))
            else:
                e = (v - 1) // 2
                aa = a.a / 2 ** (e + 1)
                cc = c.a / 2 ** (e + 1)
                par = ((aa.numerator * aa.denominator) % 2) * ((cc.numerator * cc.denominator) % 2)
                blocks.append(("P", e, par))
            newA = {}
            for k in rest:
                for l in rest:
                    # A_kl - (A_ki, A_kj) M^{-1} (A_il, A_jl)^T、M^{-1} = dinv [[c, -B],[-Bb, a]]
                    t1 = A[k][i] * c + A[k][j] * (-Bb)
                    t2 = A[k][i] * (-B) + A[k][j] * a
                    val = A[k][l] - (t1 * A[i][l] + t2 * A[j][l]) * dinv
                    newA[(k, l)] = val
            for (k, l), val in newA.items():
                A[k][l] = val
            idx = rest
    return blocks


def normal_blocks(bl):
    """Normalise by the equivalences of blocks (the count is only as correct as these equivalences):
    the sign of P: P(e,1)^2 = P(e,0)^2; d: <2^e 3>^2 = <2^e>^2.

    ブロックの同値で正規化する（この同値が誤っていれば数え方も誤る）。
    P の符号：P(e,1)^2 = P(e,0)^2。d：<2^e 3>^2 = <2^e>^2。"""
    c = Counter(bl)
    out = []
    es = sorted({b[1] for b in c if b[0] == "P"})
    for e in es:
        n1 = c[("P", e, 1)]
        n0 = c[("P", e, 0)]
        out += [("P", e, 0)] * (n0 + n1 - (n1 % 2))
        if n1 % 2:
            out.append(("P", e, 1))
    for e in sorted({b[1] for b in c if b[0] == "H"}):
        out += [("H", e)] * c[("H", e)]
    for e in sorted({b[1] for b in c if b[0] == "d"}):
        n3 = c[("d", e, 3)]
        n1 = c[("d", e, 1)]
        out += [("d", e, 1)] * (n1 + n3 - (n3 % 2))
        if n3 % 2:
            out.append(("d", e, 3))
    return tuple(out)


def d_of(bl):
    return sum(block_rank(b) for b in bl)


def direct_gram(bl):
    n = d_of(bl)
    A = [[ZERO] * n for _ in range(n)]
    pos = 0
    for b in bl:
        g = block_gram(b)
        for i in range(len(g)):
            for j in range(len(g)):
                A[pos + i][pos + j] = g[i][j]
        pos += len(g)
    return A


def glue_cyclic(bl, x):
    """Blocks of the overlattice S + <x> of S = direct sum of bl (x in S^#, coordinates in the basis of the blocks), and the exponent c of the index of C = <x>.

    格子 S = 直和(bl) に x（S^# の元。座標はブロックの基底について）を足した上格子のブロック列と、C = <x> の指数 c を返す。"""
    A = direct_gram(bl)
    n = len(A)
    c = 0
    i0 = None
    for i, z in enumerate(x):
        if not z.iszero():
            vz = -vpi(z)
            if vz > c:
                c = vz
                i0 = i
    if c == 0:
        return normal_blocks(bl), 0
    # basis {x} + {e_j, j != i0} / 基底 {x} + {e_j, j != i0}
    basis = []
    for j in range(n):
        if j == i0:
            basis.append(list(x))
        else:
            basis.append([ONE if k == j else ZERO for k in range(n)])
    Gm = [[herm(basis[a], A, basis[b]) for b in range(n)] for a in range(n)]
    return normal_blocks(jordan(Gm)), c


# ---------------------------------------------------------------- Statistics of isotropic submodules (forward dynamic programming)
# 等方部分加群の統計（前向きの動的計画）


def compositions_for_block(blocktype, mult):
    """Patterns of element types in mult copies of the same block: yield (multiplicities of the classes, weight, classes).

    同じブロック mult 個の中の元の型のパターン：(類の多重度, 重み, 類) を返す。"""
    cl = list(block_elements(blocktype).items())

    def rec(i, left):
        if i == len(cl) - 1:
            yield [left]
            return
        for k in range(left + 1):
            for rest in rec(i + 1, left - k):
                yield [k] + rest

    for pat in rec(0, mult):
        w = factorial(mult)
        for k, (key, (cnt, rep)) in zip(pat, cl):
            w = w // factorial(k) * cnt ** k
        yield pat, w, cl


def element_patterns(state):
    """Enumerate the patterns of elements y of a state (a list of blocks): yield (sorted state, coordinates, ann, q, weight).

    状態（ブロック列）の元 y のパターンを列挙する：(並べた状態, 全座標, ann, q, 重み) を返す。"""
    groups = Counter(state)
    keys = list(groups.keys())
    # keep the order of state / 順序は state と揃える
    order = []
    for k in keys:
        order.append((k, groups[k]))

    def rec(i):
        if i == len(order):
            yield [], 0, Fr(0), 1
            return
        bt, mult = order[i]
        for pat, w, cl in compositions_for_block(bt, mult):
            coords = []
            ann = 0
            q = Fr(0)
            for k, (key, (cnt, rep)) in zip(pat, cl):
                for _ in range(k):
                    coords += rep
                    ann = max(ann, key[0])
                    q += key[1]
                    if False:
                        pass
                # q has received k times the value above / q には上で k 倍が足された
            for rest, a2, q2, w2 in rec(i + 1):
                yield coords + rest, max(ann, a2), (q + q2) % 2, w * w2

    state_sorted = []
    for bt, mult in order:
        state_sorted += [bt] * mult
    for coords, ann, q, w in rec(0):
        yield state_sorted, coords, ann, q, w


# cache of transitions: (state, block of z, power of pi of z) / 遷移のキャッシュ：(state, z のブロック, z の pi の冪)
_trans = {}


def trans_rank1(state, b, j):
    """Add a block b of rank 1 to D(state), gluing along Z = pi^j D(b) (j = s means Z = 0).
    Returns {new state: {c: weight}}.

    D(state) に階数 1 のブロック b を足す。Z = pi^j D(b) への貼り合わせ（j = s は Z = 0）。
    返り値：{新しい状態: {c: 重み}}。"""
    key = (state, b, j)
    if key in _trans:
        return _trans[key]
    s = block_denom_pi(b)
    out = defaultdict(lambda: defaultdict(int))
    if j == s:
        ns = normal_blocks(state + (b,))
        out[ns][0] += 1
    else:
        z = pi_pow(j - s)
        Ab = block_gram(b)
        qz = herm([z], Ab, [z]).a % 2
        cz = s - j  # exponent of the order of z / z の位数の指数
        total_bl = state + (b,)
        for st_sorted, coords, ann, q, w in element_patterns(state):
            if ann > cz:
                continue
            if (q + qz) % 2 != 0:
                continue
            x = coords + [z]
            ns, c = glue_cyclic(tuple(st_sorted) + (b,), x)
            assert c == cz
            out[ns][c] += w
    res = {k: dict(v) for k, v in out.items()}
    _trans[key] = res
    return res


_rank2_cache = {}


def rank2_submodules(b):
    """Generators of the submodules of D of a block of rank 2 (cyclic or generated by two elements).

    階数 2 のブロックの D の部分加群の生成元の組を列挙する（巡回か 2 元生成）。"""
    if b in _rank2_cache:
        return _rank2_cache[b]
    s = block_denom_pi(b)
    R = residues_pi(s)
    pinv = pi_pow(-s)
    els = [(c1, c2) for c1 in R for c2 in R]

    def tok(c):
        # index in the system of representatives / 代表系での番号
        return c

    # elements as residues: addition and O-multiples are computed in O/pi^s / 加法と O 倍は O/pi^s の演算
    idx = {}
    Rl = [(z.a, z.b) for z in R]

    pinv_s = pi_pow(-s)

    def reduce(z):
        for w in R:
            q = (z - w) * pinv_s
            if q.a.denominator == 1 and q.b.denominator == 1:
                return w
        raise AssertionError

    redcache = {}

    def red(z):
        k = (z.a, z.b)
        if k not in redcache:
            redcache[k] = reduce(z)
        return redcache[k]

    def addv(u, v):
        return (red(u[0] + v[0]), red(u[1] + v[1]))

    def mulv(a, u):
        return (red(a * u[0]), red(a * u[1]))

    def code(u):
        return (u[0].a, u[0].b, u[1].a, u[1].b)

    units = [ONE, G(0, 1), PI]  # generate O (with addition) / 加法と合わせて O を生成

    def span(gens):
        S = {code((G(0), G(0))): (G(0), G(0))}
        frontier = [(G(0), G(0))]
        changed = True
        gens = list(gens)
        while frontier:
            new = []
            for u in frontier:
                for g in gens:
                    w = addv(u, g)
                    if code(w) not in S:
                        S[code(w)] = w
                        new.append(w)
            frontier = new
        return S

    def closure_gens(g):
        out = [g]
        out.append(mulv(G(0, 1), g))
        return out

    # cyclic submodules / 巡回部分加群
    subs = {}
    allv = els
    cyc = {}
    for g in allv:
        gens = closure_gens(g)
        # O g = Z g + Z (i g)
        S = span(gens)
        k = frozenset(S.keys())
        if k not in cyc:
            cyc[k] = g
    for k, g in cyc.items():
        subs[k] = (g,)
    # submodules generated by two elements / 2 元生成
    cl = list(cyc.items())
    for (k1, g1), (k2, g2) in itertools.combinations(cl, 2):
        S = span(closure_gens(g1) + closure_gens(g2))
        k = frozenset(S.keys())
        if k not in subs:
            subs[k] = (g1, g2)
    res = [tuple([c * pinv for c in g] for g in gens) for gens in subs.values()]
    _rank2_cache[b] = res
    return res


def all_elements(state):
    """All elements of a state (coordinates); for small states only.

    状態の元を全部（座標）で列挙する。小さい状態専用。"""
    lists = []
    for b in state:
        s = block_denom_pi(b)
        R = residues_pi(s)
        pinv = pi_pow(-s)
        if block_rank(b) == 1:
            lists.append([[c * pinv] for c in R])
        else:
            lists.append([[c1 * pinv, c2 * pinv] for c1 in R for c2 in R])
    for combo in itertools.product(*lists):
        x = []
        for c in combo:
            x += c
        yield x




def glue_multi(bl, xs):
    """Add several x one after the other (one overlattice at a time). Returns (blocks, total exponent of the index).

    複数の x を順に足す（上格子を 1 個ずつ）。返り値は (ブロック列, 指数の合計)。"""
    A = direct_gram(bl)
    n = len(A)
    basis = [[ONE if k == j else ZERO for k in range(n)] for j in range(n)]

    def coords_in_basis(x, basis):
        # coordinates through the inverse of basis / basis の逆行列で座標を求める
        M = [row[:] + [ONE if i == j else ZERO for j in range(n)] for i, row in enumerate(basis)]
        for col in range(n):
            piv = min((i for i in range(col, n) if not M[i][col].iszero()), key=lambda i: vpi(M[i][col]))
            M[col], M[piv] = M[piv], M[col]
            inv = M[col][col].inv()
            M[col] = [e * inv for e in M[col]]
            for i in range(n):
                if i != col and not M[i][col].iszero():
                    f = M[i][col]
                    M[i] = [e - f * m for e, m in zip(M[i], M[col])]
        Binv = [row[n:] for row in M]
        out = [ZERO] * n
        for i in range(n):
            if x[i].iszero():
                continue
            for j in range(n):
                out[j] = out[j] + x[i] * Binv[i][j]
        return out

    total = 0
    for x in xs:
        co = coords_in_basis(x, basis)
        c = 0
        i0 = None
        for i, z in enumerate(co):
            if not z.iszero():
                vz = -vpi(z)
                if vz > c:
                    c = vz
                    i0 = i
        if c == 0:
            continue
        basis[i0] = list(x)
        total += c
    Gm = [[herm(basis[a], A, basis[b]) for b in range(n)] for a in range(n)]
    return normal_blocks(jordan(Gm)), total


_trans2 = {}


_pair_cache = {}


def block_pair_classes(tau):
    """Split the pairs (y1, y2) of elements of the block tau into orbits of the diagonal action of (a subgroup of) Aut(tau).
    Returns [(size, y1, y2, ann1, ann2, q1, q2, real part of b12, imaginary part)], b12 = h(y1, y2) mod O.

    ブロック tau の元の組 (y1, y2) を、tau の自己同型（の部分群）の対角の作用の軌道に分ける。
    返り値：[(個数, y1, y2, ann1, ann2, q1, q2, b12 の実部, 虚部)]。b12 は h(y1, y2) mod O。"""
    if tau in _pair_cache:
        return _pair_cache[tau]
    A = block_gram(tau)
    els, perms = block_autos(tau)
    N = len(els)
    perms = perms[:200]
    parent = list(range(N * N))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for perm in perms:
        for i in range(N):
            pi_ = perm[i]
            for j in range(N):
                ra, rb = find(i * N + j), find(pi_ * N + perm[j])
                if ra != rb:
                    parent[ra] = rb
    cnt = Counter(find(k) for k in range(N * N))
    out = []
    ann = [ann_of(x) for x in els]
    q = [herm(x, A, x).a % 2 for x in els]
    for root, c in cnt.items():
        i, j = root // N, root % N
        bb = herm(els[i], A, els[j])
        out.append((c, els[i], els[j], ann[i], ann[j], q[i], q[j], bb.a % 1, bb.b % 1))
    _pair_cache[tau] = out
    return out


def trans_rank2(state, b):
    """Add a block b of rank 2 to D(state). y runs over patterns; (y1, y2) over patterns of orbits of pairs in each block.

    D(state) に階数 2 のブロック b を足す。y はパターン、(y1, y2) はブロックごとの組の軌道のパターンで数える。"""
    key = (state, b)
    if key in _trans2:
        return _trans2[key]
    out = defaultdict(lambda: defaultdict(int))
    Ab = block_gram(b)
    s_b = block_denom_pi(b)
    R_b = _res(s_b)
    pats = list(element_patterns(state))
    st_sorted = tuple(pats[0][0]) if pats else ()
    groups = list(Counter(st_sorted).items())
    # st_sorted is grouped by block type / st_sorted の順はブロックの型ごとにまとまっている
    for gens in rank2_submodules(b):
        if len(gens) == 1:
            z = gens[0]
            cz = ann_of(z)
            if cz == 0:
                ns = normal_blocks(state + (b,))
                out[ns][0] += 1
                continue
            qz = herm(z, Ab, z).a % 2
            for _, y, an, q, w in pats:
                if an > cz or (q + qz) % 2 != 0:
                    continue
                ns, c = glue_multi(st_sorted + (b,), [y + z])
                out[ns][c] += w
            continue
        z1, z2 = gens
        c1, c2 = ann_of(z1), ann_of(z2)
        q1 = herm(z1, Ab, z1).a % 2
        q2 = herm(z2, Ab, z2).a % 2
        b12 = herm(z1, Ab, z2)
        # relations a1 z1 + a2 z2 = 0 between the generators of Z / Z の生成元の間の関係
        rel = []
        for a1 in R_b:
            for a2 in R_b:
                w = [a1 * u + a2 * v for u, v in zip(z1, z2)]
                if all(t.a.denominator == 1 and t.b.denominator == 1 for t in w):
                    rel.append((a1, a2))

        def allowed(cl):
            _, y1, y2, an1, an2 = cl[:5]
            if an1 > c1 or an2 > c2:
                return False
            for a1, a2 in rel:
                w = [a1 * u + a2 * v for u, v in zip(y1, y2)]
                if not all(t.a.denominator == 1 and t.b.denominator == 1 for t in w):
                    return False
            return True

        lists = [(tau, m, [cl for cl in block_pair_classes(tau) if allowed(cl)]) for tau, m in groups]

        def rec(i):
            if i == len(lists):
                yield [], [], Fr(0), Fr(0), Fr(0), Fr(0), 1
                return
            tau, m, cls = lists[i]
            if not cls and m > 0:
                return

            def comp(k, left):
                if k == len(cls) - 1:
                    yield [left]
                    return
                for t in range(left + 1):
                    for rest in comp(k + 1, left - t):
                        yield [t] + rest

            for pat in comp(0, m):
                w = factorial(m)
                ya, yb = [], []
                sq1 = sq2 = sre = sim = Fr(0)
                for t, cl in zip(pat, cls):
                    if t == 0:
                        continue
                    w = w // factorial(t) * cl[0] ** t
                    for _ in range(t):
                        ya += cl[1]
                        yb += cl[2]
                    sq1 += t * cl[5]
                    sq2 += t * cl[6]
                    sre += t * cl[7]
                    sim += t * cl[8]
                for ra, rb, q1r, q2r, rer, imr, wr in rec(i + 1):
                    yield ya + ra, yb + rb, sq1 + q1r, sq2 + q2r, sre + rer, sim + imr, w * wr

        for ya, yb, sq1, sq2, sre, sim, w in rec(0):
            if (sq1 + q1) % 2 != 0 or (sq2 + q2) % 2 != 0:
                continue
            if ((sre + b12.a) % 1) != 0 or ((sim + b12.b) % 1) != 0:
                continue
            ns, c = glue_multi(st_sorted + (b,), [ya + z1, yb + z2])
            out[ns][c] += w
    res = {k: dict(v) for k, v in out.items()}
    _trans2[key] = res
    return res






# ---- closed formula for states made of P blocks only (keeps the number of patterns small)
# ---- P だけの状態の閉じた式（元の列挙のパターンの個数を抑える）


def gauss_binom(n, k, q):
    if k < 0 or k > n:
        return 0
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (i + 1) - 1
    return num // den


def singular_subspaces(m, eps, k):
    """Number of totally singular k-dimensional subspaces of O^eps_{2m}(2), eps = +1, -1.

    O^eps_{2m}(2) の全特異 k 次元部分空間の個数。eps = +1, -1。"""
    if k > m or (eps == -1 and k == m):
        return 0
    out = 1
    for i in range(k):
        # number of nonzero singular vectors: (2^{m-i} - eps)(2^{m-i-1} + eps) / 特異な非零ベクトルの数
        out *= (2 ** (m - i) - eps) * (2 ** (m - i - 1) + eps)
    den = 1
    for i in range(k):
        den *= 2 ** (i + 1) - 1
    assert out % den == 0
    return out // den


_tick = lambda: None


def stats_blocks(blocks):
    """From a multiset of 2-adic blocks, compute {(ell, state): weight} by forward dynamic programming.
    Order: P, H, then the blocks of rank 1.

    2 進ブロックの列（多重集合）から {(ell, 状態): 重み} を前向きの動的計画で求める。
    順序：P、H、次に階数 1。"""
    blocks = list(blocks)
    Ps = [b for b in blocks if b[0] == "P"]
    Hs = [b for b in blocks if b[0] == "H"]
    Ds = [b for b in blocks if b[0] == "d"]
    # initial state: closed formula for the P blocks / 初期状態：P だけなら閉じた式
    dist = defaultdict(int)
    assert all(b[1] == 0 for b in Ps), "only P with e = 0 occurs in root lattices"
    m = len(Ps)
    n1 = sum(1 for b in Ps if b[2] == 1)
    eps = -1 if n1 % 2 else 1
    for k in range(m + 1):
        cnt = singular_subspaces(m, eps, k)
        if cnt == 0:
            continue
        mm = m - k
        st = []
        if mm > 0:
            st = [("P", 0, 0)] * (mm - 1) + ([("P", 0, 1)] if eps == -1 else [("P", 0, 0)])
        dist[(k, normal_blocks(tuple(st)))] += cnt
    for b in Hs:
        new = defaultdict(int)
        for (ell, st), w in dist.items():
            for ns, d in trans_rank2(st, b).items():
                for c, w2 in d.items():
                    new[(ell + c, ns)] += w * w2
            _tick()
        dist = new
    for b in Ds:
        new = defaultdict(int)
        s = block_denom_pi(b)
        for (ell, st), w in dist.items():
            for j in range(s + 1):
                for ns, d in trans_rank1(st, b, j).items():
                    for c, w2 in d.items():
                        new[(ell + c, ns)] += w * w2
            _tick()
        dist = new
    return dist


_stats_store = {}


def stats_cached(blocks):
    if blocks in _stats_store:
        return _stats_store[blocks]
    dist = stats_blocks(blocks)
    out = defaultdict(int)
    for (ell, st), w in dist.items():
        out[(ell, d_of(st))] += w
    res = dict(out)
    _stats_store[blocks] = res
    return res


# naive check: build all of D and enumerate the isotropic submodules / 素朴な検算：D を全部作って等方部分加群を列挙する


def brute_stats(blocks):
    """{(ell, d): count} by brute force; for small D only. Elements are tuples of Fraction coordinates.

    {(ell, d): 個数} を全列挙で。小さい D 専用。元を座標の Fraction の組で持つ。"""
    blocks = list(blocks)
    A = direct_gram(blocks)
    n = len(A)
    coords = []
    for b in blocks:
        s = block_denom_pi(b)
        coords.append((block_rank(b), s))
    els = [tuple(x) for x in all_elements(blocks)]

    def key(x):
        # index of the representative modulo the lattice / 格子を法とする代表の番号
        out = []
        pos = 0
        for (rk, s) in coords:
            for i in range(rk):
                z = x[pos]
                out.append(z)
                pos += 1
        return out

    # reduce to integer coordinates: x_i = pi^{-s} c_i, c_i a representative in O/pi^s / 整数の座標に落とす
    resid = {}

    def reduce_vec(x):
        out = []
        pos = 0
        for (rk, s) in coords:
            pis = pi_pow(s)
            pinv = pi_pow(-s)
            R = _res(s)
            for i in range(rk):
                c = x[pos] * pis
                for j, w in enumerate(R):
                    q = (c - w) * pinv
                    if q.a.denominator == 1 and q.b.denominator == 1:
                        out.append(j)
                        break
                else:
                    raise AssertionError
                pos += 1
        return tuple(out)

    index = {reduce_vec(x): i for i, x in enumerate(els)}
    N = len(els)
    add_cache = {}

    def add(i, j):
        k = (i, j) if i <= j else (j, i)
        if k not in add_cache:
            add_cache[k] = index[reduce_vec([u + v for u, v in zip(els[i], els[j])])]
        return add_cache[k]

    def mul(g, i):
        return index[reduce_vec([g * u for u in els[i]])]

    Ic = [mul(G(0, 1), i) for i in range(N)]
    Pc = [mul(PI, i) for i in range(N)]
    q = [herm(list(x), A, list(x)).a % 2 for x in els]
    zero = index[reduce_vec([ZERO] * n)]

    def bil(i, j):
        z = herm(list(els[i]), A, list(els[j]))
        return z.a.denominator == 1 and z.b.denominator == 1

    good = [i for i in range(N) if q[i] == 0]
    start = frozenset([zero])
    seen = {start}
    stack = [start]
    res = Counter()
    while stack:
        C = stack.pop()
        Cl = list(C)
        perp = [i for i in range(N) if all(bil(i, c) for c in Cl)]
        pset = set(perp)
        # D' = perp / C、d = dim D'[pi]
        tor = [i for i in perp if Pc[i] in C]
        dd = 0
        cnt = len(tor) // len(C)
        while 2 ** dd < cnt:
            dd += 1
        ell = 0
        while 2 ** ell < len(C):
            ell += 1
        res[(ell, dd)] += 1
        for x in perp:
            if x in C or q[x] != 0:
                continue
            # C + O x
            newset = set(C)
            frontier = [x]
            gens = [x, Ic[x]]
            cur = set(C)
            changed = True
            while changed:
                changed = False
                for g in gens:
                    add_list = [add(c, g) for c in cur]
                    for a in add_list:
                        if a not in cur:
                            cur.add(a)
                            changed = True
            nc = frozenset(cur)
            if nc not in seen and all(q[c] == 0 for c in nc):
                seen.add(nc)
                stack.append(nc)
    return dict(res)


_res_cache = {}


def _res(s):
    if s not in _res_cache:
        _res_cache[s] = residues_pi(s)
    return _res_cache[s]


# ---------------------------------------------------------------- Types of root lattices / ルート格子の型

# Names of irreducible components: ('A', n), ('D', n) (n >= 4), ('E', n), ('G', m) (m >= 2, D_{2m}^g), ('Eg', 8) (E_8^g, rank 4)
# 既約成分の名前：上の 5 種（G は D_{2m}^g、Eg は階数 4 の E_8^g）


def rank_of(Y):
    if Y[0] == "Eg":
        return 4
    return Y[1]


def irreducibles(rmax):
    out = []
    for n in range(1, rmax + 1):
        out.append(("A", n))
    for n in range(4, rmax + 1):
        out.append(("D", n))
    for n in (6, 7, 8):
        if n <= rmax:
            out.append(("E", n))
    for m in range(2, rmax + 1):
        out.append(("G", m))
    if rmax >= 4:
        out.append(("Eg", 8))
    return out


def U_order(Y):
    """Order of the automorphism group as a Z[i]-lattice.

    Z[i] 格子としての自己同型群の位数。"""
    t, n = Y
    if t == "A":
        if n == 1:
            return 4
        return 2 * 2 * factorial(n + 1)
    if t == "D":
        if n == 4:
            return 2 * 1152
        return 2 * 2 ** n * factorial(n)
    if t == "E":
        return 2 * {6: 2 * 51840, 7: 2903040, 8: 696729600}[n]
    if t == "G":
        if n == 2:
            return 96
        return 4 ** n * factorial(n)
    if t == "Eg":
        return 46080
    raise ValueError(Y)


def det_of(Y):
    """Determinant over Z[i] (a positive rational) = |D|^{1/2}.

    Z[i] 上の行列式（正の有理数）= |D|^{1/2}。"""
    t, n = Y
    if t == "A":
        return n + 1
    if t == "D":
        return 4
    if t == "E":
        return {6: 3, 7: 2, 8: 1}[n]
    if t == "G":
        return 2
    return 1


def blocks2(Y):
    """Blocks of the 2-adic non-unimodular part of Y.

    Y の 2 進の非ユニモジュラー部分のブロック。"""
    t, n = Y
    if t == "A":
        m = n + 1
        e = 0
        while m % 2 == 0:
            m //= 2
            e += 1
        if e == 0:
            return []
        # generator y_2 = m y of D = O/2^e, q = m^2 n/(n+1) = n m / 2^e; compare with q(y) = u/2^e of <2^e u>
        # D = O/2^e の生成元 y_2 = m y、q = n m / 2^e。<2^e u> の q(y) = u/2^e と比べる
        u = (n * m) % 4
        return [("d", e, u)]
    if t == "D":
        if n % 2:
            return [("d", 2, n % 4)]
        if n % 4 == 0:
            return [("H", 1)]
        return [("d", 1, (n // 2) % 4), ("d", 1, (n // 2) % 4)]
    if t == "E":
        return {6: [], 7: [("d", 1, 3)], 8: []}[n]
    if t == "G":
        if n % 2:
            return [("d", 1, n % 4)]
        return [("P", 0, 1 if n % 4 == 2 else 0)]
    return []


def oddprimes_of(Y):
    """{p: list of e}: the odd-prime part of A_n and E_6.

    {p: e のリスト}：A_n、E_6 の奇素数部分。"""
    t, n = Y
    out = defaultdict(list)
    if t == "A":
        m = n + 1
        p = 3
        while m > 1 and p <= m:
            if m % p == 0:
                e = 0
                while m % p == 0:
                    m //= p
                    e += 1
                out[p].append(e)
            p += 2
    elif t == "E" and n == 6:
        out[3].append(1)
    return out


# ---------------------------------------------------------------- Local factors at odd primes / 奇素数の局所因子


def chi(p):
    return 1 if p % 4 == 1 else -1


def f_odd(p, a, b, t):
    """f_p for D_p = (O/p)^a + (O/p^2)^b; b is 0, 1 or 2 (b = 2 only for p = 3 from A_8).

    D_p = (O/p)^a 直和 (O/p^2)^b の f_p。b は 0, 1, 2（2 は p = 3 の A_8 のみ）。"""
    def lf(d):
        out = Fr(1)
        for j in range(t - d + 1, t + 1):
            s = chi(p) ** j if j >= 0 else 1
            if j == 0:
                return Fr(0)
            out *= 1 - Fr(s, p ** j) if j > 0 else Fr(1)
        return out

    def term(size_exp, d):
        return Fr(1, p ** (size_exp * t)) * lf(d)

    tot = Fr(0)
    if chi(p) == 1:
        assert b == 0
        for k1 in range(a + 1):
            for k2 in range(a - k1 + 1):
                cnt = gauss_binom(a, k1, p) * gauss_binom(a - k1, k2, p)
                tot += cnt * term(k1 + k2, a - k1 - k2)
        return tot

    def N(a, k):
        """Number of totally isotropic k-dimensional subspaces of a unitary space of dimension a over F_{p^2}.

        F_{p^2} のユニタリ空間（次元 a）の全等方 k 次元部分空間の個数。"""
        if 2 * k > a:
            return 0
        num = 1
        den = 1
        for i in range(k):
            s = a - 2 * i
            num *= (p ** s - (-1) ** s) * (p ** (s - 1) - (-1) ** (s - 1))
            den *= p ** (2 * (i + 1)) - 1
        assert num % den == 0
        return num // den

    def iso_vec(a):
        return (p ** a - (-1) ** a) * (p ** (a - 1) - (-1) ** (a - 1)) if a >= 1 else 0

    if b == 0:
        for k in range(a // 2 + 1):
            tot += N(a, k) * term(2 * k, a - 2 * k)
        return tot
    if b == 1:
        # |C| = p^{2k} * (glue part): Z = 0, or Z = p D(b) matched with an isotropic w (w = 0 included)
        # |C| = p^{2k} *（貼り合わせの部分）。Z = 0、または Z = p D(b) に等方な w を対応させる（w = 0 も含む）
        for k in range(a // 2 + 1):
            ap = a - 2 * k
            c0 = N(a, k)
            # Z = 0: D' = F^{ap} + O/p^2, d = ap + 1
            tot += c0 * term(2 * k, ap + 1)
            # w = 0: add p y0 to C; D' = F^{ap}, d = ap / C に p y0 を足す
            tot += c0 * term(2 * k + 2, ap)
            # w != 0 isotropic: D' = O/p^2 + F^{ap-2}, d = ap - 1 / w は 0 でない等方
            if ap >= 2:
                tot += c0 * iso_vec(ap) * term(2 * k + 2, ap - 1)
        return tot
    if b == 2:
        assert a == 0 and p == 3
        return f_odd_brute_9x9(t)
    raise ValueError


def f_odd_brute_9x9(t):
    """f_3 for D = (O/9)^2 (two copies of <u/9>, p = 3) by full enumeration; b(y, y) = 1/9 (unit differences are absorbed by norms).

    D = (O/9)^2（<u/9> 2 個、p = 3）の f_3 を全列挙で。b(y, y) = 1/9（単数の違いはノルムで吸収される）。"""
    # elements of O/9 as pairs (a, b) mod 9 / O/9 の元は (a, b) mod 9
    p = 3
    mod = 9
    els = [(a, b, c, d) for a in range(9) for b in range(9) for c in range(9) for d in range(9)]

    def hval(x, y):
        # h(x, y) = (x1 conj(y1) + x2 conj(y2)) / 9 mod O; returns the numerator / 分子を返す
        re = x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3]
        im = x[1] * y[0] - x[0] * y[1] + x[3] * y[2] - x[2] * y[3]
        return re % 9, im % 9

    def add(x, y):
        return tuple((u + v) % 9 for u, v in zip(x, y))

    def mi(x):
        return ((-x[1]) % 9, x[0], (-x[3]) % 9, x[2])

    def scal(k, x):
        return tuple(k * u % 9 for u in x)

    isotropic = [x for x in els if hval(x, x) == (0, 0)]
    # enumerate the submodules / 部分加群の列挙
    zero = (0, 0, 0, 0)
    start = frozenset([zero])
    seen = {start}
    stack = [start]
    res = Fr(0)

    def lf(d):
        out = Fr(1)
        for j in range(t - d + 1, t + 1):
            if j == 0:
                return Fr(0)
            if j > 0:
                out *= 1 - Fr((-1) ** j, p ** j)
        return out

    while stack:
        C = stack.pop()
        perp = [x for x in els if all(hval(x, c) == (0, 0) for c in C)]
        # d = dim (perp/C)[3]
        tor = [x for x in perp if scal(3, x) in C]
        cnt = len(tor) // len(C)
        d = 0
        while 9 ** d < cnt:
            d += 1
        res += Fr(1, len(C) ** t) * lf(d)
        for x in perp:
            if x in C or hval(x, x) != (0, 0):
                continue
            cur = set(C)
            gens = [x, mi(x)]
            changed = True
            while changed:
                changed = False
                for g in gens:
                    for c in list(cur):
                        w = add(c, g)
                        if w not in cur:
                            cur.add(w)
                            changed = True
            nc = frozenset(cur)
            if nc not in seen:
                seen.add(nc)
                stack.append(nc)
    return res


def f_two(blocks, t):
    st = stats_cached(tuple(sorted(blocks)))
    tot = Fr(0)
    for (ell, d), cnt in st.items():
        v = Fr(cnt, 2 ** (ell * t))
        zero = False
        for j in range(t - d + 2, t + 1):
            if j % 2 == 0:
                if j == 0:
                    zero = True
                    break
                if j > 0:
                    v *= 1 - Fr(1, 2 ** j)
        if not zero:
            tot += v
    return tot


def A_value(T, r):
    """A_r(T); T is {irreducible: multiplicity}.

    A_r(T)。T は {既約: 多重度}。"""
    n = sum(rank_of(Y) * m for Y, m in T.items())
    t = r - n
    if t < 0:
        return Fr(0)
    det = 1
    bl = []
    odd = defaultdict(lambda: [0, 0])
    for Y, m in T.items():
        det *= det_of(Y) ** m
        bl += blocks2(Y) * m
        for p, es in oddprimes_of(Y).items():
            for e in es:
                if e == 1:
                    odd[p][0] += m
                elif e == 2:
                    odd[p][1] += m
                else:
                    raise ValueError
    val = Sprime(t) * det ** t
    if t == 0 and False:
        pass
    for p, (a, b) in odd.items():
        val *= f_odd(p, a, b, t)
        if val == 0:
            return val
    val *= f_two(bl, t)
    return val


# ---------------------------------------------------------------- Moebius function (root sublattices) / Moebius 関数（部分ルート格子）

_c_table = {}


def _egf_coeff_vals(rmax):
    """c(A_n), c(D_n), c(G_m) by recursions of exponential generating functions.

    c(A_n)、c(D_n)、c(G_m) を指数型母関数の漸化で求める。"""
    cA = {0: Fr(1)}
    for n in range(1, rmax + 1):
        cA[n] = Fr((-1) ** n * factorial(n))

    # root sublattices of D_n: set partitions of the vertices; a block of size v >= 2 carries A_{v-1} (2^{v-1} ways) or D_v, singletons are unused
    # D_n の部分ルート格子：頂点の集合分割で、大きさ v >= 2 のブロックに A_{v-1}（2^{v-1} 通り）か D_v、単点は未使用
    # D_2 = A_1^2, D_3 = A_3
    cD = {2: Fr(1), 3: cA[3]}

    def exp_coeffs(block, N):
        """Coefficients of x^m/m! (m <= N) of exp(sum block[v] x^v / v!) (block is a dict on 1..N).

        exp(sum block[v] x^v / v!) の x^m/m! の係数（block は 1..N の辞書）。"""
        e = [Fr(0)] * (N + 1)
        e[0] = Fr(1)
        for m in range(1, N + 1):
            s = Fr(0)
            for v in range(1, m + 1):
                if v in block:
                    s += comb(m - 1, v - 1) * block[v] * e[m - v]
            e[m] = s
        return e

    for n in range(4, rmax + 1):
        block = {1: Fr(1)}
        for v in range(2, n):
            block[v] = 2 ** (v - 1) * cA[v - 1] + cD[v]
        block[n] = 2 ** (n - 1) * cA[n - 1]
        e = exp_coeffs(block, n)
        cD[n] = -e[n]
    # consistency check for D_2, D_3 (the sum is 0) / D_2、D_3 の整合性（総和が 0）の検査
    for n in (2, 3):
        block = {1: Fr(1)}
        for v in range(2, n + 1):
            block[v] = 2 ** (v - 1) * cA[v - 1] + cD[v]
        e = exp_coeffs(block, n)
        assert e[n] == 0, (n, e[n])
    cG = {1: Fr(-1)}
    for m in range(2, rmax + 1):
        block = {1: Fr(0)}
        for v in range(2, m):
            block[v] = 4 ** (v - 1) * cA[v - 1] + 2 ** (v - 1) * cD[v] + cG[v]
        block[m] = 4 ** (m - 1) * cA[m - 1] + 2 ** (m - 1) * cD[m]
        e = exp_coeffs(block, m)
        cG[m] = -e[m]
    return cA, cD, cG


# ---- exceptional types (count embeddings with explicit root systems) / 例外型（ルート系を明示して埋め込みを数える）


def e8_roots():
    """The roots of E_8 in integer coordinates (scaled by 2).

    E_8 のルートを整数座標（2 倍）で。"""
    roots = set()
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (2, -2):
                for sj in (2, -2):
                    v = [0] * 8
                    v[i] = si
                    v[j] = sj
                    roots.add(tuple(v))
    for signs in itertools.product((1, -1), repeat=8):
        if sum(1 for s in signs if s < 0) % 2 == 0:
            roots.add(tuple(signs))
    return sorted(roots)


def ip(u, v):
    return sum(a * b for a, b in zip(u, v)) // 4


def e_roots(n):
    R = e8_roots()
    if n == 8:
        return R
    # E_7: orthogonal to a fixed root a; E_6: orthogonal to an A_2 / E_7 は a に直交、E_6 は A_2 に直交
    a = R[0]
    E7 = [x for x in R if ip(x, a) == 0]
    if n == 7:
        return E7
    b = next(x for x in R if ip(x, a) == -1)
    E6 = [x for x in E7 if ip(x, b) == 0]
    return E6


def cartan_basis(Y):
    """Gram matrix of the simple roots of an irreducible ADE root lattice over Z.

    既約の Z 型 ADE の単純ルートの Gram 行列。"""
    t, n = Y
    G_ = [[0] * n for _ in range(n)]
    for i in range(n):
        G_[i][i] = 2
    if t == "A":
        edges = [(i, i + 1) for i in range(n - 1)]
    elif t == "D":
        edges = [(i, i + 1) for i in range(n - 2)] + [(n - 3, n - 1)]
    elif t == "E":
        edges = [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]
    for i, j in edges:
        G_[i][j] = G_[j][i] = -1
    return G_


def count_embeddings_real(roots, gram):
    """Number of ordered sequences of roots with Gram matrix gram, counted by orbits under reflections in orthogonal roots.

    ルートの列の Gram が gram になる順序付き列の個数。直交するルートの鏡映で軌道に分けて数える。"""
    import numpy as np

    M = np.array(roots, dtype=np.int64)
    IP = (M @ M.T) // 4
    N = len(roots)
    neg = {}
    index = {tuple(r): i for i, r in enumerate(roots)}
    # permutations of the reflections, computed lazily / 鏡映の置換を遅延計算
    refl = {}

    def reflection(i):
        if i not in refl:
            v = M - (IP[:, i][:, None] * M[i][None, :])
            refl[i] = np.array([index[tuple(x)] for x in v.tolist()], dtype=np.int32)
        return refl[i]

    n = len(gram)
    memo = {}

    def rec(prefix, level):
        if level == n:
            return 1
        cand = np.ones(N, dtype=bool)
        for k, p in enumerate(prefix):
            cand &= IP[p] == gram[k][level]
        cand_idx = np.nonzero(cand)[0]
        if len(cand_idx) == 0:
            return 0
        # roots orthogonal to the prefix / prefix に直交するルート
        orth = np.ones(N, dtype=bool)
        for p in prefix:
            orth &= IP[p] == 0
        orth_idx = np.nonzero(orth)[0]
        key = (level, tuple(sorted(orth_idx.tolist()))[:0], cand_idx.tobytes())
        # orbit decomposition / 軌道分解
        seen = set()
        total = 0
        gens = [reflection(i) for i in orth_idx.tolist()]
        for c in cand_idx.tolist():
            if c in seen:
                continue
            orb = {c}
            st = [c]
            while st:
                x = st.pop()
                for g in gens:
                    y = int(g[x])
                    if y not in orb:
                        orb.add(y)
                        st.append(y)
            seen |= orb
            total += len(orb) * rec(prefix + [c], level + 1)
        return total

    return rec([], 0)


def ade_types_upto(n, allowed=None):
    """Multisets of ADE components of rank <= n (dicts {('A', k): m}).

    階数 <= n の ADE の多重集合（辞書 {('A', k): m}）。"""
    comps = []
    for k in range(1, n + 1):
        comps.append(("A", k))
    for k in range(4, n + 1):
        comps.append(("D", k))
    for k in (6, 7, 8):
        if k <= n:
            comps.append(("E", k))

    def rec(i, left):
        if i == len(comps):
            yield {}
            return
        Y = comps[i]
        for m in range(left // Y[1] + 1):
            for rest in rec(i + 1, left - m * Y[1]):
                d = dict(rest)
                if m:
                    d[Y] = m
                yield d

    return list(rec(0, n))


def aut_real(Y):
    t, n = Y
    if t == "A":
        return 2 if n == 1 else 2 * factorial(n + 1)
    if t == "D":
        return 1152 if n == 4 else 2 ** n * factorial(n)
    return {6: 2 * 51840, 7: 2903040, 8: 696729600}[n]


def block_gram_real(T):
    comps = []
    for Y, m in sorted(T.items()):
        comps += [Y] * m
    n = sum(Y[1] for Y in comps)
    G_ = [[0] * n for _ in range(n)]
    pos = 0
    for Y in comps:
        c = cartan_basis(Y)
        for i in range(Y[1]):
            for j in range(Y[1]):
                G_[pos + i][pos + j] = c[i][j]
        pos += Y[1]
    return G_


@lru_cache(maxsize=None)
def c_exceptional(n):
    """c(E_n[i]) = -sum_{S subset E_n} c(S)。"""
    cA, cD, cG = _egf_coeff_vals(8)
    roots = e_roots(n)
    assert len(roots) == {6: 72, 7: 126, 8: 240}[n]
    total = Fr(0)
    detail = {}
    for T in ade_types_upto(n):
        if T == {("E", n): 1}:
            continue
        cs = Fr(1)
        aut = 1
        for Y, m in T.items():
            if Y[0] == "A":
                cy = cA[Y[1]]
            elif Y[0] == "D":
                cy = cD[Y[1]]
            else:
                cy = c_exceptional(Y[1])
            cs *= cy ** m
            aut *= aut_real(Y) ** m * factorial(m)
        if T:
            emb = count_embeddings_real(roots, block_gram_real(T))
        else:
            emb = 1
        assert emb % aut == 0
        total += cs * (emb // aut)
        detail[tuple(sorted(T.items()))] = emb // aut
    return -total


# ---- E_8^g (G_31): 60 lines / E_8^g（G_31）：60 本の直線




def build_eg():
    """E_8^g as a Z[i]-lattice: the roots of E_8 with g acting on the coordinates.
    On the roots of E_8, pairing coordinates (x1, x2) -> (-x2, x1) is an element of W(D_8) and preserves E_8.

    Z[i] 格子 E_8^g を、E_8 の座標に g を入れて作る。
    E_8 のルートに g を与える：座標を 2 個ずつ組にした (x1, x2) -> (-x2, x1) は W(D_8) の元で、E_8 を保つ。"""
    R = e8_roots()

    def g(v):
        return (-v[1], v[0], -v[3], v[2], -v[5], v[4], -v[7], v[6])

    Rs = set(R)
    assert all(g(v) in Rs for v in R)
    return R, g




def complex_embed_count(R, g, target):
    """Number of sequences of roots whose Gram matrix over Z[i] is target (a matrix of (re, im)); orbits under reflections.

    Z[i] 上の Gram が target（(re, im) の行列）となるルートの列の個数。鏡映で軌道に分ける。"""
    import numpy as np

    M = np.array(R, dtype=np.int64)
    GM = np.array([g(v) for v in R], dtype=np.int64)
    Re = (M @ M.T) // 4
    Im = (M @ GM.T) // 4
    N = len(R)
    index = {tuple(r): i for i, r in enumerate(R)}
    refl = {}

    def reflection(i):
        # s_x(z) = z - h(z, x) x (i acts as g) / i は g として働く
        if i not in refl:
            xv = M[i]
            gx = GM[i]
            # z - Re(h) x - Im(h) (i x) = z - (z,x) x - (z, g x) g x
            v = M - Re[:, i][:, None] * xv[None, :] - Im[:, i][:, None] * gx[None, :]
            refl[i] = np.array([index[tuple(y)] for y in v.tolist()], dtype=np.int32)
        return refl[i]

    n = len(target)

    def rec(prefix, level):
        if level == n:
            return 1
        cand = np.ones(N, dtype=bool)
        for k, p in enumerate(prefix):
            # h(x_k, x_level) = target[k][level]
            cand &= (Re[p] == target[k][level][0]) & (Im[p] == target[k][level][1])
        cand_idx = np.nonzero(cand)[0]
        if len(cand_idx) == 0:
            return 0
        orth = np.ones(N, dtype=bool)
        for p in prefix:
            orth &= (Re[p] == 0) & (Im[p] == 0)
        orth_idx = np.nonzero(orth)[0]
        gens = [reflection(i) for i in orth_idx.tolist()]
        seen = set()
        total = 0
        for c in cand_idx.tolist():
            if c in seen:
                continue
            orb = {c}
            st = [c]
            while st:
                x = st.pop()
                for gg in gens:
                    y = int(gg[x])
                    if y not in orb:
                        orb.add(y)
                        st.append(y)
            seen |= orb
            total += len(orb) * rec(prefix + [c], level + 1)
        return total

    return rec([], 0)


def complex_gram(T):
    """Gram matrix ((re, im)) of a Z[i]-basis of roots of T (a multiset of Z[i] root lattices).

    T（Z[i] ルート格子の多重集合）のルートからなる Z[i] 基底の Gram（(re, im)）。"""
    comps = []
    for Y, m in sorted(T.items()):
        comps += [Y] * m
    n = sum(rank_of(Y) for Y in comps)
    M = [[(0, 0)] * n for _ in range(n)]
    pos = 0
    for Y in comps:
        k = rank_of(Y)
        c = irr_gram_complex(Y)
        for i in range(k):
            for j in range(k):
                M[pos + i][pos + j] = c[i][j]
        pos += k
    return M


def irr_gram_complex(Y):
    t, n = Y
    if t in ("A", "D", "E"):
        c = cartan_basis(Y)
        return [[(c[i][j], 0) for j in range(n)] for i in range(n)]
    if t == "G":
        # basis pi e_1, e_1 - e_2, e_2 - e_3, ...; h(x, y) = sum x_p conj(y_p) / 基底と h
        vecs = []
        v = [(0, 0)] * n
        v[0] = (1, 1)
        vecs.append(v)
        for k in range(n - 1):
            w = [(0, 0)] * n
            w[k] = (1, 0)
            w[k + 1] = (-1, 0)
            vecs.append(w)

        def h(x, y):
            a = b = 0
            for (x1, x2), (y1, y2) in zip(x, y):
                a += x1 * y1 + x2 * y2
                b += x2 * y1 - x1 * y2
            return (a, b)

        return [[h(vecs[i], vecs[j]) for j in range(n)] for i in range(n)]
    raise ValueError(Y)


@lru_cache(maxsize=None)
def c_Eg():
    """c(E_8^g) = -sum_{S subset E_8^g} c(S); S runs over Z[i] root lattices of rank <= 4.

    c(E_8^g) = -sum_{S subset E_8^g} c(S)。S は階数 4 以下の Z[i] ルート格子。"""
    R, g = build_eg()
    cA, cD, cG = _egf_coeff_vals(8)
    comps = [("A", 1), ("A", 2), ("A", 3), ("A", 4), ("D", 4), ("G", 2), ("G", 3), ("G", 4)]

    def rec(i, left):
        if i == len(comps):
            yield {}
            return
        Y = comps[i]
        for m in range(left // rank_of(Y) + 1):
            for rest in rec(i + 1, left - m * rank_of(Y)):
                d = dict(rest)
                if m:
                    d[Y] = m
                yield d

    tot = Fr(0)
    detail = {}
    for T in rec(0, 4):
        cs = Fr(1)
        aut = 1
        for Y, m in T.items():
            cy = {"A": cA, "D": cD, "G": cG}[Y[0]][Y[1]]
            cs *= cy ** m
            aut *= U_order(Y) ** m * factorial(m)
        emb = complex_embed_count(R, g, complex_gram(T)) if T else 1
        assert emb % aut == 0, (T, emb, aut)
        tot += cs * (emb // aut)
        detail[str(sorted(T.items()))] = emb // aut
    return -tot, detail


def c_table(rmax):
    cA, cD, cG = _egf_coeff_vals(max(rmax, 8))
    out = {}
    for Y in irreducibles(rmax):
        if Y[0] == "A":
            out[Y] = cA[Y[1]]
        elif Y[0] == "D":
            out[Y] = cD[Y[1]]
        elif Y[0] == "G":
            out[Y] = cG[Y[1]]
        elif Y[0] == "E":
            out[Y] = c_exceptional(Y[1])
        else:
            out[Y] = c_Eg()[0]
    return out


def all_types(rmax):
    irr = irreducibles(rmax)

    def rec(i, left):
        if i == len(irr):
            yield {}
            return
        Y = irr[i]
        for m in range(left // rank_of(Y) + 1):
            for rest in rec(i + 1, left - m * rank_of(Y)):
                d = dict(rest)
                if m:
                    d[Y] = m
                yield d

    return rec(0, rmax)


def _tup(x):
    """Lists back to tuples, recursively (JSON keeps tuples as lists).

    JSON はタプルをリストで持つので、再帰的にタプルに戻す。"""
    return tuple(_tup(y) for y in x) if isinstance(x, list) else x


def _cache_to_json():
    """The caches (_stats_store, _trans, _trans2) as a JSON value.

    キャッシュ（_stats_store、_trans、_trans2）を JSON の値にする。"""
    return {"stats": [[k, [[ell, d, w] for (ell, d), w in v.items()]] for k, v in _stats_store.items()],
            "trans": [[list(k), [[ns, [[c, w] for c, w in d.items()]] for ns, d in v.items()]] for k, v in _trans.items()],
            "trans2": [[list(k), [[ns, [[c, w] for c, w in d.items()]] for ns, d in v.items()]] for k, v in _trans2.items()]}


def _cache_from_json(data):
    """Inverse of _cache_to_json: (stats, trans, trans2) as dictionaries with tuple keys.

    _cache_to_json の逆。タプルをキーとする辞書 (stats, trans, trans2) を返す。"""
    stats = {_tup(k): {(ell, d): w for ell, d, w in v} for k, v in data["stats"]}
    trans = {_tup(k): {_tup(ns): {c: w for c, w in d} for ns, d in v} for k, v in data["trans"]}
    trans2 = {_tup(k): {_tup(ns): {c: w for c, w in d} for ns, d in v} for k, v in data["trans2"]}
    return stats, trans, trans2


def _save_cache(path):
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(_cache_to_json(), fh)
    os.replace(tmp, path)


def _load_cache(path):
    if os.path.exists(path):
        with open(path) as fh:
            stats, trans, trans2 = _cache_from_json(json.load(fh))
        _stats_store.update(stats)
        _trans.update(trans)
        _trans2.update(trans2)


class CpuLimit(Exception):
    pass


def solve(r, verbose=True, cpu_limit=None):
    """m(emptyset) for rank r. Progress in OUT/solve-r<r>.progress.json, resume data in OUT/cache-r<r>.json
    (transitions and statistics); the same call resumes. The result is written to OUT/solve-r<r>.json.

    階数 r の m(emptyset)。途中経過は OUT/solve-r<r>.progress.json、再開用は OUT/cache-r<r>.json（遷移と統計）。
    同じ呼び出しで再開できる。結果は OUT/solve-r<r>.json に書く。"""
    global _tick
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    prog = OUT / f"solve-r{r}.progress.json"
    cache = OUT / f"cache-r{r}.json"
    cpu_prev = 0.0
    if prog.exists():
        try:
            cpu_prev = json.load(open(prog)).get("cpu_s_total", 0.0)
        except ValueError:
            cpu_prev = 0.0
    _load_cache(cache)
    ct = c_table(r)
    if verbose:
        print("table of c:", {str(k): str(v) for k, v in ct.items()}, flush=True)
    types = list(all_types(r))
    types.sort(key=lambda T: sum(rank_of(Y) * m for Y, m in T.items()))
    ntot = len(types)
    keys = set()
    for T in types:
        bl = []
        for Y, m in T.items():
            bl += blocks2(Y) * m
        keys.add(tuple(sorted(bl)))
    state = {"i": 0, "last": 0.0, "last_ck": time.time(), "rank": 0, "stage": "start"}

    def report(status="running"):
        now = time.time()
        cpu = cpu_prev + time.process_time()
        done2 = sum(1 for k in keys if k in _stats_store)
        _write_json(prog, {"job": f"hermitian-mass-solve-r{r}", "stage": state["stage"], "state": status, "done": state["i"], "total": ntot, "unit": "root types T",
                           "local2_done": done2, "local2_total": len(keys), "unit2": "2-adic types (heartbeat counts stats, not transitions)",
                           "transitions_cached": len(_trans) + len(_trans2), "current_rank": state["rank"], "updated": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
                           "elapsed_s": round(now - t0, 1), "cpu_s_total": round(cpu, 1), "failed": 0})
        print(f"[{now - t0:.0f}s, CPU {cpu:.0f}s] types {state['i']}/{ntot}, 2-adic types {done2}/{len(keys)}, transitions {len(_trans) + len(_trans2)}, rank {state['rank']}", flush=True)

    def tick():
        now = time.time()
        if now - state["last"] > 60:
            state["last"] = now
            report()
        if now - state["last_ck"] > 600:
            state["last_ck"] = now
            _save_cache(cache)
        if cpu_limit is not None and cpu_prev + time.process_time() > cpu_limit:
            _save_cache(cache)
            report("stopped-cpu-limit")
            raise CpuLimit

    _tick = tick
    total = Fr(0)
    n0 = 0
    state["stage"] = "sum-over-T"
    for i, T in enumerate(types):
        state["i"] = i
        state["rank"] = sum(rank_of(Y) * m for Y, m in T.items())
        cs = Fr(1)
        aut = 1
        for Y, m in T.items():
            cs *= ct[Y] ** m
            aut *= U_order(Y) ** m * factorial(m)
        a = A_value(T, r)
        if a == 0:
            n0 += 1
        total += cs * a / aut
        tick()
    _save_cache(cache)
    state["i"] = ntot
    report("done")
    print(f"r={r}: {ntot} types ({n0} with A=0), m(empty) = {total}, seconds {time.time() - t0:.1f}", flush=True)
    _write_json(OUT / f"solve-r{r}.json", {"r": r, "types": ntot, "A_zero": n0, "m_empty": str(total), "m_empty_float": float(total),
                                           "local2_types": len(_stats_store), "elapsed_s": round(time.time() - t0, 1)})
    return total


# ---------------------------------------------------------------- The full table of masses per root system (triangular solve)
# ルート系ごとの質量の全表（三角に解く）


def _t_key(T):
    return tuple(sorted(T.items()))


def _t_mul(u, v):
    out = defaultdict(int)
    for a, x in u.items():
        for b, y in v.items():
            d = Counter(dict(a))
            d.update(dict(b))
            out[tuple(sorted(d.items()))] += x * y
    return dict(out)


def _single(*ys):
    return tuple(sorted(Counter(ys).items()))


def _alias(t, v):
    """Name of an irreducible component (handles the aliases D_2, D_3, G_1).

    既約の名前（D_2、D_3、G_1 の別名を処理する）。"""
    if t == "A":
        return _single(("A", v))
    if t == "D":
        if v == 2:
            return _single(("A", 1), ("A", 1))
        if v == 3:
            return _single(("A", 3))
        return _single(("D", v))
    if t == "G":
        if v == 1:
            return _single(("A", 1))
        return _single(("G", v))


@lru_cache(maxsize=None)
def sub_family(kind, n):
    """Table {type: count} of root sublattices on n vertices for kind = 'A', 'D', 'G' ('A' means A_{n-1}).

    kind = 'A', 'D', 'G' の頂点 n 個の上の部分ルート格子の表 {型: 個数}。A は A_{n-1}。"""
    if n == 0:
        return {(): 1}
    out = defaultdict(int)
    for v in range(1, n + 1):
        opts = {}
        if v == 1:
            opts[()] = 1
            if kind == "G":
                opts[_alias("G", 1)] = opts.get(_alias("G", 1), 0) + 1
        else:
            gain = {"A": 1, "D": 2, "G": 4}[kind] ** (v - 1)
            k = _alias("A", v - 1)
            opts[k] = opts.get(k, 0) + gain
            if kind in ("D", "G"):
                k = _alias("D", v)
                g2 = 1 if kind == "D" else 2 ** (v - 1)
                opts[k] = opts.get(k, 0) + g2
            if kind == "G":
                k = _alias("G", v)
                opts[k] = opts.get(k, 0) + 1
        rest = sub_family(kind, n - v)
        for a, x in _t_mul(opts, rest).items():
            out[a] += comb(n - 1, v - 1) * x
    return dict(out)


def _real_sub_exceptional(n):
    roots = e_roots(n)
    out = {}
    for T in ade_types_upto(n):
        if not T:
            out[()] = 1
            continue
        if T == {("E", n): 1}:
            out[_t_key(T)] = 1
            continue
        aut = 1
        for Y, m in T.items():
            aut *= aut_real(Y) ** m * factorial(m)
        emb = count_embeddings_real(roots, block_gram_real(T))
        if emb:
            out[_t_key(T)] = emb // aut
    return out


_sub_cache = {}


def sub_irr(Y):
    if Y in _sub_cache:
        return _sub_cache[Y]
    t, n = Y
    if t == "A":
        r = sub_family("A", n + 1)
    elif t == "D":
        r = sub_family("D", n)
    elif t == "G":
        r = sub_family("G", n)
    elif t == "E":
        r = _real_sub_exceptional(n)
    else:
        r = {}
        for k, v in c_Eg()[1].items():
            import ast

            T = tuple(sorted(ast.literal_eval(k)))
            if v:
                r[T] = v
        r[_single(("Eg", 8))] = 1
    _sub_cache[Y] = r
    return r


def roots_of(T):
    n = 0
    for Y, m in T:
        t, k = Y
        if t == "A":
            z = k * (k + 1)
        elif t == "D":
            z = 2 * k * (k - 1)
        elif t == "E":
            z = {6: 72, 7: 126, 8: 240}[k]
        if t in ("A", "D", "E"):
            n += 2 * z * m
        elif t == "G":
            n += 2 * (2 * k) * (2 * k - 1) * m
        else:
            n += 240 * m
    return n


def table(r):
    """m(R) exactly for every root system R.

    全部のルート系 R について m(R) を厳密に求める。"""
    cache = OUT / f"cache-r{r}.json"
    _load_cache(cache)
    types = [_t_key(T) for T in all_types(r)]
    rk = {T: sum(rank_of(Y) * m for Y, m in T) for T in types}
    Aval = {}
    autv = {}
    for T in types:
        Aval[T] = A_value(dict(T), r)
        a = 1
        for Y, m in T:
            a *= U_order(Y) ** m * factorial(m)
        autv[T] = a
    order = sorted(types, key=lambda T: (-rk[T], -roots_of(T)))
    M = {}
    residual = dict(Aval)  # A(T) - sum_{R done} m(R) r(T, R) / 処理済みの R について引く
    for R in order:
        val = residual[R] / autv[R]
        M[R] = val
        if val == 0:
            continue
        # table of root sublattices of R / R の部分ルート格子の表
        tb = {(): 1}
        for Y, m in R:
            for _ in range(m):
                tb = _t_mul(tb, sub_irr(Y))
        for T, cnt in tb.items():
            if T != R:
                residual[T] -= val * cnt * autv[T]
    return M
