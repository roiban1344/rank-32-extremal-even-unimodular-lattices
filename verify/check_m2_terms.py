"""Recompute the embedding averages A_16(T) of identity (C.2) (Appendix C).

Usage, from the repository root:
    python3 verify/check_m2_terms.py OUT.json

Reads certificates/gaussian-rootless-mass.json and needs only the Python
standard library.  Exit status 0 means that every check passed.

What is recomputed, in exact rational arithmetic:
  * the index of root types: all multisets of the irreducible root lattices of
    Proposition C.1 of total rank at most 16 (12934 of them, with the counts by
    rank of §C.2); every term of the certificate is one of them, each once;
  * for every irreducible Y, from the Cartan matrix of its root system over Z
    (A_n[i], D_n[i], E_n[i] by base change; G_m from D_2m with g = i; E_8^g
    from E_8): det_O(Y), the 2-adic Jordan blocks of Table C1 and the odd
    p-parts of the discriminant module D(Y);
  * S_t and S'_t of (C.6) from Bernoulli numbers and the generalized Bernoulli
    numbers of chi = chi_{-4} (Bernoulli polynomials, no table of L-values);
    the passage from (C.3) (Theorem C.4, Lemma C.5) to (C.6) is checked as an
    exact rational identity for every t (the bookkeeping E = 0 of the proof of
    Theorem C.8), and the local factors of (C.6) are checked against the
    quotients of the primitive densities of Proposition C.7;
  * the odd local factors f_p(T, t) from the closed sums of §C.6 ("Closed sums
    at odd primes"), and for (p, a, b) = (3, 0, 2) by enumerating the isotropic
    submodules of (Z_3[i]/9)^2;
  * the 2-adic factor f_2(T, t) of (C.6) from the saved distribution of
    (log_2 |C|, d(C)) for the block multiset of T;
  * A_16(T) = S'_t det(T)^t prod_p f_p(T, t) (C.6) for every term, compared with
    the saved embedding_average; c(T) and |U(T)| from the irreducible
    coefficients and Proposition C.1; the term c(T) A_16(T) / |U(T)| and the
    sum (C.2), compared with term and rootless_mass;
  * as a consistency check, (C.2) for r = 4, 8, 12 from the same data
    (Table C2: 0, 0, 1/2012774400);
  * that local_statistics has one row for every block multiset the terms need
    (rows that no term uses are reported), and light invariants of each row.

What remains input:
  * the local 2-adic statistics (local_statistics): the integer counts of
    isotropic submodules C of D_2(T) by (log_2 |C|, d(C)), computed by the
    finite transition of Appendix C.6 (Lemma C.10); they are used, not rederived;
  * the irreducible coefficients c(Y): those of A_n, D_n, G_m are compared with
    the closed formulas of §C.3 (whose derivation by set-partition generating
    functions is not redone), those of E_6[i], E_7[i], E_8[i], E_8^g (counted
    from explicit roots in the paper) are compared with the values stated in
    §C.3 only;
  * |Aut E_n| (textbook values) and the centralizer orders 96 (G_2) and
    46080 (E_8^g) of Proposition C.1.
"""
import itertools
import json
try:
    import resource                      # POSIX only; used to report the memory
except ImportError:
    resource = None
import sys
import time
from collections import defaultdict
from fractions import Fraction as Q
from math import comb, factorial, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'certificates' / 'gaussian-rootless-mass.json'
R = 16                      # the rank of the theta lattices (M_2 = mu_16^empty)


# ------------------------------------------------------------- L-values, S_t, S'_t

def bernoulli_numbers(n):
    """B_0, ..., B_n with B_1 = -1/2 (sum_{k<m+1} binom(m+1, k) B_k = 0)."""
    B = [Q(1)] + [Q(0)] * n
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


BERN = bernoulli_numbers(R + 1)


def bernoulli_poly(n, x):
    return sum(comb(n, k) * BERN[k] * x ** (n - k) for k in range(n + 1))


def chi(a):
    """chi_{-4}(a), the character of Q(i)."""
    return 0 if a % 2 == 0 else (1 if a % 4 == 1 else -1)


def L_one_minus(j):
    """L(1 - j, chi^j) for j >= 1 (chi^j = chi for odd j, trivial for even j).

    L(1 - j, psi) = -B_{j, psi} / j with the generalized Bernoulli number
    B_{j, psi} = f^{j-1} sum_{a=1}^{f} psi(a) B_j(a / f) of conductor f;
    for the trivial character this is zeta(1 - j) = -B_j / j (j >= 2)."""
    if j % 2 == 0:
        return -BERN[j] / j
    f = 4
    b = Q(f) ** (j - 1) * sum(chi(a) * bernoulli_poly(j, Q(a, f)) for a in range(1, f + 1))
    return -b / j


def euler_numbers(n):
    """Euler numbers E_0, E_2, ... (secant), only for a cross-check."""
    E = [Q(0)] * (n + 1)
    E[0] = Q(1)
    for m in range(2, n + 1, 2):
        E[m] = -sum(comb(m, k) * E[k] for k in range(0, m, 2))
    return E


def S(t):
    """S_t = 2^{1-t} prod_{j=1}^t |L(1-j, chi^j)| (Theorem C.8)."""
    v = Q(2) ** (1 - t)
    for j in range(1, t + 1):
        v *= abs(L_one_minus(j))
    return v


def S_prime(t):
    """S'_t = 2^{-t [t odd]} S_t, S'_0 = 1 (Theorem C.8)."""
    if t == 0:
        return Q(1)
    return S(t) / (2 ** t if t % 2 else 1)


def check_L_and_C6_constant():
    """Checks of the archimedean bookkeeping.  Returns a list of failures."""
    bad = []
    E = euler_numbers(R)
    # L(0, chi) = 1/2 (= 2 h / w for Q(i)); L(1 - j, chi) = E_{j-1} / 2 for odd j.
    if L_one_minus(1) != Q(1, 2):
        bad.append('L(0, chi) != 1/2')
    for j in range(1, R + 1, 2):
        if L_one_minus(j) != E[j - 1] / 2:
            bad.append(f'L(1-{j}, chi) != E_{j - 1}/2')
    # From (C.3) to (C.6).  With mu_r = S'_r, the functional equations
    #   L(j, chi^j) = |L(1-j, chi^j)| (2 pi)^j / (2 (j-1)! 2^{(2j-1)[j odd]})
    # and prod_p alpha_p(T) = prod_p f_p(T, t) * prod_{j=t+1}^r L(j, chi^j)^{-1},
    # (C.3) reads A_r(T) = K * det(T)^t prod_p f_p with K given below; the powers
    # of pi cancel since sum_{j=t+1}^r j = n (2r - n + 1) / 2.  (C.6) says K = S'_t.
    for r in (4, 8, 12, 16):
        for n in range(0, r):
            t = r - n
            if sum(range(t + 1, r + 1)) * 2 != n * (2 * r - n + 1):
                bad.append(f'pi powers r={r} n={n}')
            K = S_prime(r) * Q(2) ** n
            for j in range(t + 1, r + 1):          # 1 / prod_{j=0}^{n-1} Gamma(r - j)
                K /= factorial(j - 1)
            for j in range(t + 1, r + 1):          # 1 / L(j, chi^j) without (2 pi)^j
                K *= Q(2 * factorial(j - 1) * 2 ** ((2 * j - 1) * (j % 2)), 2 ** j) / abs(L_one_minus(j))
            if K != S_prime(t):
                bad.append(f'(C.3) -> (C.6) constant r={r} t={t}')
    return bad


# ------------------------------------------------------------- local factors

def ell_odd(p, d, t):
    """prod_{j=t-d+1}^t (1 - chi(p)^j p^{-j}), zero if the range contains j = 0."""
    v = Q(1)
    for j in range(t - d + 1, t + 1):
        if j == 0:
            return Q(0)
        if j > 0:
            v *= 1 - Q(chi(p) ** j, p ** j)
    return v


def ell_two(d, t):
    """prod_{j even, t-d+2 <= j <= t} (1 - 2^{-j}) (J_2 of (C.6)), zero if j = 0 occurs."""
    v = Q(1)
    for j in range(t - d + 2, t + 1):
        if j % 2 == 0:
            if j == 0:
                return Q(0)
            if j > 0:
                v *= 1 - Q(1, 2 ** j)
    return v


def alpha_prim_two(n, d, r):
    """Proposition C.7 at p = 2: prod_{k=0}^{(n+d)/2-1} (1 - 4^k 2^{-r})."""
    v = Q(1)
    for k in range((n + d) // 2):
        v *= 1 - Q(4 ** k, 2 ** r)
    return v


def alpha_prim_odd(p, n, d, r):
    """Proposition C.7 at odd p: prod_{j=r-n-d+1}^r (1 - chi(p)^j p^{-j})."""
    lo = r - n - d + 1
    if lo <= 0:                  # the factor j = 0 vanishes: no primitive embeddings
        return Q(0)
    v = Q(1)
    for j in range(lo, r + 1):
        v *= 1 - Q(chi(p) ** j, p ** j)
    return v


def check_local_factors_against_C7():
    """f_p of (C.6) divides alpha_{p,prim}(T_C) by the value for p not dividing det T."""
    bad = []
    for n in range(1, R + 1):
        t = R - n
        for d in range(0, R + 1):
            if (n + d) % 2 == 0 and n + d <= 2 * R:
                # the unramified value: d = 0 factors j = r - 2k, t + 1 <= j <= r, j even
                base = Q(1)
                for j in range(t + 1, R + 1):
                    if j % 2 == 0:
                        base *= 1 - Q(1, 2 ** j)
                if alpha_prim_two(n, d, R) / base != ell_two(d, t):
                    bad.append(f'p=2 n={n} d={d}')
            for p in (3, 5, 7, 11, 13, 17):
                if d <= n:
                    base = alpha_prim_odd(p, n, 0, R)
                    if alpha_prim_odd(p, n, d, R) / base != ell_odd(p, d, t):
                        bad.append(f'p={p} n={n} d={d}')
    return bad


def gauss_binom(n, k, q):
    """Number of k-subspaces of F_q^n."""
    if k < 0 or k > n:
        return 0
    num = den = 1
    for i in range(k):
        num *= q ** (n - i) - 1
        den *= q ** (i + 1) - 1
    assert num % den == 0
    return num // den


def unitary_isotropic(p, a, k):
    """N(a, k): totally isotropic k-subspaces of a non-degenerate Hermitian a-space over F_{p^2}."""
    if 2 * k > a:
        return 0
    num = den = 1
    for j in range(k):
        s = a - 2 * j
        num *= (p ** s - (-1) ** s) * (p ** (s - 1) - (-1) ** (s - 1))
        den *= p ** (2 * (j + 1)) - 1
    assert num % den == 0
    return num // den


def distribution_9x9():
    """Isotropic O-submodules C of (Z[i]/9)^2 (p = 3, (a, b) = (0, 2), from A_8[i]^2),
    as {(log_3 |C|, d(C)): count}.

    D_3(A_8[i]) = O/9 with h(y, y) = 8/9; on two copies h(x, y) = 8 (x_1 conj y_1 +
    x_2 conj y_2) / 9 mod O, and the conditions h = 0 mod O do not see the unit 8.
    As in §C.6: start from C = 0, adjoin x and i x for isotropic x in C^perp \\ C and
    take the additive span, discarding duplicates.  d(C) is read off from
    |(C^perp / C)[3]| = 9^{d(C)}."""
    els = list(itertools.product(range(9), repeat=4))      # (Re x1, Im x1, Re x2, Im x2)

    def h(x, y):  # numerator of x_1 conj(y_1) + x_2 conj(y_2), as (re, im) mod 9
        return ((x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3]) % 9,
                (x[1] * y[0] - x[0] * y[1] + x[3] * y[2] - x[2] * y[3]) % 9)

    def times_i(x):
        return ((-x[1]) % 9, x[0], (-x[3]) % 9, x[2])

    def span(base, gens):
        cur = set(base)
        frontier = list(cur)
        while frontier:
            new = []
            for c in frontier:
                for g in gens:
                    w = tuple((u + v) % 9 for u, v in zip(c, g))
                    if w not in cur:
                        cur.add(w)
                        new.append(w)
            frontier = new
        return frozenset(cur)

    zero = (0, 0, 0, 0)
    start = (frozenset([zero]), ())
    seen = {start[0]}
    stack = [start]
    dist = defaultdict(int)
    while stack:
        C, gens = stack.pop()
        perp = [x for x in els if all(h(x, g) == (0, 0) for g in gens)]
        tor = sum(1 for x in perp if tuple(3 * u % 9 for u in x) in C)
        assert tor % len(C) == 0 and len(perp) * len(C) == 9 ** 4
        size, d = tor // len(C), 0
        while 9 ** d < size:
            d += 1
        assert 9 ** d == size
        k = 0
        while 3 ** k < len(C):
            k += 1
        assert 3 ** k == len(C)
        dist[(k, d)] += 1
        for x in perp:
            if x in C or h(x, x) != (0, 0):
                continue
            nC = span(C, [x, times_i(x)])
            if nC not in seen:
                seen.add(nC)
                stack.append((nC, gens + (x, times_i(x))))
    return dict(dist)


_dist9 = None


def f_odd(p, a, b, t):
    """f_p(T, t) of (C.6) for D_p(T) = (R_p/p)^a + (R_p/p^2)^b (§C.6, closed sums)."""
    def term(log_c, d):
        return Q(1, p ** (log_c * t)) * ell_odd(p, d, t)

    total = Q(0)
    if p % 4 == 1:
        # two dual a-spaces over F_p; C = (C_+, C_-), C_- in C_+^perp
        assert b == 0
        for k1 in range(a + 1):
            for k2 in range(a - k1 + 1):
                total += gauss_binom(a, k1, p) * gauss_binom(a - k1, k2, p) * term(k1 + k2, a - k1 - k2)
        return total
    if b == 0:
        for k in range(a // 2 + 1):
            total += unitary_isotropic(p, a, k) * term(2 * k, a - 2 * k)
        return total
    if b == 1:
        # graph transition with the cyclic summand B = R_p / p^2 (table of §C.6)
        for k in range(a // 2 + 1):
            ap = a - 2 * k
            c0 = unitary_isotropic(p, a, k)
            total += c0 * term(2 * k, ap + 1)                       # Z = 0
            total += c0 * term(2 * k + 2, ap)                       # Z = pB, zero map
            if ap >= 2:                                             # Z = pB, isotropic image
                total += c0 * unitary_isotropic(p, ap, 1) * (p ** 2 - 1) * term(2 * k + 2, ap - 1)
        return total
    assert (p, a, b) == (3, 0, 2)
    global _dist9
    if _dist9 is None:
        _dist9 = distribution_9x9()
    return sum((cnt * term(k, d) for (k, d), cnt in _dist9.items()), Q(0))


def f_two(stats, t):
    """f_2(T, t) of (C.6) from the counts {(ell, d): number of C} with |C| = 2^ell."""
    return sum((Q(cnt, 2 ** (ell * t)) * ell_two(d, t) for (ell, d), cnt in stats.items()), Q(0))


# ------------------------------------------------------------- root lattices from their Cartan matrices

def cartan(kind, n):
    """Cartan matrix (Gram matrix of the simple roots, norm 2) of A_n, D_n, E_n (Bourbaki numbering)."""
    if kind == 'A':
        edges = [(k, k + 1) for k in range(n - 1)]
    elif kind == 'D':
        edges = [(k, k + 1) for k in range(n - 2)] + [(n - 3, n - 1)]
    elif kind == 'E':
        edges = [(0, 2), (2, 3), (3, 4), (1, 3)] + [(k, k + 1) for k in range(4, n - 1)]
    else:
        raise ValueError(kind)
    C = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
    for i, j in edges:
        C[i][j] = C[j][i] = -1
    return C


def det_and_inverse(M):
    n = len(M)
    A = [[Q(x) for x in row] + [Q(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    det = Q(1)
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return det, [row[n:] for row in A]


def discriminant_form(kind, n):
    """D(C) = C^#/C of the Z-root lattice: its elements (root coordinates mod 1) generated by
    the fundamental weights, with q(x) = x^T C x mod 2 and b(x, y) = x^T C y mod 1."""
    C = cartan(kind, n)
    det, Cinv = det_and_inverse(C)
    assert det.denominator == 1 and det > 0

    def red(x):
        return tuple(v - (v.numerator // v.denominator) for v in x)

    def add(x, y):
        return red(tuple(u + v for u, v in zip(x, y)))

    def form(x, y):
        return sum(x[i] * C[i][j] * y[j] for i in range(n) for j in range(n))

    zero = tuple(Q(0) for _ in range(n))
    group = {zero}
    for w in Cinv:                       # rows of C^{-1}: the fundamental weights
        w = red(w)
        frontier = list(group)
        while frontier:
            new = []
            for g in frontier:
                s = add(g, w)
                if s not in group:
                    group.add(s)
                    new.append(s)
            frontier = new
    assert len(group) == det, (kind, n)

    def order(x):
        k, y = 1, x
        while y != zero:
            y, k = add(y, x), k + 1
        return k

    def q(x):
        return form(x, x) % 2

    def b(x, y):
        return form(x, y) % 1

    return int(det), group, order, q, b, add


def jordan_two_and_odd(kind, n):
    """For the base change C[i] of the Z-root lattice C = kind_n: det_O, the 2-adic blocks of
    Table C1 and the odd p-parts {p: exponent} of D(C[i]) = D(C) (x) O."""
    det, group, order, q, b, add = discriminant_form(kind, n)
    elements = sorted(group, key=order)
    blocks, odd = [], {}
    gen = elements[-1]
    if order(gen) == det:                # cyclic Z/N: p-parts Z/p^e, base change O/p^e
        N = det
        e2, m = 0, N
        while m % 2 == 0:
            m //= 2
            e2 += 1
        if e2:
            y = gen
            for _ in range(m - 1):
                y = add(y, gen)          # y = m * gen generates the 2-part Z/2^e
            u = q(y) * 2 ** e2           # q(y) = u / 2^e mod 2, u odd, u mod 4 (Lemma C.9(ii))
            assert u.denominator == 1 and u % 2 == 1
            blocks.append(('d', e2, int(u) % 4))
        p = 3
        while m > 1:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            if e:
                odd[p] = e
            p += 2
    else:                                # (Z/2)^2 (D_n, n even)
        assert det == 4 and all(order(x) <= 2 for x in group)
        s, c = [x for x in elements if order(x) == 2][:2]
        if all(q(x).denominator == 1 for x in group):
            # q integral on D: 2-adically C = unimodular + 2 (even plane); over O_2 the
            # plane of scale pi^2 is 2H (Lemma C.9(iii)).
            assert b(s, c) == Q(1, 2)
            blocks.append(('H', 1))
        else:
            # find an orthogonal basis of two elements of q in 1/2 + Z: <2u> + <2u'>
            pair = next((x, y) for x in group for y in group
                        if x != y and order(x) == order(y) == 2 and b(x, y) == 0
                        and q(x).denominator == 2 and q(y).denominator == 2)
            for x in pair:
                blocks.append(('d', 1, int(2 * q(x)) % 4))
    return det, blocks, odd


def gaussian_D2m(m):
    """G_m: the real lattice D_2m in Z^2m, z_k = x_{2k-1} + i x_{2k}, g = multiplication by i.

    G_m = {z in O^m : sum z_k in pi O} is exactly D_2m, since a + b i = a + b mod pi.
    D(D_2m) = {0, v, s, c}; returns the 2-adic block of D(G_m) as an O-module."""
    dim = 2 * m
    half = Q(1, 2)

    def in_D(x):
        return all(v.denominator == 1 for v in x) and sum(x) % 2 == 0

    def g(x):
        out = []
        for k in range(m):
            out += [-x[2 * k + 1], x[2 * k]]
        return tuple(out)

    def sub(x, y):
        return tuple(u - v for u, v in zip(x, y))

    def q(x):
        return sum(v * v for v in x) % 2

    zero = tuple(Q(0) for _ in range(dim))
    v = tuple(Q(int(k == 0)) for k in range(dim))
    s = tuple(half for _ in range(dim))
    c = tuple(half - (1 if k == 0 else 0) for k in range(dim))
    reps = [zero, v, s, c]
    assert all(not in_D(sub(x, y)) for x, y in itertools.combinations(reps, 2))

    def cls(x):
        return next(i for i, y in enumerate(reps) if in_D(sub(x, y)))

    gperm = [cls(g(x)) for x in reps]
    if gperm == [0, 1, 2, 3]:
        # g = 1 on D, so pi = 1 + g kills D: D = (O/pi)^2, a block P(0, t) (scale pi, rank 2,
        # Lemma C.9(iv)); t = 0 iff q has a non-zero isotropic element.
        t = 0 if any(q(x) == 0 for x in reps[1:]) else 1
        return [('P', 0, t)]
    # g exchanges s and c: D = O/2 = O/pi^2 generated by s, a block <2u> with q(s) = u/2.
    assert gperm == [0, 1, 3, 2]
    u = 2 * q(s)
    assert u.denominator == 1 and u % 2 == 1
    return [('d', 1, int(u) % 4)]


AUT_E = {6: 2 * 51840, 7: 2903040, 8: 696729600}      # |Aut E_n| (|W(E_6)| = 51840, -1 not in W(E_6))
EXCEPTIONAL_C = {'E6': 7680, 'E7': -233280, 'E8': -28183680, 'Eg8': -5184}   # values stated in §C.3


def irreducibles(rmax):
    """Proposition C.1, as (name, n) with Gaussian rank."""
    out = [('A', n) for n in range(1, rmax + 1)]
    out += [('D', n) for n in range(4, rmax + 1)]
    out += [('E', n) for n in (6, 7, 8) if n <= rmax]
    out += [('G', m) for m in range(2, rmax + 1)]
    if rmax >= 4:
        out.append(('Eg', 8))
    return out


def rank(Y):
    return 4 if Y[0] == 'Eg' else Y[1]


def unitary_order(Y):
    """|U(Y)| of Proposition C.1."""
    kind, n = Y
    if kind == 'A':
        return 4 if n == 1 else 4 * factorial(n + 1)
    if kind == 'D':
        return 2304 if n == 4 else 2 ** (n + 1) * factorial(n)
    if kind == 'E':
        return 2 * AUT_E[n]
    if kind == 'G':
        return 96 if n == 2 else 4 ** n * factorial(n)
    return 46080


def closed_coefficient(Y):
    """c(Y) from the closed formulas of §C.3 (None for the exceptional types)."""
    kind, n = Y
    if kind == 'A':
        return (-1) ** n * factorial(n)
    if kind == 'D':
        return (-1) ** n * factorial(n - 1) * (2 ** (n - 1) - 1)
    if kind == 'G':
        return (-1) ** n * 2 ** (n - 1) * factorial(n - 1)
    return None


def local_data(Y):
    """(det_O(Y), 2-adic blocks, {odd p: exponent}) of an irreducible root lattice over O."""
    kind, n = Y
    if kind in ('A', 'D', 'E'):
        # C[i] = C (x) O: its Hermitian Gram in the basis of C is the Cartan matrix of C.
        return jordan_two_and_odd(kind, n)
    if kind == 'G':
        det_real = discriminant_form('D', 2 * n)[0]
        blocks = gaussian_D2m(n)
    else:                                  # E_8^g: the real lattice is E_8
        det_real = discriminant_form('E', 8)[0]
        blocks = []
    # det_Z of the trace form is det_O^2 (realification of a Hermitian Gram matrix).
    det = isqrt(det_real)
    assert det * det == det_real
    return det, blocks, {}


def all_types(rmax, irr):
    """All multisets {Y: multiplicity} of total rank <= rmax, as sorted tuples."""
    def rec(i, left):
        if i == len(irr):
            yield ()
            return
        Y = irr[i]
        for m in range(left // rank(Y) + 1):
            for rest in rec(i + 1, left - m * rank(Y)):
                yield (((Y, m),) if m else ()) + rest
    return [tuple(sorted(t)) for t in rec(0, rmax)]


# ------------------------------------------------------------- A_r(T) by (C.6)

class Averages:
    def __init__(self, stats, local):
        self.stats = stats             # {sorted block tuple: {(ell, d): count}}
        self.local = local             # {Y: (det, blocks, odd)}
        self.used = set()
        self.missing = set()
        self._f2 = {}
        self._fp = {}

    def A(self, T, r):
        """A_r(T) = S'_t det(T)^t prod_p f_p(T, t) (C.6), t = r - n; T = ((Y, mult), ...)."""
        n = sum(rank(Y) * m for Y, m in T)
        t = r - n
        assert t >= 0
        det, blocks, odd = 1, [], defaultdict(lambda: [0, 0])
        for Y, m in T:
            dY, bY, oY = self.local[Y]
            det *= dY ** m
            blocks += bY * m
            for p, e in oY.items():
                assert e in (1, 2)
                odd[p][e - 1] += m
        value = S_prime(t) * Q(det) ** t
        for p, (a, b) in sorted(odd.items()):
            key = (p, a, b, t)
            if key not in self._fp:
                self._fp[key] = f_odd(p, a, b, t)
            value *= self._fp[key]
        key = tuple(sorted(blocks))
        if key not in self.stats:
            self.missing.add(key)
            return None
        self.used.add(key)
        if (key, t) not in self._f2:
            self._f2[(key, t)] = f_two(self.stats[key], t)
        return value * self._f2[(key, t)]


def main():
    start = time.time()
    out_path = Path(sys.argv[1])
    with open(CERT) as fh:
        data = json.load(fh)
    failures = []

    def fail(msg):
        failures.append(msg)

    for msg in check_L_and_C6_constant() + check_local_factors_against_C7():
        fail(msg)

    irr = irreducibles(R)
    local = {Y: local_data(Y) for Y in irr}

    # local statistics
    stats = {}
    for row in data['local_statistics']:
        key = tuple(sorted(tuple(b) for b in row['blocks']))
        cnt = {}
        for ell, d, c in row['counts']:
            if (ell, d) in cnt or not (isinstance(c, int) and c > 0):
                fail(f'local_statistics {key}: bad entry {(ell, d, c)}')
            cnt[(ell, d)] = c
        if key in stats:
            fail(f'local_statistics {key}: duplicate')
        # C = 0 is the only submodule of order 1, with d = d(D_2(T)) = sum of block ranks
        d0 = sum(1 if b[0] == 'd' else 2 for b in key)
        if {k: v for k, v in cnt.items() if k[0] == 0} != {(0, d0): 1}:
            fail(f'local_statistics {key}: ell = 0 row is not {{(0, {d0}): 1}}')
        stats[key] = cnt

    # coefficients
    coeff = {}
    for name, value in data['irreducible_coefficients'].items():
        coeff[name] = Q(value)
    names = {f'{k}{n}': (k, n) for k, n in irr}
    if set(coeff) != set(names):
        fail('irreducible_coefficients: names differ from Proposition C.1')
    for name, Y in names.items():
        closed = closed_coefficient(Y)
        expect = closed if closed is not None else EXCEPTIONAL_C[name]
        if coeff.get(name) != expect:
            fail(f'coefficient c({name}) = {coeff.get(name)}, expected {expect}')

    # the index of root types
    index = all_types(R, irr)
    by_rank = [0] * (R + 1)
    for T in index:
        by_rank[sum(rank(Y) * m for Y, m in T)] += 1
    if len(index) != 12934 or [by_rank[k] for k in (0, 4, 8, 12, 16)] != [1, 12, 113, 802, 4663]:
        fail(f'type index: {len(index)} types, by rank {by_rank}')
    if by_rank != data['counts_by_rank']:
        fail('counts_by_rank differ from the type index')

    avg = Averages(stats, local)
    seen = set()
    mismatch_avg = mismatch_coef = mismatch_order = mismatch_term = 0
    total = Q(0)
    nonzero = 0
    for row in data['terms']:
        T = tuple(sorted(((tuple(Y)[0], tuple(Y)[1]), m) for Y, m in row['type']))
        if T in seen or any(Y not in local for Y, _ in T):
            fail(f'type {T}: repeated or unknown')
            continue
        seen.add(T)
        if row['rank'] != sum(rank(Y) * m for Y, m in T):
            fail(f'type {T}: rank')
        c, order = Q(1), 1
        for Y, m in T:
            c *= coeff[f'{Y[0]}{Y[1]}'] ** m
            order *= unitary_order(Y) ** m * factorial(m)
        a = avg.A(T, R)
        if a is None:
            fail(f'type {T}: no local statistics for its blocks')
            continue
        nonzero += a != 0
        if a != Q(row['embedding_average']):
            mismatch_avg += 1
            fail(f'type {T}: A_16 = {a}, certificate {row["embedding_average"]}')
        if c != Q(row['coefficient']):
            mismatch_coef += 1
            fail(f'type {T}: c(T) = {c}, certificate {row["coefficient"]}')
        if order != row['group_order']:
            mismatch_order += 1
            fail(f'type {T}: |U(T)| = {order}, certificate {row["group_order"]}')
        term = c * a / order
        if term != Q(row['term']):
            mismatch_term += 1
            fail(f'type {T}: term {term}, certificate {row["term"]}')
        total += term
    if seen != set(index):
        fail(f'terms cover {len(seen & set(index))} of {len(index)} types')
    if total != Q(data['rootless_mass']):
        fail(f'sum (C.2) = {total}, certificate {data["rootless_mass"]}')
    # (The 11780 zero masses of §C.7 are masses of exact root systems, after Moebius
    # inversion over the sublattice poset, not zero averages; they are not checked here.)

    # Table C2 for r = 4, 8, 12 from the same data
    table_c2 = {}
    for r, expect in ((4, Q(0)), (8, Q(0)), (12, Q(1, 2012774400))):
        types = all_types(r, irreducibles(r))
        m = Q(0)
        for T in types:
            a = avg.A(T, r)
            if a is None:
                fail(f'r={r} type {T}: no local statistics')
                continue
            c, order = Q(1), 1
            for Y, k in T:
                c *= coeff[f'{Y[0]}{Y[1]}'] ** k
                order *= unitary_order(Y) ** k * factorial(k)
            m += c * a / order
        table_c2[r] = {'types': len(types), 'rootless_mass': str(m)}
        if m != expect or len(types) != {4: 22, 8: 254, 12: 2035}[r]:
            fail(f'Table C2 r={r}: {len(types)} types, mass {m}')

    used_16 = {tuple(sorted(sum((list(local[Y][1]) * m for Y, m in T), []))) for T in index}
    unused = sorted(set(stats) - used_16)
    if avg.missing:
        fail(f'{len(avg.missing)} block multisets missing from local_statistics')

    report = {
        'certificate': str(CERT.relative_to(ROOT)),
        'terms': len(data['terms']),
        'types_in_index': len(index),
        'embedding_averages_recomputed': len(seen),
        'embedding_averages_nonzero': nonzero,
        'mismatches': {'embedding_average': mismatch_avg, 'coefficient': mismatch_coef,
                       'group_order': mismatch_order, 'term': mismatch_term},
        'rootless_mass': str(total),
        'rootless_mass_agrees': total == Q(data['rootless_mass']),
        'table_C2': table_c2,
        'local_statistics_rows': len(stats),
        'local_statistics_used_rank16': len(used_16 & set(stats)),
        'local_statistics_unused': [list(map(list, k)) for k in unused],
        'local_statistics_missing': [list(map(list, k)) for k in sorted(avg.missing)],
        'odd_factor_cases': len(avg._fp),
        'distribution_9x9': {f'{k},{d}': v for (k, d), v in sorted((_dist9 or {}).items())},
        'irreducible_local_data': {f'{k}{n}': {'det': local[(k, n)][0], 'blocks2': local[(k, n)][1],
                                               'odd': local[(k, n)][2]} for k, n in irr},
        'failures_count': len(failures),
        'failures': failures[:50],
        'ok': not failures,
        'runtime_seconds': round(time.time() - start, 2),
        'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource else None,
        'inputs_not_recomputed': ['local_statistics (finite transition of Appendix C.6)',
                                  'c(E_6[i]), c(E_7[i]), c(E_8[i]), c(E_8^g) (compared with §C.3)',
                                  'closed formulas for c(A_n), c(D_n), c(G_m) (§C.3)',
                                  '|Aut E_n|, |U(G_2)| = 96, |U(E_8^g)| = 46080 (Proposition C.1)'],
    }
    tmp = out_path.with_name(out_path.name + '.tmp')
    tmp.write_text(json.dumps(report, indent=2) + '\n')
    tmp.replace(out_path)
    if failures:
        print(f'FAIL: {len(failures)} problem(s); first: {failures[0]}', flush=True)
        sys.exit(1)
    print(f'All {len(seen)} embedding averages A_16(T) recomputed and equal to the certificate; '
          f'sum (C.2) = {total}; {len(unused)} unused local rows; '
          f'{report["runtime_seconds"]} s', flush=True)


if __name__ == '__main__':
    main()
