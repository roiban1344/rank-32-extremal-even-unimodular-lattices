"""A joint linear program on the theta series of F and Z that excludes types p-(z, d)-s in rank 32 (p = 3, 5, 7).

The variables are the theta series of F and Z, in M_{n/2}(Gamma0(p), chi); the Fricke involution maps
them to theta_{sqrt(p) X^#} as in rank32.theta_lp_basis. Constraints:
  (a) theta_F, theta_Z = 1 + 0 q + (>= 0) (minimum >= 4); the coefficients of theta_{sqrt(p) F^#} and
      theta_{sqrt(p) Z^#} are >= 0, with constant term 1.
  (b) min Z^# >= 2(p - 1)/p: the coefficients q^1 .. q^(p-2) of theta_{sqrt(p) Z^#} vanish.
  (c) gluing: the vectors of norm 2j/p of a class a of F^# (q(a) = 2j/p mod 2) and the vectors of norm
      2 - 2j/p of the corresponding class of Z^# cannot both exist (their sum would be a root). One class
      contains at most m(alpha) = 1 + floor(alpha/(2 - alpha)) vectors of norm alpha (alpha >= 1), or 1
      (alpha < 1). Hence c_F[j]/m(2j/p) + c_Z[p - j]/m(2(p - j)/p) <= N(j), where N(j) is the number
      of classes with q(a) = 2j/p.
Truncation only relaxes the constraints, so infeasibility is a correct conclusion (the basis is checked by
its rank and by the linear relations, as in rank32.theta_lp_basis). The program is solved exactly over
the rationals (PPL).
  types_for(p, n)               the types that pass the genus conditions and the density sieve
  class_counts(p, d, s)         N(j)
  basis(pairs)                  the linearly independent pairs
  solve(p, d, n2, s, thF, thZ, N, extra)   None if infeasible, otherwise the ranges of some coefficients
Used by computation C2 of the first paper (code/produce/element_types/prime_order_type_lp.py).

階数 32 の型 p-(z, d)-s を、F と Z のテータ級数の同時線形計画で除く（p = 3, 5, 7）。

変数は F と Z のテータ級数（M_{n/2}(Gamma0(p), chi) の中）。Fricke 対合で theta_{sqrt(p) X^#} に
写ることは rank32.theta_lp_basis と同じ。制約：
  (a) theta_F, theta_Z = 1 + 0 q + (>= 0)（最小 >= 4）。theta_{sqrt(p) F^#}, theta_{sqrt(p) Z^#} の係数は
      >= 0 で、定数項は 1。
  (b) min Z^# >= 2(p - 1)/p：theta_{sqrt(p) Z^#} の q^1 .. q^(p-2) が 0。
  (c) 貼り合わせ：F^# の類 a（q(a) = 2j/p mod 2）のノルム 2j/p のベクトルと、対応する Z^# の類の
      ノルム 2 - 2j/p のベクトルは同時には存在しない（和がルートになる）。1 つの類のノルム alpha の
      ベクトルは高々 m(alpha) = 1 + floor(alpha/(2 - alpha)) 本（alpha >= 1）か 1 本（alpha < 1）。
      そこで c_F[j]/m(2j/p) + c_Z[p - j]/m(2(p - j)/p) <= N(j)。N(j) は q(a) = 2j/p の類の数。
打ち切りは制約を緩めるだけなので、実行不可能という結論は正しい（基底の正しさは、rank32.theta_lp_basis
と同じく階数と一次関係で確かめる）。線形計画は有理数で厳密に解く（PPL）。
  types_for(p, n)               種の条件と密度の篩を通る型
  class_counts(p, d, s)         N(j)
  basis(pairs)                  一次独立な組
  solve(p, d, n2, s, thF, thZ, N, extra)   実行不可能なら None、そうでなければいくつかの係数の範囲
第 1 の論文の計算 C2 が使う（code/produce/element_types/prime_order_type_lp.py）。
"""
from sage.all import QQ, matrix, MixedIntegerLinearProgram, kronecker_symbol

from . import genus_mass as ptmass
from . import theta_lp_basis as tb

PREC = 40  # the largest power of q used in the linear program / 線形計画で使う q の冪の上限

# Values of the field "verdict" of the output (code/produce/element_types/prime_order_type_lp.py) and of the
# ranges of solve().
# 出力の欄 verdict（code/produce/element_types/prime_order_type_lp.py）と solve() の範囲の値。
EXCLUDED = "excluded"           # the program is infeasible / 実行不可能
REMAINS = "remains"             # the program is feasible / 実行可能
UNDECIDED = "undecided"         # prefix: the blocks do not span the space / 接頭辞。ブロックが空間を張らない
NO_GENUS = "no genus"           # the character is not a character mod p / 指標が mod p の指標にならない
CHECK_FAILED = "check failed"   # the g_i do not satisfy the relations of the theta_i / g_i が一次関係を満たさない
UNBOUNDED = "unbounded"         # a range of solve() without bound / solve() の有界でない範囲


def m_alpha(alpha):
    """The largest number of vectors of norm alpha in one class: 1 + floor(alpha/(2 - alpha)) for alpha >= 1, 1 below.

    1 つの類のノルム alpha のベクトルの最大数：alpha >= 1 なら 1 + floor(alpha/(2 - alpha))、それ未満は 1。"""
    if alpha < 1:
        return 1
    return 1 + int(alpha / (2 - alpha))


def eta(a, p):
    """The Legendre symbol (a / p).

    Legendre 記号 (a / p)。"""
    return kronecker_symbol(a, p)


def class_counts(p, d, s):
    """N(j): the number of nonzero a in F^#/F with q(a) = 2j/p mod 2.

    Q = (p/2) q is a non-degenerate quadratic form on F_p^s whose discriminant is fixed by the unimodular
    component of sqrt(p) F^#.

    N(j)：F^#/F の 0 でない元 a のうち q(a) = 2j/p mod 2 のものの数。
    Q = (p/2) q は F_p^s 上の非退化 2 次形式で、その判別式は sqrt(p) F^# の単模成分から決まる。"""
    e2 = ptmass.eps2(p, d, s)
    detQ = e2 * eta(2, p) ** s  # eta(det Q)
    N = {}
    for j in range(p):
        if s == 0:
            N[j] = 0
            continue
        if s % 2 == 0:
            e = eta(-1, p) ** (s // 2) * e2 * 1  # eta(2)^s = 1
            if j == 0:
                N[j] = p ** (s - 1) + e * (p ** (s // 2) - p ** (s // 2 - 1)) - 1
            else:
                N[j] = p ** (s - 1) - e * p ** (s // 2 - 1)
        else:
            if j == 0:
                N[j] = p ** (s - 1) - 1
            else:
                N[j] = p ** (s - 1) + eta(j, p) * eta(-1, p) ** ((s - 1) // 2) * detQ * p ** ((s - 1) // 2)
    assert sum(N.values()) == p ** s - 1, (p, d, s, N)
    return N


def basis(pairs):
    """The linearly independent pairs (theta_i, g_i), in order.

    一次独立な組 (theta_i, g_i) を、順に選ぶ。"""
    out = []
    rows = []
    rank = 0
    for f, g, c, si in pairs:
        row = [f[j] for j in range(tb.P)]
        r = matrix(QQ, rows + [row]).rank()
        if r > rank:
            rows.append(row)
            out.append((f, g))
            rank = r
    return out


def solve(p, d, n2, s, thF, thZ, N, extra=None):
    """The joint linear program; None if infeasible, otherwise the ranges of some coefficients.

    extra is a list of additional equations (F or Z, power, value). A range that is unbounded is written
    "unbounded".

    同時線形計画。実行不可能なら None、実行可能なら係数の範囲を返す。
    extra は追加の等式 (F か Z, 冪, 値) の列。上下に有界でない範囲は "unbounded" と書く。"""
    def problem(obj, maximize):
        lp_ = MixedIntegerLinearProgram(maximization=maximize, solver="PPL")
        a = lp_.new_variable(nonnegative=False)
        b = lp_.new_variable(nonnegative=False)
        Ft = lambda j: sum(a[i] * thF[i][0][j] for i in range(len(thF)))
        Fg = lambda j: sum(a[i] * thF[i][1][j] for i in range(len(thF)))
        Zt = lambda j: sum(b[i] * thZ[i][0][j] for i in range(len(thZ)))
        Zg = lambda j: sum(b[i] * thZ[i][1][j] for i in range(len(thZ)))
        if thF:
            lp_.add_constraint(Ft(0) == 1)
            lp_.add_constraint(Ft(1) == 0)
            lp_.add_constraint(Fg(0) == 1)
            for j in range(2, PREC):
                lp_.add_constraint(Ft(j) >= 0)
            for j in range(1, PREC):
                lp_.add_constraint(Fg(j) >= 0)
        lp_.add_constraint(Zt(0) == 1)
        lp_.add_constraint(Zt(1) == 0)
        lp_.add_constraint(Zg(0) == 1)
        for j in range(2, PREC):
            lp_.add_constraint(Zt(j) >= 0)
        for j in range(1, PREC):
            if j <= p - 2:
                lp_.add_constraint(Zg(j) == 0)
            else:
                lp_.add_constraint(Zg(j) >= 0)
        if thF:
            for j in range(1, p):
                aj = QQ(2 * j) / p
                bj = 2 - aj
                lp_.add_constraint(Fg(j) * (1 / QQ(m_alpha(aj))) + Zg(p - j) * (1 / QQ(m_alpha(bj))) <= N[j])
        for part, j, v in (extra or []):
            lp_.add_constraint((Ft if part == "F" else Zt)(j) == v)
        if obj is None:
            lp_.set_objective(Zt(2))
        else:
            part, j = obj
            lp_.set_objective({"F": Ft, "Z": Zt, "Fg": Fg, "Zg": Zg}[part](j))
        return lp_
    try:
        problem(None, True).solve()
    except Exception as ex:
        if "no feasible" in str(ex).lower():
            return None
        if "unbounded" not in str(ex).lower():
            raise
    ranges = {}
    for key in ([("F", 2), ("Fg", 1), ("Fg", 2)] if thF else []) + [("Z", 2), ("Zg", p - 1)]:
        vals = []
        for mx in (False, True):
            try:
                vals.append(str(problem(key, mx).solve()))
            except Exception as ex:
                vals.append(UNBOUNDED if "unbounded" in str(ex).lower() else f"? {ex}")
        ranges[f"{key[0]}[{key[1]}]"] = vals
    return ranges



def types_for(p, n=32):
    """The types (z, d, s) for p in rank n that pass the genus conditions and the density sieve (minimum 4).

    階数 n の p の型 (z, d, s) のうち、種の条件と密度の篩（最小 4）を通るもの。"""
    out = []
    for z in range(1, n // (p - 1) + 1):
        d = n - z * (p - 1)
        for s in range(0, min(z, d) + 1):
            if (z - s) % 2:
                continue
            n2 = z * (p - 1)
            if not ptmass.genus_exists(p, n2, s):
                continue
            if d == 0 and s != 0:
                continue
            if d > 0 and not ptmass.genus_exists(p, d, s):
                continue
            if (d == 0 or tb.sieve.ok(d, p, s, 4)) and tb.sieve.ok(n2, p, s, 4):
                out.append((z, d, s))
    return out
