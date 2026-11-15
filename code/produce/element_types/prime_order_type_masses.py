"""Masses of the genera of F and Z, and the mass of pairs (L, sigma), for every type p-(z, d)-s in rank n.

Paper: computation C2 of the first paper (Table 11), the table of types of elements of odd prime order.
The formulas are in rank32.genus_mass.
Modes:
  test        checks the formulas: the values of [HKN09, Theorem 1] for p = 3; mass_U in rank 1 and 2
              against 1/|Aut| of Z[zeta_p]^z; p = 3, z = 6, 7 against Z[omega]^6 + K_12 and
              Z[omega]^7 + Z[omega] + K_12; for E8, the sum of M over the types against the number of
              elements of order p of W(E8) (GAP, through Sage); in rank 16, against E8 + E8 and D16^+;
              in rank 24, p = 13, against A_12^2 and the Leech lattice.
  table N     one row per type p-(z, d)-s with p <= 31 in rank N whose genera exist: whether it passes
              the density sieve (minimum 4), mass_F, mass_H, |O|, M (exact fractions and floats).
Output: --out DIR (default results/C2), the file mass-N.json (mode table)
Needs:  SageMath (with GAP for test). The soft limit of the address space is set to 6 GB (RLIMIT_AS).
        Time: test about 10 seconds; table 32 about 1 to 3 minutes.

階数 n の全ての型 p-(z, d)-s について、F と Z の種の質量と、対 (L, sigma) の質量を求める。

論文：第 1 の論文の計算 C2（表 11）。奇素数位数の元の型の表。式は rank32.genus_mass にある。
モード：
  test        式を確かめる：p = 3 で [HKN09, Theorem 1] の値。階数 1 と 2 の mass_U を Z[zeta_p]^z の
              1/|Aut| と。p = 3、z = 6, 7 を Z[omega]^6 + K_12 と Z[omega]^7 + Z[omega] + K_12 と。
              E8 で、型についての M の和を W(E8) の位数 p の元の数と（GAP、Sage 経由）。階数 16 で
              E8 + E8 と D16^+ と。階数 24、p = 13 で A_12^2 と Leech 格子と。
  table N     階数 N の型 p-(z, d)-s（p <= 31）のうち種があるものを 1 行ずつ：密度の篩（最小 4）を
              通るか、mass_F、mass_H、|O|、M（分数と浮動小数）。
出力：--out DIR（既定は results/C2）の mass-N.json（モード table）
必要なもの：SageMath（test では GAP も）。アドレス空間のソフトリミットを 6 GB にする（RLIMIT_AS）。
      時間：test は約 10 秒、table 32 は約 1〜3 分。
"""
import argparse
import resource
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # code/ as the import root / code/ を import の起点にする

_lim_soft, _lim_hard = resource.getrlimit(resource.RLIMIT_AS)
resource.setrlimit(resource.RLIMIT_AS, (6 << 30, _lim_hard))  # soft limit only (callers may raise it) / ソフトリミットだけ（呼び出し側が上げられるように）

from sage.all import ZZ, QQ  # noqa: E402

from rank32 import genus_mass as gm  # noqa: E402
from rank32.io import write_json_atomic_indent0  # noqa: E402

DEFAULT_OUT = "results/C2"


def test():
    """Check the mass formulas against known values (see the module docstring).

    質量の式を既知の値と照合する（冒頭の説明を参照）。"""
    # p = 3: the values of [HKN09, Theorem 1] (Eisenstein lattices = (z = 4k, s = 0))
    # p = 3：[HKN09, Theorem 1] の値（Eisenstein 格子 = (z = 4k, s = 0)）
    vals = {4: QQ(1) / (2 ** 7 * 3 ** 5 * 5), 8: QQ(1) / (2 ** 15 * 3 ** 10 * 5 ** 2),
            12: QQ(691 * 809 * 1847) / (2 ** 22 * 3 ** 17 * 5 ** 3 * 7 * 11 * 13),
            16: QQ(13 * 47 * 419 * 691 * 809 * 1847 * 3617 * 16519) / (2 ** 31 * 3 ** 22 * 5 ** 4 * 11 * 17)}
    for z, v in vals.items():
        got = gm.mass_H(3, z, 0)
        print(f"HKN09 mu_{z}: {got} {'agrees' if got == v else 'differs ' + str(v)}")
        assert got == v
    # rank 1, 2: Z[zeta_p] only (Aut = 2p and (2p)^2 * 2); p = 3, 5, 7, 11, 13, 17, 19 have class number 1
    # 階数 1, 2 は Z[zeta_p] のみ（Aut = 2p と (2p)^2 * 2）。p = 3, 5, 7, 11, 13, 17, 19 は類数 1
    for p in (3, 5, 7, 11, 13, 17, 19):
        assert gm.mass_U(p, 1) == QQ(1) / (2 * p), (p, gm.mass_U(p, 1))
    for p in (3, 5):
        assert gm.mass_U(p, 2) == QQ(1) / (8 * p * p), (p, gm.mass_U(p, 2))
    print("rank 1 (p <= 19), rank 2 (p = 3, 5): agrees with 1/|Aut| of Z[zeta_p]^z")
    # p = 3, z = 6, 7: Feit. Z[omega]^6, K_12 / Z[omega]^7, Z[omega] + K_12
    # p = 3, z = 6, 7：Feit。Z[omega]^6, K_12 / Z[omega]^7, Z[omega] + K_12
    K12 = QQ(1) / 39191040
    assert gm.mass_U(3, 6) == QQ(1) / (6 ** 6 * 720) + K12
    assert gm.mass_U(3, 7) == QQ(1) / (6 ** 7 * 5040) + QQ(1) / 6 * K12
    print("p = 3, z = 6, 7: Z[omega]^6 + K_12 and Z[omega]^7 + Z[omega] (+) K_12 fill the mass exactly")
    # E8 (n = 8): orders of centralizers (type A_{p-1}^k for each k); the actual number is the number of elements of W(E8), |W|/|C|.
    # Type A2: 1120 subsystems x 2 elements = 2240 elements. Type A4^2: |C| = 600.
    # E8（n = 8）：中心化群の位数（各 k で A_{p-1}^k の型）。実際の数は W(E8) 内の元の数 = |W|/|C|。
    # A2 型：部分系 1120 個 x 2 元 = 2240 元。A4^2 型：|C| = 600。
    assert gm.M_pair(3, 8, 1, 1) == QQ(2240) / 696729600
    assert gm.M_pair(5, 8, 2, 0) == QQ(1) / 600
    print("E8: 3-(1,6)-1 has 2240 elements, 5-(2,0)-0 has |C| = 600: agrees")
    # the numbers of elements of order 3, 5, 7 of W(E8) (the sum of M over the types = number of elements/|W|)
    # E8 の位数 3, 5, 7 の元の総数（型ごとの M の和 = 元の数/|W|）
    W = ZZ(696729600)
    for p, tot in ((3, sum_elements_E8(3)), (5, sum_elements_E8(5)), (7, sum_elements_E8(7))):
        s = sum(gm.M_pair(p, 8, z, s_) for z in range(1, 8 // (p - 1) + 1) for s_ in range(0, z + 1))
        print(f"E8, p = {p}: sum over types of M = {s}, elements/|W| = {tot / W}, {'agrees' if s == tot / W else 'differs'}")
        assert s == tot / W
    # n = 16: E8^2 and D16^+ / n = 16：E8^2 と D16^+
    for p in (3, 5, 7):
        check_16(p)
    # n = 24, p = 13, type 13-(2,0)-0: the Niemeier lattice A_12^2 (|Aut| = (13!)^2 * 4) and the Leech lattice
    # (order 312 of the centralizer of an element of order 13 in Co_0 = 2.Co_1: the class 13a of the character table
    # "2.Co1" of the GAP character table library). This value is used only in this self-test.
    # n = 24、p = 13、型 13-(2,0)-0：Niemeier 格子 A_12^2（|Aut| = (13!)^2 * 4）と Leech 格子（Co_0 の位数 13 の
    # 元の中心化群の位数 312。GAP の指標表のライブラリの指標表 "2.Co1" の類 13a）。この値はこの自己試験でだけ使う。
    from sage.all import factorial
    v = QQ(factorial(12) ** 2) / (factorial(13) ** 2 * 4) + QQ(1) / 312
    assert gm.M_pair(13, 24, 2, 0) == v, (gm.M_pair(13, 24, 2, 0), v)
    print(f"n = 24, p = 13: M(13-(2,0)-0) = {v} agrees with the elements of A_12^2 and Leech")
    print("all passed")


def sum_elements_E8(p):
    """The number of elements of order p of W(E8).

    W(E8) の位数 p の元の数。"""
    return ZZ(_e8_order_counts()[p])


_e8_counts = None


def _e8_order_counts():
    global _e8_counts
    if _e8_counts is None:
        from sage.all import WeylGroup, gap  # noqa: F401
        G = WeylGroup(["E", 8], implementation="permutation")
        cc = G.gap().ConjugacyClasses()
        counts = {}
        for c in cc:
            o = int(c.Representative().Order())
            counts[o] = counts.get(o, 0) + int(c.Size())
        _e8_counts = counts
    return _e8_counts


def check_16(p):
    """n = 16: sum over types of M(p; z, s) = #{elements of order p of Aut E8^2}/|Aut E8^2| + #{... of Aut D16^+}/|Aut D16^+|.

    n = 16：sum_型 M(p; z, s) = #{Aut E8^2 の位数 p の元}/|Aut E8^2| + #{Aut D16^+ の位数 p の元}/|Aut D16^+|。"""
    from sage.all import factorial
    counts = _e8_order_counts()
    wE8 = ZZ(696729600)
    # Aut(E8^2) = W(E8) wr 2; the elements of odd prime order p lie in W(E8)^2: (1 + c)^2 - 1
    # Aut(E8^2) = W(E8) wr 2、奇素数 p の元は W(E8)^2 の中：(1 + c)^2 - 1
    c = ZZ(counts[p])
    e = ((1 + c) ** 2 - 1) / (2 * wE8 ** 2)
    # Aut(D16^+) = W(D16) (order 2^15 16!): the elements of odd prime order p are counted by k p-cycles and the signs.
    # The signs: only those with product +1 on each cycle (product -1 gives order 2p): 2^(p-1) per cycle; fixed points have sign +.
    # Aut(D16^+) = W(D16)（位数 2^15 16!）：奇素数 p の元は、k 個の p 巡回置換と符号の取り方で数える。
    # 符号は巡回置換ごとに積が +1 のものだけ（積 -1 だと位数 2p）：各巡回 2^(p-1) 通り、固定点は符号 + のみ。
    a = QQ(0)
    for k in range(1, 16 // p + 1):
        nperm = factorial(16) / (factorial(16 - k * p) * ZZ(p) ** k * factorial(k))
        a += nperm * ZZ(2) ** ((p - 1) * k)
    d = a / (ZZ(2) ** 15 * factorial(16))
    tot = sum(gm.M_pair(p, 16, z, s) for z in range(1, 16 // (p - 1) + 1) for s in range(0, z + 1))
    ok = tot == e + d
    print(f"n = 16, p = {p}: sum over types of M = {tot}, E8^2 + D16^+ = {e + d}, {'agrees' if ok else 'differs'}")
    assert ok


def table(n, out):
    """Write mass-n.json in the directory out: one row per type whose genera exist.

    種のある型ごとに 1 行の mass-n.json を、ディレクトリ out に書く。"""
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        for z in range(1, n // (p - 1) + 1):
            d = n - z * (p - 1)
            for s in range(0, min(z, d) + 1):
                if (z - s) % 2:
                    continue
                n2 = z * (p - 1)
                gex = gm.genus_exists(p, n2, s) and (d == 0 and s == 0 or d > 0 and gm.genus_exists(p, d, s))
                if not gex:
                    continue
                dens = (gm.sieve.ok(d, p, s, 4) if d else True) and gm.sieve.ok(n2, p, s, 4)
                M = gm.M_pair(p, n, z, s)
                mF = gm.mass_Z_genus(p, d, s) if d else QQ(1)
                mH = gm.mass_H(p, z, s)
                rec = dict(p=p, z=z, d=d, s=s, density_ok=bool(dens), mass_F=str(mF), mass_H=str(mH),
                           O=str(gm.order_O(p, s, gm.q_type(p, d, s))), M=str(M), M_float=float(M), mass_F_float=float(mF), mass_H_float=float(mH))
                rows.append(rec)
                if dens:
                    print(f"{p}-({z},{d})-{s}: mass of F {float(mF):.3e}, mass of Z {float(mH):.3e}, |O| {float(gm.order_O(p, s, gm.q_type(p, d, s))):.3e}, M {float(M):.3e}", flush=True)
    write_json_atomic_indent0(out / f"mass-{n}.json", rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mode", nargs="?", default="test", choices=["test", "table"], help="test or table / test か table")
    ap.add_argument("n", nargs="?", type=int, help="the rank (mode table) / 階数（モード table）")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output directory / 出力のディレクトリ")
    a = ap.parse_args()
    if a.mode == "test":
        test()
    else:
        table(a.n, Path(a.out))


if __name__ == "__main__":
    main()
