"""Calls to PARI qfauto / qfisom that avoid stack overflow and stalls, and a fast short-vector list.

Why the wrapper:
  - The default call searches for the depth of its pruning table automatically, and that is
    where the stack overflows and the memory grows. We pass fl = [1, 0] (depth fixed at 1).
    This only changes the strength of the pruning, never the answer.
    fl = [0, 0] turns the pruning off and is not used (a lattice of the 16-dimensional
    genus (8, 1) goes from 0.06 seconds to more than 30 seconds).
  - The search tree depends strongly on the basis (1 second and more than 150 seconds for the
    same lattice). We stop after a short time limit, change the basis at random and retry with
    a longer limit.
  - One long vector in the basis makes PARI enumerate everything up to its norm, so a new basis
    never increases the largest diagonal entry.
Not usable for: the extremal lattices of rank 32 (they are 7-designs and the pruning does not
work; use rank32.bliss_automorphisms), lattices with a huge group (A_1^24 and the like; build
the group from the root system), and some pairs of forms over the Eisenstein integers.
When the retries run out, QfTimeout is raised; callers record it as undecided.

Input: an integer matrix of Sage, or a list of matrices [G, F, ...] (the first positive definite;
qfauto then finds the maps that keep all of them). Conventions as in PARI: a generator g satisfies
g^T G g = G, and qfisom returns S with A = S^T B S. All returned matrices are in the original
basis and are checked before they are returned.
Uses the alarm of cysignals, so it cannot be combined with an alarm of the caller.

PARI の qfauto / qfisom を、スタック溢れと無応答を避けて呼ぶ。短いベクトルの速い一覧も置く。

包む理由：
  - 既定の呼び方は、枝刈りの表の深さを自動で探す。ここでスタックが溢れ、メモリが膨らむ。
    そこで fl = [1, 0]（深さを 1 に固定）を渡す。変わるのは枝刈りの強さだけで、答えは変わらない。
    fl = [0, 0] は枝刈りを切るので使わない（16 次元の種 (8, 1) の格子で 0.06 秒が 30 秒超になる）。
  - 探索のツリーは基底で桁違いに変わる（同じ格子で 1 秒と 150 秒超）。短い上限で切り、基底を乱択で
    替えて、上限を伸ばしながら再試行する。
  - 基底に長いベクトルが 1 本あると、そのノルムまで全部列挙する。基底を替えるときも、対角の
    最大値を増やさない。
使えない場合：階数 32 の極値格子（7-デザインで枝刈りが効かない。rank32.bliss_automorphisms を
使う）、群が巨大な格子（A_1^24 など。ルート系から群を作る）、Eisenstein 整数上の 2 形式の一部。
再試行を使い切ると QfTimeout を送出する。呼び出し側は「未判定」として記録する。

入力：Sage の整数行列、または行列の列 [G, F, ...]（先頭が正定値。qfauto はすべてを保つ変換を
求める）。約束は PARI と同じ：生成元 g は g^T G g = G、qfisom の S は A = S^T B S。返す行列は
すべて元の基底についてのもので、返す前に検算する。
cysignals の alarm を使うので、呼び出し側の alarm とは併用できない。
"""
import gc
import random
import warnings

import numpy as np
from cysignals.alarm import AlarmInterrupt, alarm, cancel_alarm
from sage.all import ZZ, identity_matrix, matrix, pari

FL = [1, 0]


class QfTimeout(TimeoutError):
    """The retries ran out.

    再試行の回数を使い切った。"""


def _forms(G):
    return list(G) if isinstance(G, (list, tuple)) else [G]


def _pari(forms):
    return pari(forms[0]) if len(forms) == 1 else pari(list(forms))


def _apply(forms, T):
    return [T.transpose() * F * T for F in forms]


def _maxdiag(G):
    return max(G[i, i] for i in range(G.nrows()))


def _random_basis(forms, rng):
    """A random unimodular change of basis T: a signed permutation followed by LLL; it never increases the largest diagonal entry of the first form.

    基底変換 T（ユニモジュラー）を乱択する。符号つき置換のあとに LLL。先頭の形式の対角の最大値は増やさない。"""
    G = forms[0]
    n = G.nrows()
    bound = _maxdiag(G)
    for _ in range(20):
        idx = list(range(n))
        rng.shuffle(idx)
        P = matrix(ZZ, n, n)
        for i, j in enumerate(idx):
            P[j, i] = rng.choice((1, -1))
        H = P.transpose() * G * P
        T = P * H.LLL_gram()
        if _maxdiag(T.transpose() * G * T) <= bound:
            return T
    return P  # keeps the diagonal; only the order of the branches changes / 対角を保つ。枝の順序が変わるだけで、悪化はしない


def _retry(f, forms, limit, tries, cap, seed):
    """Try f(forms in a changed basis) under a time limit; on timeout change the basis and multiply the limit by 1.5. Returns (result, T).

    cysignals keeps the interrupt exception until the next interrupt (cysigs.exc_value) and
    releases it inside sig_on at that time. If the traceback of the exception still holds this
    frame (the local r is a result on the PARI stack), the release happens inside sig_on and
    cypari2 raises SystemError('calling remove_from_pari_stack() inside sig_on()'). So we drop
    the traceback and collect cyclic garbage before the next try. A SystemError that still
    occurs counts as one try and is retried; if it repeats, QfTimeout is raised.

    f(基底を替えた forms) を上限つきで試し、切れたら基底を替えて上限を 1.5 倍にする。(結果, T) を返す。

    cysignals は割り込みの例外を次の割り込みまで握り（cysigs.exc_value）、次の割り込みのとき
    sig_on の中で手放す。例外のトレースバックがこの関数のフレーム（局所変数 r は PARI スタック上の
    結果）を抱えたままだと、その解放が sig_on の中で起き、cypari2 が
    SystemError('calling remove_from_pari_stack() inside sig_on()') を出す。そこでトレースバックを
    外し、次の試行の前に循環ゴミも回収する。それでも出た SystemError は 1 回の試行として数えて
    再試行し、続けば QfTimeout にする。"""
    rng = random.Random(seed)
    T = identity_matrix(ZZ, forms[0].nrows())
    errors = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="cypari2 leaked")  # the stack left after an interrupt returns at the next completion / 中断の後に残るスタックは次の完了で戻る
        for i in range(tries):
            if i:
                gc.collect()  # free cyclic garbage holding PARI Gen objects outside sig_on / PARI の Gen を含む循環ゴミを sig_on の外で解放する
            try:
                alarm(limit)
                r = f(_apply(forms, T))
                cancel_alarm()
                return r, T
            except AlarmInterrupt as e:
                e.__traceback__ = None
                limit = min(limit * 1.5, cap)
                T = _random_basis(forms, rng)
            except SystemError as e:
                e.__traceback__ = None
                errors.append(str(e))
                if len(errors) >= 3:
                    raise QfTimeout(f"SystemError of cypari2 repeated {len(errors)} times: {errors[-1]}") from None
            finally:
                cancel_alarm()
    raise QfTimeout(f"did not finish in {tries} tries (last limit {limit:.0f} seconds, SystemError {len(errors)} times)")


def qfauto(G, limit=4, tries=30, cap=600, seed=0):
    """Return (order, list of generators); each generator g satisfies g^T F g = F in the original basis for every form F.

    (位数, 生成元のリスト) を返す。生成元は元の基底で g^T F g = F（すべての形式 F）。"""
    forms = _forms(G)
    r, T = _retry(lambda fs: _pari(fs).qfauto(FL), forms, limit, tries, cap, seed)
    Ti = T.inverse().change_ring(ZZ)
    gens = [T * matrix(ZZ, g) * Ti for g in r[1]]
    for g in gens:
        if any(g.transpose() * F * g != F for F in forms):
            raise ArithmeticError("check of the generators of qfauto failed")
    return int(r[0]), gens


def qfisom(A, B, limit=4, tries=30, cap=600, seed=0):
    """Return an integer matrix S with A = S^T B S, or None if A and B are not isometric; only the basis of A is changed (it fixes the order of the search).

    A = S^T B S となる S（整数行列）を返す。同型でなければ None。基底を替えるのは A だけ（探索の順序は A で決まる）。"""
    fa, fb = _forms(A), _forms(B)
    r, T = _retry(lambda fs: pari.qfisom(_pari(fs), _pari(fb), FL), fa, limit, tries, cap, seed)
    if r == 0:
        return None
    S = matrix(ZZ, r) * T.inverse().change_ring(ZZ)
    if any(S.transpose() * Fb * S != Fa for Fa, Fb in zip(fa, fb)):
        raise ArithmeticError("check of the transformation matrix of qfisom failed")
    return S


def short_vectors(G, bound):
    """The vectors of norm <= bound, one of each pair +-v, as rows of an int64 array.

    matrix(ZZ, qfminim(...)[2]) builds one Python object per entry (21 seconds and 524 MB for
    70000 vectors of dimension 32); going through the string takes 1 second and 11 MB.
    qfminim is called with flag = 0 (enough for Gram matrices with small integer entries; a large
    rounding error raises an error; flag = 2 is about 30 times slower).

    ノルム <= bound のベクトルを、+-v の片方ずつ、行ベクトルの int64 配列で返す。

    matrix(ZZ, qfminim(...)[2]) は 1 要素ごとに Python のオブジェクトを作る（32 次元の 7 万本で
    21 秒・524 MB）。文字列を経由すれば 1 秒・11 MB である。
    qfminim は flag = 0 で呼ぶ（小さい整数の Gram 行列はこれでよく、丸め誤差が大きければエラーに
    なる。flag = 2 は約 30 倍遅い）。"""
    n = G.nrows()
    r = pari(G).qfminim(bound, None, 0)
    if int(r[0]) == 0:
        return np.zeros((0, n), dtype=np.int64)
    s = str(r[2].mattranspose())
    if s.startswith("Mat("):  # with a single row PARI prints Mat([a, b, ...]) / 1 行だけのとき PARI は Mat([a, b, ...]) と表示する
        s = s[4:-1]
    V = np.fromstring(s[1:-1].replace(";", ","), sep=",", dtype=np.int64).reshape(-1, n)
    if 2 * len(V) != int(r[0]):
        raise ArithmeticError("the number of vectors from qfminim does not match")
    return V
