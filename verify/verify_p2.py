"""Check the data of the paper on the automorphism groups of the 99 lattices with the largest groups (P2).

Usage, from the repository root:
    python3 verify/verify_p2.py [report.json]

Needs only the Python standard library and the files of this repository.  Exit status 0 means that every
check passed.  All arithmetic is exact integer arithmetic.

What is checked (data/largest-groups.jsonl, one line per lattice; field definitions in data/README.md):
  A  the 99 labels: distinct, of the form 32.<|Aut L|>.<i>, in the order (|Aut L|, i), and each with the same
     |Aut L| as the row of data/classified-lattices.jsonl (the data of [P1]) with the same label;
  B  each Gram matrix G: integer, symmetric, even, det G = 1 (Bareiss) and positive definite (leading minors);
  C  each generator g of Aut L: an integer matrix with g^T G g = G and det g = +-1;
  D  each presentation (P2, §2.3), on its integer matrices S (rows: S G S^T = G):
       - the matrices of the presentation preserve G;
       - every relation required by the lemma is present, its right-hand side is a word in the later stages,
         and it holds when both sides are evaluated as integer matrices;
       - the transpose of every generator of Aut L equals the product of its stored word (so the matrices of
         the presentation generate Aut L);
       - the bound of the lemma equals |Aut L|: the product of the relative orders (Lemma 2.1), |R| times the
         order of the presentation of Q = G/R (Lemma 2.3, radical and quotient), or the product of r_t and
         |T_t|^(k_t) over the stages (Lemma 2.3);
  E  the columns of Table A1 in the data equal paper/p2-lattices-99.json; |Aut L| equals its factorization in the
     table, and the lemma in the table is the lemma of the presentation;
  F  the profiles: for every conjugacy class, size * centralizer = |Aut L|, and its frame shape has degree 32,
     multiplicity of the eigenvalue 1 equal to the rank of the fixed lattice and a_1 equal to the trace; the class
     sizes add up to |Aut L|; the orbits on the minimal vectors add up to 146,880.
Not recomputed (the assumptions of P2, §2.5): the orders |Aut L| themselves; the orders of the finitely presented
groups of the block stages (the quotient Q and the simple sections T), which were counted by coset enumeration,
directly or along the chain of Lemma 2.4.  These orders are inputs; the checks print them.  For the chains of
Lemma 2.4 the stored subsets, relator counts, indices and last order are checked for consistency with each other,
with the stored relators and with |T|; the indices themselves are inputs.
"""
import hashlib
import json
import re
import sys
from fractions import Fraction
from operator import mul
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/largest-groups.jsonl'
CLASSIFIED = ROOT / 'data/classified-lattices.jsonl'
TABLE = ROOT / 'paper/p2-lattices-99.json'
TABLE_COLUMNS = ['aut_factor', 'fstar', 'g_mod_fstar', 'chief', 'presentation', 'origin']
N_LATTICES = 99
RANK = 32


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# ------------------------------------------------------------------ integer matrices (tuples of tuples)

def mat(a):
    return tuple(tuple(int(x) for x in row) for row in a)


def is_int_matrix(a, n=RANK):
    return (isinstance(a, list) and len(a) == n and all(isinstance(r, list) and len(r) == n for r in a)
            and all(type(x) is int for r in a for x in r))


def matmul(a, b):
    bt = tuple(zip(*b))
    return tuple(tuple(sum(map(mul, row, col)) for col in bt) for row in a)


def transpose(a):
    return tuple(zip(*a))


def identity(n=RANK):
    return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))


def det(a):
    """Exact determinant by fraction-free elimination (Bareiss)."""
    m = [list(r) for r in a]
    n, sign, prev = len(m), 1, 1
    for k in range(n - 1):
        if m[k][k] == 0:
            for i in range(k + 1, n):
                if m[i][k]:
                    m[k], m[i] = m[i], m[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
        prev = m[k][k]
    return sign * m[n - 1][n - 1]


def leading_minors(a):
    """The leading principal minors of a, by Bareiss elimination without row exchanges (stops at a zero minor)."""
    m = [list(r) for r in a]
    n, prev, minors = len(m), 1, []
    for k in range(n):
        minors.append(m[k][k])
        if m[k][k] == 0 or k == n - 1:
            break
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
        prev = m[k][k]
    return minors


def inverse_integer(a):
    """The inverse of a unimodular integer matrix, by Gauss-Jordan over the rationals; None if not integral."""
    n = len(a)
    m = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        p = next(i for i in range(k, n) if m[i][k] != 0)
        m[k], m[p] = m[p], m[k]
        piv = m[k][k]
        m[k] = [x / piv for x in m[k]]
        for i in range(n):
            if i != k and m[i][k] != 0:
                f = m[i][k]
                m[i] = [x - f * y for x, y in zip(m[i], m[k])]
    inv = [row[n:] for row in m]
    if any(x.denominator != 1 for row in inv for x in row):
        return None
    return tuple(tuple(int(x) for x in row) for row in inv)


class Group:
    """Integer matrices S of a presentation (rows: S G S^T = G), with cached powers and inverses.
    The inverse of S is G S^T G^{-1}, which is checked to be the inverse."""

    def __init__(self, G, Ginv, mats):
        self.G, self.Ginv = G, Ginv
        self.m = [mat(S) for S in mats]
        self.one = identity(len(G))
        self.preserves = all(matmul(matmul(S, G), transpose(S)) == G for S in self.m)
        self.inv = []
        self.inverse_ok = True
        for S in self.m:
            T = matmul(matmul(G, transpose(S)), Ginv)
            self.inverse_ok &= matmul(S, T) == self.one
            self.inv.append(T)
        self.cache = {}

    def power(self, a, e):
        """S_a^e for an integer e (negative: a power of the inverse)."""
        if e == 0:
            return self.one
        key = (a, e)
        if key not in self.cache:
            base = self.m[a] if e > 0 else self.inv[a]
            r = base
            for _ in range(abs(e) - 1):
                r = matmul(r, base)
            self.cache[key] = r
        return self.cache[key]

    def pairs(self, word):
        """A word given as a list of [generator, exponent] pairs, evaluated from left to right."""
        W = self.one
        for a, e in word:
            if e:
                W = matmul(W, self.power(a, e))
        return W

    def exponents(self, e, offset=0):
        """A word given by an exponent vector: S_{offset}^{e_0} S_{offset+1}^{e_1} ..."""
        return self.pairs([(offset + k, x) for k, x in enumerate(e) if x])

    def letters(self, rep, offset=0):
        """A word in GAP's letter representation: a > 0 is generator a, a < 0 its inverse (1-based, from offset)."""
        return self.pairs([(offset + abs(x) - 1, 1 if x > 0 else -1) for x in rep])


class Checks:
    def __init__(self):
        self.results = []

    def __call__(self, section, name, ok, detail=''):
        self.results.append({'section': section, 'check': name, 'ok': bool(ok), 'detail': str(detail)})
        print('PASS' if ok else 'FAIL', section, name, detail)


def factorization(n):
    f, p = [], 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            f.append((p, e))
        p += 1
    if n > 1:
        f.append((n, 1))
    return ' '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in f)


# ------------------------------------------------------------------ the three kinds of presentation

def check_pc(L, G, Ginv, gens, pr):
    """Lemma 2.1: x_i^{r_i} = w_ii and x_i^{-1} x_j x_i = w_ij (i < j), w in x_{i+1}, ..., x_m."""
    H = Group(G, Ginv, pr['generators'])
    m, ro = len(H.m), pr['relative_orders']
    shape = len(ro) == m and all(type(r) is int and r >= 1 for r in ro)
    want = {('power', i, i) for i in range(m)} | {('conjugate', i, j) for i in range(m) for j in range(i + 1, m)}
    seen, holds = set(), True
    for rel in pr['relations']:
        i = rel['i']
        j = rel.get('j', i)
        w = rel['word']
        seen.add((rel['kind'], i, j))
        shape &= len(w) == m and not any(w[:i + 1])
        if rel['kind'] == 'power':
            shape &= rel['r'] == ro[i]
            holds &= H.power(i, ro[i]) == H.exponents(w)
        else:
            holds &= matmul(matmul(H.inv[i], H.m[j]), H.m[i]) == H.exponents(w)
    shape &= seen == want and len(pr['relations']) == len(want)
    words = len(pr['generator_words']) == len(gens) and all(len(e) == m and transpose(g) == H.exponents(e)
                                                           for g, e in zip(gens, pr['generator_words']))
    bound = 1
    for r in ro:
        bound *= r
    return dict(preserve=H.preserves and H.inverse_ok, shape=shape, holds=holds, words=words, bound=bound,
                n_relations=len(pr['relations']), n_generators=m, inputs=[])


def check_radical(L, G, Ginv, gens, pr):
    """Lemma 2.3 with one block stage (k = 1) for Q = G/R, generated by the lifts y, followed by the cyclic stages
    of a polycyclic presentation of the solvable radical R: (i) the relations of Lemma 2.1 for R, (ii) y^{-1} x_j y =
    a word in x, (iii) each relator of the presentation of Q, evaluated on the lifts, = a word in x."""
    X = Group(G, Ginv, pr['radical_generators'])
    Y = Group(G, Ginv, pr['lifts'])
    m, ro, nl = len(X.m), pr['radical_relative_orders'], len(Y.m)
    shape = len(ro) == m and all(type(r) is int and r >= 1 for r in ro)
    want = {('power', i, i) for i in range(m)} | {('conjugate', i, j) for i in range(m) for j in range(i + 1, m)}
    seen, holds = set(), True
    for rel in pr['radical_relations']:
        i = rel['i']
        j = rel.get('j', i)
        w = rel['word']
        seen.add((rel['kind'], i, j))
        shape &= len(w) == m and not any(w[:i + 1])
        if rel['kind'] == 'power':
            shape &= rel['r'] == ro[i]
            holds &= X.power(i, ro[i]) == X.exponents(w)
        else:
            holds &= matmul(matmul(X.inv[i], X.m[j]), X.m[i]) == X.exponents(w)
    shape &= seen == want and len(pr['radical_relations']) == len(want)
    seen = set()
    for rel in pr['lift_conjugates']:
        a, j, w = rel['lift'], rel['j'], rel['word']
        seen.add((a, j))
        shape &= len(w) == m
        holds &= matmul(matmul(Y.inv[a], X.m[j]), Y.m[a]) == X.exponents(w)
    shape &= seen == {(a, j) for a in range(nl) for j in range(m)} and len(pr['lift_conjugates']) == nl * m
    shape &= len(pr['quotient_relators']) == len(pr['quotient_relator_values'])
    for rep, w in zip(pr['quotient_relators'], pr['quotient_relator_values']):
        shape &= len(w) == m and all(1 <= abs(x) <= nl for x in rep)
        holds &= Y.letters(rep) == X.exponents(w)
    words = len(pr['generator_words']) == len(gens)
    for g, gw in zip(gens, pr['generator_words']):
        words &= len(gw['radical']) == m and all(1 <= abs(x) <= nl for x in gw['quotient'])
        words &= transpose(g) == matmul(X.exponents(gw['radical']), Y.letters(gw['quotient']))
    nR = 1
    for r in ro:
        nR *= r
    q = pr['quotient_presentation_order']
    nQ = q['order']
    shape &= nQ == pr['quotient_order']
    inputs = [f"|Q_fp| = {nQ} ({q['method']})"]
    if 'lemma_2_4' in q:
        b = q['lemma_2_4']
        # Lemma 2.4 with X_0 = {y}: R_0 must contain y^r, and |Q_fp| <= index * r
        shape &= b['relator'] in pr['quotient_relators'] and set(b['relator']) == {b['generator']} and b['index'] * len(b['relator']) == nQ
        inputs.append(f"[Q_fp : <y_{b['generator']}>] = {b['index']} (input), relator y_{b['generator']}^{len(b['relator'])} present")
    return dict(preserve=X.preserves and Y.preserves and X.inverse_ok and Y.inverse_ok, shape=shape, holds=holds, words=words,
                bound=nR * nQ, n_relations=len(pr['radical_relations']) + len(pr['lift_conjugates']) + len(pr['quotient_relators']),
                n_generators=m + nl, inputs=inputs, radical_order=nR)


def check_subpresentation_chain(comp, T_order):
    """The stored chain of Lemma 2.4 for one simple section: consistent with itself, with the stored relators and with |T|."""
    rels = comp['relators']
    chain = comp['presentation_order']['lemma_2_4_chain']
    ok = chain[0]['generators'] == list(range(1, len(comp['generators']) + 1))
    product = 1
    for s, nxt in zip(chain, chain[1:] + [None]):
        S = set(s['generators'])
        ok &= s['relators'] == sum(1 for r in rels if {abs(x) for x in r} <= S)
        if nxt is None:
            ok &= 'order' in s and 'index' not in s
            product *= s['order']
        else:
            ok &= set(s['subgroup']) < S and nxt['generators'] == s['subgroup'] and s['index'] >= 1
            product *= s['index']
    ok &= product == comp['presentation_order']['order'] == T_order
    indices = ' * '.join(str(s.get('index', s.get('order'))) for s in chain)
    return ok, indices


def check_chain(L, G, Ginv, gens, pr):
    """Lemma 2.3: cyclic stages (x^r = word in later stages) and block stages (each relator of each component and the
    commutators of generators of different components = words in later stages), and x^{-1} y x = word in the stages
    after that of x, for every generator x and every generator y of a later stage."""
    H = Group(G, Ginv, pr['generators'])
    n = len(H.m)
    stage_of, own = {}, []
    shape = True
    a = 0
    for t, st in enumerate(pr['stages']):
        gs = st['generators'] if st['kind'] == 'cyclic' else [x for c in st['components'] for x in c['generators']]
        shape &= gs == list(range(a, a + len(gs)))
        if st['kind'] == 'block':
            shape &= len(st['components']) == st['k'] and all(1 <= abs(x) <= len(c['generators'])
                                                               for c in st['components'] for r in c['relators'] for x in r)
        a += len(gs)
        own.append(gs)
        for x in gs:
            stage_of[x] = t
    shape &= a == n
    want = set()
    for t, st in enumerate(pr['stages']):
        later = range(own[t][-1] + 1, n)
        if st['kind'] == 'cyclic':
            want.add(('power', t, own[t][0], None))
        else:
            for j, c in enumerate(st['components']):
                want |= {('relator', t, j, q) for q in range(len(c['relators']))}
            for j, c in enumerate(st['components']):
                for d in st['components'][j + 1:]:
                    want |= {('commutator', t, x, y) for x in c['generators'] for y in d['generators']}
        want |= {('conjugate', t, x, y) for x in own[t] for y in later}
    seen, holds = set(), True
    for rel in pr['relations']:
        t, w, kind = rel['stage'], rel['word'], rel['kind']
        shape &= all(stage_of[x] > t for x, _ in w)
        if kind == 'power':
            x = rel['a']
            seen.add(('power', t, x, None))
            shape &= pr['stages'][t]['kind'] == 'cyclic' and rel['r'] == pr['stages'][t]['r'] and own[t] == [x]
            holds &= H.power(x, rel['r']) == H.pairs(w)
        elif kind == 'relator':
            j, q = rel['component'], rel['relator']
            seen.add(('relator', t, j, q))
            c = pr['stages'][t]['components'][j]
            holds &= H.letters(c['relators'][q], c['generators'][0]) == H.pairs(w)
        elif kind == 'commutator':
            x, y = rel['a'], rel['b']
            seen.add(('commutator', t, x, y))
            holds &= matmul(matmul(H.inv[x], H.inv[y]), matmul(H.m[x], H.m[y])) == H.pairs(w)
        else:
            x, y = rel['a'], rel['b']
            seen.add(('conjugate', t, x, y))
            shape &= stage_of[x] == t
            holds &= matmul(matmul(H.inv[x], H.m[y]), H.m[x]) == H.pairs(w)
    shape &= seen == want and len(pr['relations']) == len(want)
    words = len(pr['generator_words']) == len(gens) and all(transpose(g) == H.pairs(w) for g, w in zip(gens, pr['generator_words']))
    bound, inputs = 1, []
    for st in pr['stages']:
        if st['kind'] == 'cyclic':
            bound *= st['r']
            continue
        orders = [c['presentation_order']['order'] for c in st['components']]
        shape &= all(o == st['T_order'] for o in orders)
        bound *= st['T_order'] ** st['k']
        for c in st['components']:
            po = c['presentation_order']
            if 'lemma_2_4_chain' in po:
                ok, indices = check_subpresentation_chain(c, st['T_order'])
                shape &= ok
                inputs.append(f"|T_fp| = {po['order']} = {indices} (Lemma 2.4; indices input; {len(c['generators'])} generators, "
                              f"{len(c['relators'])} relators)")
            else:
                inputs.append(f"|T_fp| = {po['order']} ({po['method']})")
    shape &= bound == pr['order_bound']
    return dict(preserve=H.preserves and H.inverse_ok, shape=shape, holds=holds, words=words, bound=bound,
                n_relations=len(pr['relations']), n_generators=n, inputs=inputs)


CHECKERS = {'pc': check_pc, 'radical': check_radical, 'chain': check_chain}


def main():
    c = Checks()
    with open(DATA) as f:
        recs = [json.loads(line) for line in f]
    table = json.loads(TABLE.read_text())['lattices']
    with open(CLASSIFIED) as f:
        p1 = {}
        for line in f:
            r = json.loads(line)
            p1[r['label']] = r

    # A. labels
    labels = [r['label'] for r in recs]
    form = all(re.fullmatch(r'32\.([0-9]+)\.([1-9][0-9]*)', s) and int(s.split('.')[1]) == r['aut'] and type(r['aut']) is int
               for s, r in zip(labels, recs))
    c('A', f'{N_LATTICES} lattices with distinct labels 32.<|Aut L|>.<i>', len(recs) == N_LATTICES and len(set(labels)) == len(labels) and form,
      f'{len(recs)} lines')
    keys = [(r['aut'], int(r['label'].split('.')[2])) for r in recs]
    c('A', 'lines in the order (|Aut L|, i)', keys == sorted(keys))
    bad = [s for s, r in zip(labels, recs) if s not in p1 or int(p1[s]['aut']) != r['aut']]
    parts = sorted({p1[s]['part'] for s in labels if s in p1})
    c('A', 'each label is a row of data/classified-lattices.jsonl with the same |Aut L|', not bad,
      f'parts {parts}; failures {bad[:5]}')

    # B. Gram matrices
    bad, inverses = [], {}
    for r in recs:
        G = r['gram']
        ok = is_int_matrix(G)
        if ok:
            G = mat(G)
            ok = G == transpose(G) and all(G[i][i] % 2 == 0 for i in range(RANK)) and det(G) == 1
            ok = ok and all(x > 0 for x in leading_minors(G)) and len(leading_minors(G)) == RANK
            Ginv = inverse_integer(G) if ok else None
            ok = ok and Ginv is not None and matmul(G, Ginv) == identity()
            inverses[r['label']] = Ginv
        if not ok:
            bad.append(r['label'])
    c('B', 'Gram matrices: integer, symmetric, even, det 1, positive definite (32 leading minors > 0)', not bad,
      f'{len(recs)} matrices; failures {bad[:5]}')

    # C. generators of Aut L
    bad, n_gens = [], 0
    for r in recs:
        G = mat(r['gram'])
        ok = r['generators_convention'] == 'columns: x -> g x, g^T G g = G' and len(r['generators']) >= 1
        for g in r['generators']:
            n_gens += 1
            ok &= is_int_matrix(g)
            if ok:
                g = mat(g)
                ok &= matmul(matmul(transpose(g), G), g) == G and det(g) in (1, -1)
        if not ok:
            bad.append(r['label'])
    c('C', 'generators of Aut L: integer, g^T G g = G, det g = +-1', not bad, f'{n_gens} generators; failures {bad[:5]}')

    # D. presentations
    results, inputs = {}, []
    for r in recs:
        pr = r['presentation']
        kind = pr['kind']
        G, Ginv = mat(r['gram']), inverses.get(r['label'])
        if Ginv is None or kind not in CHECKERS:
            results[r['label']] = dict(kind=kind, preserve=False, shape=False, holds=False, words=False, bound=None, n_relations=0)
            continue
        gens = [mat(g) for g in r['generators']]
        try:
            res = CHECKERS[kind](r, G, Ginv, gens, pr)
        except (KeyError, IndexError, TypeError, ValueError, StopIteration) as e:
            # malformed data (a missing field, a generator index out of range): the lattice fails every check of D
            print(f"FAIL D {r['label']}: malformed presentation ({type(e).__name__}: {e})")
            res = dict(preserve=False, shape=False, holds=False, words=False, bound=None, n_relations=0, inputs=[])
        res['kind'] = kind
        res['shape'] &= pr['convention'].startswith('rows: v -> v S, S G S^T = G') and pr['lemma'] == {'pc': '2.1'}.get(kind, '2.3')
        results[r['label']] = res
        for x in res['inputs']:
            inputs.append(f"{r['label']}: {x}")
    for kind, lemma in (('pc', 'Lemma 2.1'), ('radical', 'Lemma 2.3, radical and quotient'), ('chain', 'Lemma 2.3, stages along a chief series')):
        sel = {s: x for s, x in results.items() if x['kind'] == kind}
        nrel = sum(x['n_relations'] for x in sel.values())
        for key, what in (('preserve', 'the matrices of the presentation preserve G (S G S^T = G)'),
                          ('shape', 'every required relation is present, with its right-hand side in the later stages'),
                          ('holds', 'every relation holds on the integer matrices'),
                          ('words', 'the transpose of every generator of Aut L is the product of its word')):
            bad = [s for s, x in sel.items() if not x[key]]
            c('D', f'{lemma}: {what}', sel and not bad, f'{len(sel)} lattices, {nrel} relations; failures {bad[:5]}')
        bad = [s for s, x in sel.items() if x['bound'] != recs[labels.index(s)]['aut']]
        c('D', f'{lemma}: the bound of the lemma equals |Aut L|', sel and not bad, f'{len(sel)} lattices; failures {bad[:5]}')
    by_kind = {k: sum(1 for x in results.values() if x['kind'] == k) for k in CHECKERS}
    c('D', 'every lattice has a presentation of one of the three kinds', sum(by_kind.values()) == len(recs) == N_LATTICES, by_kind)
    c('D', 'orders of finitely presented groups used as inputs (coset enumeration; P2 §2.5)', True, '; '.join(inputs))

    # E. Table A1
    same = len(table) == len(recs) and all(t['label'] == r['label'] and t['id'] == r['id'] and t['aut'] == r['aut']
                                           and all(t[k] == r['table_a1'][k] for k in TABLE_COLUMNS) for t, r in zip(table, recs))
    c('E', 'the columns of Table A1 in the data equal paper/p2-lattices-99.json, row by row', same, f'{len(table)} rows')
    bad = [r['label'] for r in recs if r['table_a1']['aut_factor'] != factorization(r['aut'])]
    c('E', 'Table A1: |Aut L| equals its factorization', not bad, f'failures {bad[:5]}')

    def lemma_column(r):
        pr = r['presentation']
        if pr['kind'] == 'radical':     # Lemma 2.4 bounds the quotient of 12-20-c63:236 (stored as lemma_2_4)
            return '2.3, 2.4' if 'lemma_2_4' in pr['quotient_presentation_order'] else '2.3'
        if pr['kind'] != 'chain':
            return pr['lemma']
        uses_2_4 = any('lemma_2_4_chain' in comp['presentation_order'] for st in pr['stages'] if st['kind'] == 'block' for comp in st['components'])
        return '2.3, 2.4' if uses_2_4 else '2.3'
    bad = [r['label'] for r in recs if r['table_a1']['presentation'] != lemma_column(r)]
    c('E', 'Table A1: the column "pres." is the lemma of the stored presentation', not bad, f'failures {bad[:5]}')

    # F. The profiles (§2.2): consistency of the stored action on L.  For each class: size * centralizer = |Aut L|; the
    # frame shape prod (x^k - 1)^(a_k) has degree sum k a_k = 32, its multiplicity of the eigenvalue 1, sum a_k, is the
    # rank of the fixed lattice, and a_1 is the trace.  The class sizes add up to |Aut L|, with one class of the
    # identity, and the orbits on the minimal vectors add up to 146,880.  The classes themselves (GAP) are inputs.
    def frame_exponents(frame):
        out = {}
        for tok in frame.split():
            k, _, a = tok.partition('^')
            out[int(k)] = out.get(int(k), 0) + (int(a) if a else 1)
        return out
    bad = []
    for r in recs:
        pf, aut = r['profile'], r['aut']
        cl = pf['action_on_L']['classes']
        ok = sum(x['size'] for x in cl) == aut and sum(1 for x in cl if x['order'] == 1) == 1
        for x in cl:
            fr = frame_exponents(x['frame'])
            ok &= x['size'] * x['centralizer'] == aut and sum(k * a for k, a in fr.items()) == RANK
            ok &= sum(fr.values()) == x['fixed']['rank'] and fr.get(1, 0) == x['trace']
        ok &= sum(pf['action_on_L']['orbits_on_minimal_vectors']) == 146880
        if not ok:
            bad.append(r['label'])
    c('F', 'profiles: class sizes add up to |Aut L| (size * centralizer = |Aut L|); frame shapes of degree 32, with fixed rank '
           'and trace; orbits on the minimal vectors add up to 146,880', not bad,
      f"{sum(len(r['profile']['action_on_L']['classes']) for r in recs)} classes; failures {bad[:5]}")

    failures = [x for x in c.results if not x['ok']]
    report = {'checks': c.results, 'failures': len(failures), 'verifier_sha256': sha256(Path(__file__)),
              'data_sha256': sha256(DATA)}
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(json.dumps(report, indent=1) + '\n')
    print(f'{len(c.results) - len(failures)} of {len(c.results)} checks passed')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
