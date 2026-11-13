"""Check the certificates and data of this repository against the paper.

Usage, from the repository root:
    python3 verify/verify.py [report.json]

Needs only the Python standard library and the files of this repository (no
release assets).  Exit status 0 means that every check passed.

What is checked (details in verify/README.md):
  A  digests of all files listed in MANIFEST.json;
  B  every computation cited in the paper is in its Table 11;
  C  the six rational infeasibility certificates of Appendix B.1.1;
  D  identity (5.1): King's mass and the classified lattices give an integer,
     the totals of Theorems A and B, and Table 8;
  E  the Gaussian rootless mass M_2: every embedding average A_16(T) recomputed
     (verify/check_m2_terms.py), every term of (C.2) from its coefficient,
     group order and embedding average, the coefficients from the closed
     formulas of Appendix C.3 and the unitary orders of Proposition C.1;
  F  integer certificates: generators of the two largest groups preserve their
     Gram matrices, and the 106 identifying isometries satisfy T^T G T = G'
     with det T = +-1;
  G  the classification counts of §6 against the classified lattices, and the
     row-level witnesses of Lemmas 6.3 and 6.6 (general-part-witnesses.json):
     the identification classes from the saved rows and edges, a witness row of
     an exact type for 12,519 classes, and the subtraction as set differences of
     individual list entries and Gaussian keys; the 849 quaternion triples in
     integers (quaternion-triples.json): Gram matrices even, positive definite
     and unimodular, W G W^T = G and W^2 = -1, and two isometries T^T G T = G'
     with det T = +-1 joining the three classes of each triple;
  H  the finite checks of the 9-torsion reduction (verify/check_exponent9.py) and of
     the orbits of elementary quadratic 2-forms (verify/check_quadratic2_orbits.py);
  I  the completeness identities from their terms: (4.1) for the 83 pairs of
     involution components and for the types of odd prime order, the mass
     identities of Table B4 and (B.2), identity (B.1) for type (8, 24), and
     sum f/|U| = M_2 over the Gaussian list (Proposition 4.6) and the counts of
     the methods of its weights, with the counts of Tables 5, 6, B3 and B4.
It does not recompute automorphism groups, lists, masses of genera or isometry
decisions; those rest on the statements and assumptions listed in Tables 7 and
10 of the paper.
"""
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'paper'
CERT = ROOT / 'certificates'


def load(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# Numbers in the Typst source of the papers: integers of 4 or more digits are grouped by commas, written 11,044,268 in
# text and "11,044,268" (a text in math) inside $...$.
def fm(n):
    return f'"{n:,}"' if abs(n) >= 1000 else str(n)


def ft(n):
    return f'{n:,}'


NUM = r'"?([0-9][0-9,]*)"?'


def fmf(x):
    # a decimal in math, 3 places: "34,313.567" when the integer part has 4 or more digits
    t = f'{float(x):,.3f}'
    return f'"{t}"' if abs(float(x)) >= 1000 else t


def num(s):
    return int(s.replace(',', ''))


def matmul(a, b):
    bt = list(zip(*b))
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def transpose(a):
    return [list(r) for r in zip(*a)]


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


def order_GO(n, q, e=None):
    """|GO^e_n(q)|, the full orthogonal group of a non-degenerate quadratic form (n odd only for odd q)."""
    if n == 0:
        return 1
    if n % 2:
        m = (n - 1) // 2
        r = 2 * q ** (m * m)
        for i in range(1, m + 1):
            r *= q ** (2 * i) - 1
        return r
    m = n // 2
    r = 2 * q ** (m * (m - 1)) * (q ** m - e)
    for i in range(1, m):
        r *= q ** (2 * i) - 1
    return r


def order_Sp(n, q):
    m = n // 2
    r = q ** (m * m)
    for i in range(1, m + 1):
        r *= q ** (2 * i) - 1
    return r


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def order_O_odd(p, d, s):
    """Lemma 2.1(i): |O(q)| for a p-elementary lattice of rank d with D = F_p^s, p odd."""
    if s % 2:
        return order_GO(s, p)
    # chi(disc Q) = e^{2 pi i d/8} eps_p^{-s}, eps_p = 1 or i; as an 8th root of unity, index d - 2s [p = 3 mod 4]
    k = (d - (2 * s if p % 4 == 3 else 0)) % 8
    if k not in (0, 4):
        return None
    chi_disc = 1 if k == 0 else -1
    return order_GO(s, p, legendre((-1) ** (s // 2), p) * chi_disc)


def quadratic_value_counts(p, d, s):
    """For a p-elementary lattice of rank d with D = F_p^s (p odd): the number N(j) of non-zero elements with
    Q(a) = j, j = 0, ..., p - 1, where Q = (p/2) q (Table 3); the type of Q is that of Lemma 2.1(i)."""
    k = (d - (2 * s if p % 4 == 3 else 0)) % 8
    if k not in (0, 4):
        return None
    chi_disc = 1 if k == 0 else -1
    counts = {}
    if s % 2 == 0:
        eps = legendre((-1) ** (s // 2), p) * chi_disc
        counts[0] = p ** (s - 1) + eps * (p ** (s // 2) - p ** (s // 2 - 1)) - 1
        for j in range(1, p):
            counts[j] = p ** (s - 1) - eps * p ** (s // 2 - 1)
    else:
        eta = legendre((-1) ** ((s - 1) // 2), p) * chi_disc
        counts[0] = p ** (s - 1) - 1
        for j in range(1, p):
            counts[j] = p ** (s - 1) + legendre(j, p) * eta * p ** ((s - 1) // 2)
    return counts


def lp_rows(case):
    """The inequalities of a linear program of Appendix B.1.1, rebuilt from the stored bases of theta series."""
    p, z, d, s = case['type']
    spF, spZ = case['spaces']
    F = [(list(map(Q, b['theta'])), list(map(Q, b['dual_scaled']))) for b in spF['basis']]
    Z = [(list(map(Q, b['theta'])), list(map(Q, b['dual_scaled']))) for b in spZ['basis']]
    rows = []

    def coeff(part, which, j):
        if part == 'F':
            return [v[which][j] for v in F] + [Q(0)] * len(Z)
        return [Q(0)] * len(F) + [v[which][j] for v in Z]

    def add(a, b, label, equal=False):
        rows.append({'a': [str(x) for x in a], 'b': str(b), 'label': label})
        if equal:
            add([-x for x in a], -b, label + ' (reverse)')

    def m_alpha(alpha):
        return 1 if alpha < 1 else 1 + int(alpha / (2 - alpha))

    for part in ['F', 'Z']:
        for which, j, val in [(0, 0, 1), (0, 1, 0), (1, 0, 1)]:
            add(coeff(part, which, j), val, f'{part}:{which}:{j}={val}', True)
        for j in range(2, 40):
            add([-x for x in coeff(part, 0, j)], 0, f'{part}:theta:{j}>=0')
        for j in range(1, 40):
            add([-x for x in coeff(part, 1, j)], 0, f'{part}:dual:{j}>=0', equal=part == 'Z' and j <= p - 2)
    counts = quadratic_value_counts(p, d, s)
    for j in range(1, p):
        mF, mZ = m_alpha(Q(2 * j, p)), m_alpha(Q(2 * (p - j), p))
        add([u / mF + v / mZ for u, v in zip(coeff('F', 1, j), coeff('Z', 1, p - j))], counts[j], f'gluing:{j}')
    return rows, counts


def order_O_two(d, s, delta):
    """Lemma 2.1(ii): |O(q)| for a 2-elementary lattice of rank d with invariants s, delta."""
    if delta == 0:
        return order_GO(s, 2, 1 if d % 8 == 0 else -1) if d % 4 == 0 else None
    if s % 2:
        return order_GO(s - 1, 2, 1 if d % 8 in (1, 7) else -1)
    if d % 4 == 2:
        return 2 * order_Sp(s - 2, 2)
    return 2 ** (s - 2) * order_GO(s - 2, 2, 1 if d % 8 == 0 else -1)


def inverse_fraction(a):
    """Exact inverse by Gauss-Jordan elimination over the rationals."""
    n = len(a)
    m = [[Q(x) for x in row] + [Q(int(i == j)) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        p = next(i for i in range(k, n) if m[i][k] != 0)
        m[k], m[p] = m[p], m[k]
        piv = m[k][k]
        m[k] = [x / piv for x in m[k]]
        for i in range(n):
            if i != k and m[i][k] != 0:
                f = m[i][k]
                m[i] = [x - f * y for x, y in zip(m[i], m[k])]
    return [row[n:] for row in m]


class Checks:
    def __init__(self):
        self.results = []

    def __call__(self, section, name, ok, detail=''):
        self.results.append({'section': section, 'check': name, 'ok': bool(ok), 'detail': str(detail)})
        print('PASS' if ok else 'FAIL', section, name, detail)


def main():
    c = Checks()
    paper_files = sorted(PAPER.glob('*.typ'))
    body = ''.join(p.read_text() for p in paper_files if p.name not in ('computations.typ', 'style.typ', 'assignment-count.typ'))
    main_text = (PAPER / 'extremal-lattices-rank-32.typ').read_text()
    ledger = (PAPER / 'computations.typ').read_text()

    # A. digests
    manifest = ROOT / 'MANIFEST.json'
    if manifest.exists():
        m = load(manifest)
        bad = [f['path'] for f in m['files'] if not (ROOT / f['path']).is_file() or sha256(ROOT / f['path']) != f['sha256']]
        c('A', 'files match MANIFEST.json', not bad, f"{len(m['files'])} files; mismatches {bad[:5]}")
    else:
        c('A', 'MANIFEST.json present', False)

    # B. ledger
    cited = set(re.findall(r'C[0-9]+', ' '.join(re.findall(r'#calc\(([^)]*)\)', body))))
    defined = set(re.findall(r'\(id: "(C[0-9]+)"', ledger))
    c('B', 'every cited computation is in Table 11', cited <= defined, f'missing {sorted(cited - defined)}')

    # C. LP certificates
    lp = load(CERT / 'type-exclusion-lp.json')
    ok = len(lp['cases']) == 6
    for case in lp['cases']:
        rows, y = case['inequalities'], [Q(v) for v in case['witness']]
        ok &= len(y) == len(rows) and all(v >= 0 for v in y)
        for j in range(case['variables']):
            ok &= sum(v * Q(r['a'][j]) for v, r in zip(y, rows)) == 0
        ok &= sum(v * Q(r['b']) for v, r in zip(y, rows)) == -1
    c('C', 'six rational certificates: y >= 0, yA = 0, yb = -1', ok, ', '.join(str(case['type']) for case in lp['cases']))
    # The rows themselves, rebuilt from the stored bases of theta series and dual theta series: the normalizations, the
    # non-negativity of the coefficients, and the gluing equations with the numbers N(j) of classes with q = 2j/p,
    # which are recomputed from (p, d, s) (Lemma 2.1(i)).  The bases (spaces of modular forms) are inputs.
    bad = []
    for case in lp['cases']:
        try:
            rows, counts = lp_rows(case)
        except (KeyError, TypeError, ZeroDivisionError):
            rows, counts = None, None
        p_ = case['type'][0]
        if rows != case['inequalities'] or counts is None or {str(j): str(counts[j]) for j in range(p_)} != case['class_counts']:
            bad.append(case['type'])
    c('C', 'the inequalities of the six systems rebuilt from the stored theta series, with the class counts N(j) from (p, d, s)',
      not bad, f'failures {bad}')

    # D. identity (5.1)
    rows = [json.loads(l) for l in open(ROOT / 'data/classified-lattices.jsonl')]
    m = Q(load(CERT / 'king-rootless-mass.json')['mass'])
    shown = re.search(r'\$ m = ' + NUM + ' / ' + NUM + ' approx', main_text)
    c('D', "King's mass in the paper equals the certificate", shown is not None and Q(num(shown[1]), num(shown[2])) == m)
    labels = [r['label'] for r in rows]
    c('D', 'labels distinct; orders even and > 2', len(set(labels)) == len(labels) and all(r['aut'] > 2 and r['aut'] % 2 == 0 for r in rows))
    total = 2 * m + sum((1 - Q(2, r['aut']) for r in rows), Q(0))
    c('D', 'N = 2m + sum (1 - 2/|Aut L|) is an integer', total.denominator == 1, total)
    for name, value in (('N', int(total)), ('|S|', len(rows)), ('|Aut L| = 2', int(total) - len(rows))):
        c('D', f'paper states {name}', f'${fm(value)}$' in main_text, fm(value))
    # Table 9 read from the paper: every (|Aut L|, number) cell equals the data, |Aut L| = 2 included.
    t9 = main_text[main_text.index('caption: [Table 9.'):]
    t9 = t9[:t9.index('\n)\n')]
    cells = re.findall(r'\[\$' + NUM + r'\$\], \[\$' + NUM + r'\$\]', t9)
    table9 = {num(a): num(b) for a, b in cells}
    data9 = Counter(r['aut'] for r in rows)
    data9[2] = int(total) - len(rows)
    c('D', 'Table 9 (lattices by |Aut L|) agrees with the data, cell by cell', len(cells) == len(table9) and table9 == dict(data9),
      f'{len(table9)} orders in the table, {len(data9)} in the data')
    parts = Counter(r['part'] for r in rows)
    c('D', 'parts general / c2xc2 / cyclic-quaternion', (parts['general'], parts['c2xc2'], parts['cyclic-quaternion']) == (12546, 67284, 63499), dict(parts))
    cq = Counter(r['group']['structure'] for r in rows if r['part'] == 'cyclic-quaternion')
    expected = {'C4': 62468, 'C8': 175, 'C16': 7, 'Q8': 793, 'Q16': 48, 'Q32': 6, 'Q64': 2}
    c('D', 'Table 8 (the lattices of S_CQ) agrees with the data',
      all(cq.get(k) == v and f'[${k[0]}_{k[1:]}$], [${fm(v)}$]' in main_text for k, v in expected.items()) and sum(cq.values()) == 63499, dict(cq))

    # E. M_2
    g = load(CERT / 'gaussian-rootless-mass.json')
    coeff = {k: Q(v) for k, v in g['irreducible_coefficients'].items()}
    aut_e = {6: 103680, 7: 2903040, 8: 696729600}

    def unitary(name, n):
        return {'A': lambda: 4 if n == 1 else 4 * math.factorial(n + 1),
                'D': lambda: 2304 if n == 4 else 2 ** (n + 1) * math.factorial(n),
                'E': lambda: 2 * aut_e[n],
                'G': lambda: 96 if n == 2 else 4 ** n * math.factorial(n),
                'Eg': lambda: 46080}[name]()

    ok_terms = True
    for t in g['terms']:
        cf, order = Q(1), 1
        for (name, n), mult in t['type']:
            cf *= coeff[f'{name}{n}'] ** mult
            order *= unitary(name, n) ** mult * math.factorial(mult)
        ok_terms &= cf == Q(t['coefficient']) and order == t['group_order'] and Q(t['term']) == cf * Q(t['embedding_average']) / order
    c('E', 'every term of (C.2): coefficient, |U(T)| and term', ok_terms, f"{len(g['terms'])} terms")
    value = sum((Q(t['term']) for t in g['terms']), Q(0))
    shown = re.search(r'M_2 = ' + NUM + ' slash ' + NUM, main_text)
    c('E', 'sum of the terms is M_2 of Theorem 4.5', value == Q(g['rootless_mass']) and shown is not None and value == Q(num(shown[1]), num(shown[2])), value)
    closed = {'A': lambda n: (-1) ** n * math.factorial(n), 'D': lambda n: (-1) ** n * math.factorial(n - 1) * (2 ** (n - 1) - 1),
              'G': lambda m: (-1) ** m * 2 ** (m - 1) * math.factorial(m - 1)}
    ok = all(Q(closed[k[0]](int(k[1:]))) == v for k, v in coeff.items() if k[0] in closed and k[1:].isdigit())
    c('E', 'coefficients of A_n, D_n, G_m from the closed formulas of Appendix C.3', ok)
    rc = load(CERT / 'root-coefficients.json')
    by_rank, nonzero = [0] * 17, [0] * 17
    for t in g['terms']:
        by_rank[t['rank']] += 1
        nonzero[t['rank']] += Q(t['embedding_average']) != 0
    c('E', 'root-coefficients.json: its counts of root types by rank agree with the term table',
      rc['counts_by_rank'] == by_rank and rc['nonzero_by_rank'] == nonzero and rc['total_types'] == rc['terms'] == len(g['terms'])
      and rc['terms_with_nonzero_average'] == sum(nonzero) and not rc['coefficient_mismatches'] and not rc['term_mismatches'],
      f'{len(g["terms"])} types, {sum(nonzero)} with non-zero average')
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, '-I', str(ROOT / 'verify/check_m2_terms.py'), str(Path(tmp) / 'm2.json')], capture_output=True, text=True)
        m2 = load(Path(tmp) / 'm2.json') if r.returncode == 0 else {}
    c('E', 'every A_16(T) of (C.2) recomputed from the root data and the local 2-adic statistics (verify/check_m2_terms.py); '
      'the sum is M_2', r.returncode == 0, (r.stdout + r.stderr)[-300:])

    # F. integer certificates
    ex = load(CERT / 'exceptional-full-groups.json')
    ok = True
    for lat in ex['lattices']:
        G = lat['gram']
        ok &= all(matmul(matmul(transpose(h), G), h) == G for h in lat['generators'])
        ok &= int(lat['order']) == 2 * lat['orbit_lines'] * int(lat['stabilizer_order'])
        ok &= any(r['label'] == lat['label'] and r['aut'] == int(lat['order']) for r in rows)
    c('F', 'generators of the two largest groups preserve the Gram matrix; orders agree with the data', ok)
    ident = load(CERT / 'general-part-identification.json')
    bad = []
    for iso in ident['isometries']:
        T, G1, G2 = iso['T'], iso['gram_first'], iso['gram_second']
        if matmul(matmul(transpose(T), G1), T) != G2 or abs(det(T)) != 1:
            bad.append((iso['first'], iso['second']))
    c('F', 'identifying isometries: T^T G T = G\' and det T = +-1', not bad and len(ident['isometries']) == 106, f"{len(ident['isometries'])} isometries; failures {bad[:3]}")

    # G. classification counts
    general = {r['label']: r['aut'] for r in rows if r['part'] == 'general'}
    classes = ident['classes']
    c('G', '12,546 identification classes = the general part', len(classes) == 12546 and {k['label']: k['aut'] for k in classes} == general)
    dist = load(CERT / 'general-part-distinctness.json')
    prof_ok = all(len({tuple(dist['unlabelled_profiles'][l]['key']) for l in b}) == len(b) for b in dist['signature_buckets'])
    c('G', 'Lemma 6.3: 12,519 classes with exact identifiers and 27 with distinct profiles',
      dist['classes_with_exact_identifier'] + dist['classes_without_exact_identifier'] == 12546
      and dist['classes_without_exact_identifier'] == len(dist['unlabelled_profiles']) == 27 and prof_ok
      and set(dist['unlabelled_profiles']) <= set(general) and '$"12,519"$' in body)
    c2 = load(CERT / 'c2xc2-part.json')
    c('G', 'the 67,284 lattices of S_(C_2^2), each with |Aut L| = 4',
      sorted(x['label'] for x in c2['lattices']) == sorted(r['label'] for r in rows if r['part'] == 'c2xc2') and all(x['aut'] == 4 for x in c2['lattices']))
    qt = load(CERT / 'quaternion-triples.json')
    quaternion = sum(v for k, v in cq.items() if k.startswith('Q'))
    c('G', 'quaternion triples: 849 lattices, 3 Gaussian classes each', len(qt['triples']) == quaternion == 849
      and len(qt['classes']) == 3 * 849 and all(len(t['keys']) == 3 for t in qt['triples']))
    sub = load(CERT / 'subtraction-records.json')
    cyclic = sum(v for k, v in cq.items() if k.startswith('C'))
    c('G', 'Lemma 6.6: 76,742 Gaussian classes minus 11,545 of S_gen = cyclic + 3 x quaternion',
      sub['list_classes'] - sub['distinct_keys'] == cyclic + 3 * quaternion == 65197
      and sub['keys_not_in_list'] == sub['hermitian_rows_without_key'] == sub['unitary_order_or_f_mismatches'] == 0
      and sub['trivial_involution_rows_without_identifier'] == sub['lemma_6_2_violations'] == 0)

    # G (continued). The row-level witnesses of Lemmas 6.3 and 6.6, replayed from general-part-witnesses.json
    # and the other certificates.  Inputs: the full groups and their complete conjugacy-class rows (types,
    # centralizer orders), the identifiers and Gaussian keys of the rows, the inverse and self-inverse
    # identifiers and the positive isometries.  Checked: the bookkeeping that turns them into the lemmas.
    wit = load(CERT / 'general-part-witnesses.json')
    inv = load(CERT / 'involution-gluing-masses.json')
    k824 = load(CERT / 'type-8-24-masses.json')
    gl = load(CERT / 'gaussian-list-weights.json')
    recs = wit['records']          # pair rows: [row, kind, type, identifier, centralizer order, self-conjugate]
    parent = {}

    def root(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def join(x, y):
        x, y = root(x), root(y)
        if x != y:
            parent[max(x, y)] = min(x, y)

    def classes_of_records():
        return len({root(r['item']) for r in recs})

    for r in recs:
        root(r['item'])
        for row in r['pair_rows']:
            if row[3] is not None:
                join(r['item'], row[3])
    steps = [classes_of_records()]
    for edges in (wit['inverse_edges'], wit['isometry_edges'][:39], wit['isometry_edges'][39:]):
        for a, b in edges:
            join(a, b)
        steps.append(classes_of_records())
    shown = re.search(r'the \$"17,736"\$ successful full-group records fall into \$"13,396"\$ classes by shared exact '
                      r'identifiers, into \$"12,585"\$ when the inverse identifiers are added, and into \$"12,546"\$ by the 39 '
                      r'isometries among the 66 unlabelled records; the 67 isometries', main_text)
    unlabelled = [r for r in recs if all(row[3] is None for row in r['pair_rows'])]
    c('G', '§6.3: the records join step by step into 13,396, 12,585, 12,546 classes (identifiers, inverse identifiers, '
      'the 39 isometries of the 66 unlabelled records); the 67 isometries of the large groups join none',
      len(recs) == 17736 and steps == [13396, 12585, 12546, 12546] and len(unlabelled) == 66
      and len(wit['isometry_edges']) == 106 and shown is not None, f'{len(recs)} records; classes after each step {steps}')
    label_at, aut_at, clash = {}, {}, []
    for s, label, aut in wit['database_sources']:
        t = root(s)
        if t in label_at:
            clash.append(label)
        label_at[t], aut_at[t] = label, aut
    roots = {root(x) for x in list(parent)}
    c('G', 'Lemma 6.3: the records, their identifiers, the inverse and the isometry edges join into 12,546 classes, '
      'each with exactly one lattice of the general part, of the same |Aut L| as its records',
      not clash and roots == set(label_at) == {root(r['item']) for r in recs} and len(roots) == 12546
      and sorted(label_at.values()) == sorted(general) and all(general[label_at[t]] == aut_at[t] for t in roots)
      and all(r['aut'] == aut_at[root(r['item'])] for r in recs),
      f"{len(recs)} records, {len(parent)} entries, {len(roots)} classes; two lattices in one class: {clash[:3]}")
    totals, missing = Counter(), Counter()
    for r in recs:
        for row in r['pair_rows']:
            totals[row[2]] += 1
            missing[row[2]] += row[3] is None
    exact = sorted(t for t in totals if not missing[t])
    witness, no_row = {}, []
    for i, r in enumerate(recs):
        lab = label_at[root(r['item'])]
        for row in r['pair_rows']:
            if row[3] is not None and not missing[row[2]]:
                witness.setdefault(lab, [i, row[0], row[2], row[3]])
        if not any(row[1] == 'odd' or row[4] > 4 or row[5] for row in r['pair_rows']):
            no_row.append(i)       # Lemma 6.2: a record of S_gen has a non-trivial row
    unanchored = set(general) - set(witness)
    prof = [tuple(dist['unlabelled_profiles'][l]['key']) for l in sorted(unanchored & set(dist['unlabelled_profiles']))]
    c('G', 'Lemma 6.3: the exact types from the rows (no row without identifier); a witness row of an exact type for '
      'each of 12,519 classes; the other 27 classes are those with saved profiles, pairwise different',
      exact == dist['exact_types'] and len(recs) == dist['successful_records']
      and {t: {'rows': totals[t], 'without_id': missing[t]} for t in totals if missing[t]} == dist['inexact_types']
      and witness == wit['component_witnesses'] and len(witness) == dist['classes_with_exact_identifier'] == 12519
      and unanchored == set(dist['unlabelled_profiles']) and len(prof) == len(set(prof)) == 27 and not no_row,
      f'{len(exact)} exact and {len(totals) - len(exact)} other types, {len(witness)} witnesses, {len(unanchored)} classes with profiles; records without a row of S_gen {no_row[:3]}')
    bad = [e for k in classes for e in k['list_entries'] if e not in parent or label_at[root(e)] != k['label']]
    for iso in ident['isometries']:
        for e in (iso['first'], iso['second']):
            if (e != iso['label']) if e.startswith('32.') else (e not in parent or label_at[root(e)] != iso['label']):
                bad.append(e)
    edges_ok = len(wit['isometry_edges']) == len(ident['isometries']) and all(
        all(e == p or (e.startswith('record:') and label_at[root(e)] == p) for e, p in zip(edge, (iso['first'], iso['second'])))
        for edge, iso in zip(wit['isometry_edges'], ident['isometries']))
    pub_inv, wit_inv = {tuple(x) for x in ident['inverse_identifiers']}, {tuple(x) for x in wit['inverse_edges']}
    high = wit['high_shell_key_rows']
    high_ok = [{'list_entry': e, 'lattice': label_at[root(recs[i]['item'])]} for i, _, e in high] == ident['high_shell_key_completions'] \
        and all(any(row[0] == j and row[3] == e for row in recs[i]['pair_rows']) for i, j, e in high)
    c('G', 'identification: every list entry and isometry endpoint of general-part-identification.json lies in the class '
      'of its lattice; its inverse identifiers and isometries are edges here; the high-shell key completions are rows of '
      'their lattice; every inverse edge lies in a class of the general part',
      not bad and edges_ok and high_ok and pub_inv <= wit_inv and len(wit_inv) == len(wit['inverse_edges']),
      f'{len(wit_inv)} inverse edges ({len(wit_inv - pub_inv)} not in general-part-identification); failures {bad[:3]}')
    trivial, bad = [], []
    for i, r in enumerate(recs):
        t = [row[3] for row in r['pair_rows'] if row[1] == 'inv' and row[4] == 4 and not row[5]]
        trivial += t
        if len(t) > 1 or (t and (r['aut'] % 8 != 4 or None in t or root(t[0]) != root(r['item']))):
            bad.append(i)
    lists = {p['list']: p for p in inv['pairs']}
    cand = {f"{p['list']}:{i}" for p in inv['pairs'] for i, (s, e) in enumerate(p['classes']) if s == 4 and e}
    cand |= {f'8-24:{i}' for i, x in enumerate(k824['lattices']) if int(x['aut']) == 2}
    selfinv = set(wit['self_inverse_entries'])
    remaining = cand - selfinv - set(parent)
    oriented = []
    for x in c2['lattices']:
        e, f = x['list_entry'], x['inverse_entry']
        name = e.rsplit(':', 1)[0]
        if name != '8-24' and name not in lists:
            bad.append(e)
            continue
        same_parts = name != '8-24' and lists[name]['F'] == lists[name]['K']
        if (f is not None) != same_parts or f == e:
            bad.append(e)
        oriented += [e] if f is None else [e, f]
    selfinv_ok = selfinv <= cand and all(lists[s.rsplit(':', 1)[0]]['F'] == lists[s.rsplit(':', 1)[0]]['K'] for s in selfinv)
    c('G', 'Lemma 6.6(i): the extremal involution entries with |Stab| = 4 and the entries of type (8, 24) with |Aut K| = 2, '
      'minus the self-inverse entries and the entries of S_gen, are exactly the oriented entries of the 67,284 lattices '
      'of S_(C_2^2) (inverse entry exactly for F = K); in S_gen at most one trivial involution row per record, '
      'only for |Aut L| = 4 mod 8',
      not bad and selfinv_ok and len(oriented) == len(set(oriented)) and set(oriented) == remaining
      and len(trivial) == sub['trivial_involution_rows'] and len(c2['lattices']) == 67284,
      f'{len(cand)} entries, {len(selfinv)} self-inverse ({len(selfinv & set(parent))} of them in S_gen), '
      f'{len(cand & set(parent))} in S_gen, {len(remaining)} remain; {len(trivial)} trivial rows ({len(set(trivial))} entries); failures {bad[:3]}')
    weight = {k: (u, f) for k, u, f in gl['classes']}
    owner, n_rows, bad = {}, 0, []
    for i, r in enumerate(recs):
        lab = label_at[root(r['item'])]
        for j, k, u, f in r['gaussian_rows']:
            n_rows += 1
            if weight.get(k) != (u, f) or r['aut'] % u or owner.setdefault(k, lab) != lab:
                bad.append((i, j))
    cq_rows = {r['label']: r for r in rows if r['part'] == 'cyclic-quaternion'}
    keys = []
    for label, aut, ks in wit['cq_lattices']:
        r, w = cq_rows.get(label), [weight.get(k) for k in ks]
        if r is None or r['aut'] != aut or None in w or any(k in owner for k in ks):
            bad.append(label)
            continue
        s = r['group']['structure']
        if len(ks) == 1:
            ok = s == f'C{aut}' and w[0] == (aut, 2)
        else:
            ok = s == f'Q{aut}' and len(ks) == 3 and all(f == 1 for _, f in w) and sorted(u for u, _ in w) == [4, 4, aut // 2]
        if not ok:
            bad.append(label)
        keys += ks
    n_q = sum(len(ks) == 3 for _, _, ks in wit['cq_lattices'])
    c('G', 'Lemma 6.6(ii), (iii): the Gaussian keys of Proposition 4.6 minus the keys of S_gen are exactly the keys of the '
      '63,499 lattices of S_CQ: one key with f = 2 and |U| = |Aut L| for C_(2^k), three keys with f = 1 and '
      '|U| = 4, 4, |Aut L|/2 for Q_(2^k); in S_gen every key has the |U| and f of the list and lies in one class',
      not bad and len(keys) == len(set(keys)) and set(keys) == set(weight) - set(owner)
      and sorted(x[0] for x in wit['cq_lattices']) == sorted(cq_rows) and n_rows == sub['hermitian_rows']
      and len(owner) == sub['distinct_keys'] and n_q == len(qt['triples']) and 3 * n_q == len(qt['classes']),
      f'{n_rows} Hermitian rows of S_gen with {len(owner)} keys; {len(keys)} keys in {len(wit["cq_lattices"])} lattices ({n_q} triples); failures {bad[:3]}')

    # G (continued). The quaternion triples of S_CQ (Proposition 6.5(ii), Lemma 6.6(iii)) in integers.
    # Inputs: that the keys and |U| belong to (G, W) and that different triples are different lattices.
    gram = {x['key']: x['gram'] for x in qt['classes']}
    bad = []
    for x in qt['classes']:
        G = x['gram']
        shape = len(G) == 32 and all(len(r) == 32 for r in G)
        if not (shape and all(G[i][j] == G[j][i] for i in range(32) for j in range(i)) and all(G[i][i] % 2 == 0 for i in range(32))):
            bad.append(x['key'])
            continue
        minors = leading_minors(G)
        if not (len(minors) == 32 and all(v > 0 for v in minors) and minors[-1] == 1):
            bad.append(x['key'])
    c('G', 'quaternion triples: every Gram matrix is symmetric, even, positive definite (leading minors > 0) and has det 1',
      not bad and len(gram) == len(qt['classes']) == 2547, f'{len(gram)} Gram matrices; failures {bad[:3]}')
    minus_one = [[-int(i == j) for j in range(32)] for i in range(32)]
    bad = [x['key'] for x in qt['classes']
           if matmul(matmul(x['W'], x['gram']), transpose(x['W'])) != x['gram'] or matmul(x['W'], x['W']) != minus_one]
    c('G', 'quaternion triples: W G W^T = G and W^2 = -1 (W is multiplication by i on row vectors)', not bad, f'failures {bad[:3]}')
    bad, owner_q = [], {}
    for t in qt['triples']:
        ks, par = t['keys'], {k: k for k in t['keys']}

        def find(x):
            while par[x] != x:
                x = par[x]
            return x
        for iso in t['isometries']:
            a, b, T = iso['first'], iso['second'], iso['T']
            if a not in par or b not in par or matmul(matmul(transpose(T), gram[a]), T) != gram[b] or abs(det(T)) != 1:
                bad.append((t['label'], a, b))
            else:
                par[find(a)] = find(b)
        if len(set(ks)) != 3 or len({find(k) for k in ks}) != 1:
            bad.append(t['label'])
        for k in ks:
            owner_q.setdefault(k, t['label'])
    n_iso = sum(len(t['isometries']) for t in qt['triples'])
    c('G', 'quaternion triples: in each triple, integer isometries T^T G_first T = G_second with det T = +-1 join the three classes',
      not bad and n_iso == 1698 and len(owner_q) == 3 * len(qt['triples']) == len(gram) and set(owner_q) == set(gram),
      f'{len(qt["triples"])} triples, {n_iso} isometries; failures {bad[:3]}')
    cq_keys = {label: (aut, sorted(ks)) for label, aut, ks in wit['cq_lattices'] if len(ks) == 3}
    c('G', 'quaternion triples: the keys of each triple are the three keys of its label in cq_lattices; |U| and f as in gaussian-list-weights',
      {t['label']: (t['aut'], sorted(t['keys'])) for t in qt['triples']} == cq_keys and len(cq_keys) == len(qt['triples'])
      and all(weight.get(x['key']) == (x['U'], x['f']) for x in qt['classes']),
      f'{len(cq_keys)} labels with three keys')

    # H. 9-torsion
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, '-I', str(ROOT / 'verify/check_exponent9.py'), str(Path(tmp) / 'e9.json')], capture_output=True, text=True)
        rerun = load(Path(tmp) / 'e9.json') if r.returncode == 0 else {}
    saved = load(CERT / 'exponent-9-reduction.json')
    c('H', 'the 9-torsion finite check reruns to the certificate',
      r.returncode == 0 and rerun.get('cases') == saved['cases'] and rerun.get('admissible_graphs') == saved['admissible_graphs'], r.stderr[-200:])
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, '-I', str(ROOT / 'verify/check_quadratic2_orbits.py'), str(Path(tmp) / 'q2.json')], capture_output=True, text=True)
        rerun = load(Path(tmp) / 'q2.json') if r.returncode == 0 else {}
    saved = load(CERT / 'quadratic-2-form-orbits.json')
    c('H', 'the orbit check of the elementary quadratic 2-forms (Appendix D.4.6) reruns to the certificate',
      r.returncode == 0 and rerun.get('cases') == saved['cases'] and rerun.get('forms_checked') == saved['forms_checked'] == 230, r.stderr[-200:])

    # I. completeness identities
    comps = inv['components']
    ok_comp, bad = True, []
    for k, comp in enumerate(comps):
        G = comp['gram']
        n = len(G)
        sym_even = all(G[i][j] == G[j][i] for i in range(n) for j in range(n)) and all(G[i][i] % 2 == 0 for i in range(n))
        inverse = inverse_fraction(G)
        two_el = all((2 * x).denominator == 1 for r in inverse for x in r)
        uses = [p for p in inv['pairs'] if k in (p['F'], p['K'])]
        s_ok = all(abs(det(G)) == 2 ** p['s'] for p in uses)
        delta = 0 if all(inverse[i][i].denominator == 1 for i in range(n)) else 1
        d_ok = all(p['delta'] == delta for p in uses)
        if not (sym_even and two_el and s_ok and d_ok):
            bad.append(k)
    c('I', 'involution components: even, 2-elementary, det 2^s, and delta from the Gram matrix', not bad,
      f'{len(comps)} components; failures {bad[:5]}')
    # Lemma 2.1: |O(q_F)| of every pair from the closed formulas (rank d of F from its Gram matrix, s and delta checked
    # above; for odd p the type p-(z,f)-s), so that the right side of (4.1) does not rest on the program for R(q).
    od = load(CERT / 'odd-prime-masses.json')
    bad = [p['list'] for p in inv['pairs'] if order_O_two(len(comps[p['F']]['gram']), int(p['s']), int(p['delta'])) != int(p['order_O_q'])]
    bad += [t['type'] for t in od['types'] if order_O_odd(t['p'], t['f'], t['s']) != int(t['order_O_q'])]
    c('I', 'Lemma 2.1: |O(q_F)| of every pair of components equals its closed formula',
      not bad, f"{len(inv['pairs'])} involution lists, {len(od['types'])} types of odd prime order; failures {bad[:5]}")
    bad = []
    for p in inv['pairs']:
        lhs = sum((Q(1, s) for s, _ in p['classes']), Q(0))
        kF, kK = p['kernels']
        aF, aK = int(p['aut_F']), int(p['aut_K'])
        if (lhs != Q(int(p['order_O_q']), aF * aK) or aF != int(comps[p['F']]['aut']) or aK != int(comps[p['K']]['aut'])
                or any(s % (kF * kK) for s, _ in p['classes'])):
            bad.append(p['list'])
    c('I', 'identity (4.1) for the 83 pairs of involution components', not bad and len(inv['pairs']) == 83,
      f"{sum(len(p['classes']) for p in inv['pairs'])} classes; failures {bad[:5]}")
    by_type = {}
    for p in inv['pairs']:
        t = tuple(p['type'])
        cl, ex, m = by_type.get(t, (0, 0, Q(0)))
        m_ext = sum((Q(1, s) for s, e in p['classes'] if e), Q(0))
        # both orientations: a fixed type is listed once; for (16, 16) a pair with F = K holds both, F != K counts twice
        weight = 2 if t != (16, 16) or p['F'] != p['K'] else 1
        by_type[t] = (cl + len(p['classes']), ex + sum(e for _, e in p['classes']), m + weight * m_ext)
    rows5 = all(f'[${fm(cl)}$], [${fm(ex)}$], [${fmf(m)}$]' in main_text for cl, ex, m in by_type.values())
    c('I', 'Table 5: classes, extremal classes and M of the four glued involution types', rows5 and len(by_type) == 4,
      {f'{t[0]}-{t[1]}': (cl, ex) for t, (cl, ex, _) in by_type.items()})
    rows6 = []
    for p in inv['pairs']:
        if p['type'] != [16, 16]:
            continue
        ext = [s for s, e in p['classes'] if e]
        m_ext = sum((Q(1, s) for s in ext), Q(0))
        m_txt = f'{fm(m_ext.numerator)} slash {fm(m_ext.denominator)}' if p['F'] == p['K'] else f'2 dot {fm(m_ext.numerator)} slash {fm(m_ext.denominator)}'
        n_cl, n_ext, n_4 = len(p['classes']), len(ext), sum(s == 4 for s in ext)
        rows6.append(f'[${fm(n_cl)}$], [${fm(n_ext)}$], [{ft(n_4)}], [${m_txt}$]' in main_text
                     or f'[{ft(n_cl)}], [{ft(n_ext)}], [{ft(n_4)}], [${m_txt}$]' in main_text
                     or f'[${fm(n_cl)}$], [${fm(n_ext)}$], [${fm(n_4)}$], [${m_txt}$]' in main_text)
    c('I', 'Table 6: the eleven pairs of type (16, 16)', len(rows6) == 11 and all(rows6), rows6.count(False))

    # list entries of the classified lattices point to extremal classes of the right order
    cls = {p['list']: p['classes'] for p in inv['pairs']}

    def entry(e):
        name, i = e.rsplit(':', 1)
        i = int(i)
        if name == '8-24':
            return (2 * int(k824['lattices'][i]['aut']), True) if i < len(k824['lattices']) else None
        return tuple(cls[name][i]) if name in cls and i < len(cls[name]) else None

    bad = [x['list_entry'] for x in c2['lattices'] if entry(x['list_entry']) != (4, True)]
    c('I', 'every lattice of S_(C_2^2) is an extremal involution class with |Stab| = 4', not bad, f'failures {bad[:3]}')
    inv_entries = [e for k in classes for e in k['list_entries'] if e.rsplit(':', 1)[0] in cls or e.startswith('8-24:')]
    bad = [e for e in inv_entries if entry(e) is None or not entry(e)[1]]
    c('I', 'every involution entry of S_gen is an extremal class of its list', not bad, f'{len(inv_entries)} entries; failures {bad[:3]}')

    # type (8, 24): identity (B.1)
    lat = k824['lattices']
    census = [dict(map(tuple, x['roots_census'])) for x in lat]
    rhs = {r: sum((Q(cs.get(r, 0), int(x['aut'])) for cs, x in zip(census, lat)), Q(0)) for r in (0, 48, 72)}
    lhs = {0: sum((Q(1, s) for s in k824['orbit_stabilizers']['Leech']), Q(0)),
           48: sum((Q(1, s) for s in k824['orbit_stabilizers']['A1^24']), Q(0)),
           72: sum((Q(g, s) for s, g in k824['orbit_terms_A2_12']), Q(0)) / 15}
    shown = re.search(r'g\(W_3\) slash \|"Stab"\(W_3\)\| = ' + NUM + ' slash ' + NUM + r'\$', body)
    c('I', 'identity (B.1) for the Leech lattice and A_1^24 (all orbits of good W) and A_2^12 (the 106,624 orbits of good W_3, divided by 15, as in Appendix B.3)',
      all(lhs[r] == rhs[r] for r in lhs) and len(k824['orbit_terms_A2_12']) == 106624 and lhs[72] == Q(k824['left_side_A2_12']) and all(sum(cs.values()) == 270 for cs in census)
      and shown is not None and Q(num(shown[1]), num(shown[2])) == lhs[72], {r: str(v) for r, v in lhs.items()})
    # The 14,482 lattices K are distinct: 14,462 by a complete invariant (the canonical form of the graph of the 720
    # vectors, which generate K^#), 16 by that form completed by the orbit of N/N', and the 4 others because each has an
    # order |Aut K| that no other K has (§B.3).
    froms = Counter(x.get('aut_from') for x in lat)
    aut_count = Counter(int(x['aut']) for x in lat)
    c('I', 'type (8, 24): the four lattices K separated by |Aut K| have an order that no other K has',
      froms == {'graph': 14462, 'orbit': 16, 'qfauto': 4} and all(aut_count[int(x['aut'])] == 1 for x in lat if x['aut_from'] == 'qfauto'),
      dict(froms))
    split = Counter(min(cs) for cs in census)
    found = sum((Q(cs.get(72, 0), int(x['aut'])) for cs, x in zip(census, lat) if min(cs) == 72), Q(0))
    m824 = sum((Q(1, int(x['aut'])) for x in lat), Q(0))
    c('I', 'Table B3 and Table 5: 1216 + 13,231 + 35 = 14,482 lattices K; the 35 fill 429583/768; M',
      (split[0], split[48], split[72], len(lat)) == (1216, 13231, 35, 14482) and found == Q(429583, 768)
      and '$"429,583" slash 768$' in body and '[$"1,244"$], [$"1,216"$]' in body and '[$"376,306"$], [$"13,231"$]' in body
      and len(k824['orbit_stabilizers']['Leech']) == 1244 and len(k824['orbit_stabilizers']['A1^24']) == 376306
      and f'[$"14,482"$], [$"14,482"$], [${fmf(m824)}$]' in main_text, dict(split))

    # types of odd prime order: identity (4.1), Table B4, identity (B.2)
    od = load(CERT / 'odd-prime-masses.json')
    gF, gZ = od['genera_F'], od['genera_Z']
    bad = []
    for g in gF:
        if g['target'] is not None and sum((Q(1, int(x['aut'])) for x in g['classes']), Q(0)) != Q(g['target']):
            bad.append(('F', g['p'], g['f'], g['s']))
        n4, n = sum(x['min4'] for x in g['classes']), len(g['classes'])
        cell = f"$({g['f']}, {g['s']})$: {ft(n4)}" + (f' [{ft(n)}]' if g['completeness'] == 'mass' else '')
        if g['completeness'] != 's = f' and g['p'] != 31 and cell not in body:
            bad.append(('F cell', cell))
    for g in gZ:
        if g['classes'] is None:        # rank 1, p = 31: every class has |U| = 62
            ok31 = Q(g['mass_genus']) == g['classes_total'] * Q(1, g['unitary_order_every_class']) and Q(g['target']) == g['rootless_classes'] * Q(1, 62)
            if not ok31:
                bad.append(('Z', 31))
            continue
        if sum((Q(1, int(x['U'])) for x in g['classes']), Q(0)) != Q(g['target']):
            bad.append(('Z', g['p'], g['z'], g['s']))
        if g['table_B4'] and g['p'] not in (17, 31) and f"$({g['z']}, {g['s']})$: {ft(len(g['classes']))}" not in body:
            bad.append(('Z cell', g['p'], g['z'], g['s']))
    n_b4 = sum(g['table_B4'] for g in gZ)
    c('I', 'Table B4: the component lists reach the mass of their genus or the rootless mass (for the lists completed by (B.2), a check against the rootless mass)', not bad and len(gF) == 16 and n_b4 == 19,
      f'{len(gF)} genera of F, {n_b4} of Z in Table B4 and {len(gZ) - n_b4} further levels of the (B.2) chains; failures {bad[:4]}')
    bad, ext_by_p, n_pairs = [], Counter(), 0
    for t in od['types']:
        O = int(t['order_O_q'])
        for q in t['pairs']:
            aF = 1 if q['aut_F'] is None else int(q['aut_F'])
            aZ = int(q['aut_Z'])
            if q['F'] is not None and gF[t['genus_F']]['classes'][q['F']]['aut'] != q['aut_F']:
                bad.append((t['type'], 'aut_F'))
            if q['Z'] is not None and gZ[t['genus_Z']]['classes'][q['Z']]['U'] != q['aut_Z']:
                bad.append((t['type'], 'aut_Z'))
            stab = q['stabilizers']
            if 'decompositions_all' in q:
                O1 = t['order_O_q1']
                ok = (q['decompositions_all'] * O1 * O1 == O and Q(q['mass_all']) == Q(O, aF * aZ)
                      and sum((Q(1, x) for x in stab), Q(0)) == Q(q['mass_good']) == Q(q['decompositions_good'] * O1 * O1, aF * aZ))
            elif stab is None:
                ok = Q(q['mass']) == Q(O, aF * aZ)
            else:
                ok = sum((Q(1, x) for x in stab), Q(0)) == Q(O, aF * aZ)
                n_pairs += 1
            if not ok:
                bad.append((t['type'], q['F'], q['Z']))
            if q['extremal']:
                ext_by_p[t['p'] == 3] += sum(q['extremal'])
    c('I', 'identity (4.1) for the pairs of components of the 19 types of odd prime order', not bad and len(od['types']) == 19,
      f'{n_pairs} pairs with all double cosets; failures {bad[:4]}')
    c('I', 'Section 4.2: 6617 extremal pairs of order 3 and 218 of orders 5, 7, 17, 31',
      (ext_by_p[True], ext_by_p[False]) == (6617, 218) and 'have $"6,617"$ extremal pairs' in body and 'have $218$ extremal pairs' in body,
      dict(ext_by_p))
    # identifiers of order 3 in the full-group records of S_gen: <z>-<f>-<s>:<i>-<j>/<k> is the k-th double coset of
    # the i-th F of minimum >= 4 and the j-th Z; it must be extremal with |Stab(phi)| = |C(h)| of the row
    by_type = {(t['z'], t['f'], t['s']): t for t in od['types'] if t['p'] == 3}
    bad, n_odd = [], 0
    for rec in wit['records']:
        for row_index, kind, rtype, ident, cent, _ in rec['pair_rows']:
            if kind != 'odd' or ident is None:
                continue
            n_odd += 1
            try:
                head, rest = ident.rsplit(':', 1)
                z, f, s = map(int, head.split('-'))
                fz, k = rest.split('/')
                i, j = map(int, fz.split('-'))
                t = by_type[z, f, s]
                min4 = [n for n, x in enumerate(gF[t['genus_F']]['classes']) if x['min4']]
                (q,) = [q for q in t['pairs'] if q['F'] == min4[i] and q['Z'] == j]
                ok = rtype == f'3-{z}-{f}-{s}' and q['stabilizers'][int(k)] == cent and q['extremal'][int(k)]
            except (KeyError, ValueError, IndexError, TypeError):
                ok = False
            if not ok:
                bad.append(ident)
    c('I', 'every identifier of order 3 in the full-group records of S_gen is an extremal double coset of its list with |Stab| = |C(h)|',
      not bad and n_odd > 0, f'{n_odd} rows; failures {bad[:3]}')
    bad = []
    for st in od['b2_steps']:
        lhs = sum((Q(st['N'], u) for u in st['lower']), Q(0))
        if st['totals_only']:
            ok = lhs == Q(st['identity_left']) == Q(st['identity_right'])
        else:
            ok = lhs == sum((Q(g, u) for u, g in st['upper']), Q(0))
        ok &= st['N'] == (st['p'] ** (st['z'] - st['s']) - 1) // (st['p'] - 1)
        if not ok:
            bad.append((st['p'], st['z'], st['s']))
    c('I', 'identity (B.2) at the saved steps, term by term', not bad and not any(st['totals_only'] for st in od['b2_steps']),
      f"{len(od['b2_steps'])} steps; failures {bad}")
    # The chains of (B.2): each step goes from the list of (z, s + 2) to that of (z, s), and every list completed by (B.2)
    # is reached by steps from a root (z, z) whose list reaches the rootless mass.
    by_genus = {(g['p'], g['z'], g['s']): g for g in gZ}
    step_to = {(st['p'], st['z'], st['s']): st for st in od['b2_steps']}
    bad, roots = [], set()
    for (p, z, s), st in step_to.items():
        up, lo = by_genus.get((p, z, s + 2)), by_genus.get((p, z, s))
        if up is None or lo is None or [u for u, _ in st['upper']] != [int(x['U']) for x in up['classes']] \
                or st['lower'] != [int(x['U']) for x in lo['classes']]:
            bad.append(('levels', p, z, s))
    for (p, z, s), g in by_genus.items():
        if '(B.2)' not in g['completeness']:
            continue
        t = s
        while (p, z, t) in step_to and t < z:
            t += 2
        root = by_genus.get((p, z, t))
        if t != z or root is None or 'rootless' not in root['completeness'] \
                or sum((Q(1, int(x['U'])) for x in root['classes']), Q(0)) != Q(root['target']):
            bad.append(('chain', p, z, s))
        else:
            roots.add((p, z))
    c('I', 'the chains of (B.2): each step joins the lists of consecutive levels, and each chain starts at a list (z, z) of rootless mass',
      not bad and len(step_to) == len(od['b2_steps']), f'roots {sorted(roots)}; failures {bad}')

    # Lemma 6.3(ii) (exhaustion of S_gen): for every list of §4.1-§4.2, the rows of the 12,546 lattices of S_gen, one
    # record per lattice, weighted 1/|C(h)| (2/|C(sigma)| for an involution entry of a list with F = K whose sigma is not
    # conjugate to -sigma), add up to the extremal mass of the list less the classes T_0 assigned to S_(C_2^2); T_0
    # consists of trivial classes (|Stab| = 4, not self-inverse) none of which, nor its partner, is the identifier of a
    # row.  With injectivity (distinct rows are distinct classes), every other extremal class is an entry of S_gen.
    recs_by_item = {r['item']: r for r in wit['records']}
    sources = [s_ for s_, _, _ in wit['database_sources']]
    selfinv = set(wit['self_inverse_entries'])
    t0 = set()
    for x in c2['lattices']:
        t0.add(x['list_entry'])
        if x['inverse_entry']:
            t0.add(x['inverse_entry'])
    lists = {q['list']: q for q in inv['pairs']}
    row_mass, row_ids, unknown = {}, set(), []
    for s_ in sources:
        for row, kind, typ, ident_, cen, selfconj in recs_by_item[s_]['pair_rows']:
            name = typ[4:] if typ.startswith('inv:') else None
            twice = name in lists and lists[name]['F'] == lists[name]['K'] and not selfconj
            row_mass[typ] = row_mass.get(typ, Q(0)) + Q(2 if twice else 1, cen)
            if ident_ is not None:
                row_ids.add(ident_)
    expected = {}
    for name, q in lists.items():
        expected['inv:' + name] = sum((Q(1, st) for i, (st, e) in enumerate(q['classes']) if e and f'{name}:{i}' not in t0), Q(0))
    expected['inv:8-24'] = sum((Q(1, 2 * int(x['aut'])) for i, x in enumerate(k824['lattices']) if f'8-24:{i}' not in t0), Q(0))
    for t in od['types']:
        key = f"3-{t['z']}-{t['f']}-{t['s']}" if t['p'] == 3 else f"odd:p={t['p']}:f={t['f']}"
        m = Q(0)
        for q in t['pairs']:
            if 'mass_good' in q:
                m += Q(q['mass_good'])
            elif q.get('stabilizers') is not None:
                m += sum((Q(1, int(st)) for st, e in zip(q['stabilizers'], q['extremal']) if e), Q(0))
        expected[key] = expected.get(key, Q(0)) + m

    def entry_class(e):          # [|Stab(phi)|, extremal] of a list entry '<list>:<index>'
        name, i = e.rsplit(':', 1)
        if name == '8-24':
            return [2 * int(k824['lattices'][int(i)]['aut']), True]
        return lists[name]['classes'][int(i)] if name in lists else None
    t0_ok = all(entry_class(e) is not None and int(entry_class(e)[0]) == 4 and entry_class(e)[1] and e not in selfinv for e in t0) \
        and not (t0 & row_ids)
    # Pairs of odd prime order saved by their mass only add nothing to the expected side; they must be pairs whose
    # gluings are not extremal (3-(6,20)-6 with Z = K_12, §4.2), as the certificate states.
    mass_only = [(t['type'], q['F']) for t in od['types'] for q in t['pairs'] if q.get('stabilizers') is None and 'mass_good' not in q]
    mass_only_ok = all(t['type'] == '3-(6,20)-6' and 'not extremal' in t['cosets_saved'] for t in od['types']
                       for q in t['pairs'] if q.get('stabilizers') is None and 'mass_good' not in q)
    # Injectivity of rows -> classes where the rows carry identifiers: an involution row stands for {x, x'} (x' the
    # inverse class for F = K), and these sets are distinct over all rows; for F = K a row that is not self-conjugate
    # has a partner and a self-conjugate one has none; the centralizer order of the row is |Stab| of its class; odd rows
    # with identifiers have distinct identifiers.  Rows without identifiers rest on the records (Table 10).
    partner = {}
    for a_, b_ in wit['inverse_edges']:
        partner.setdefault(a_, set()).add(b_)
        partner.setdefault(b_, set()).add(a_)
    seen_inv, seen_odd, inj_bad = set(), set(), []
    for s_ in sources:
        for row, kind, typ, ident_, cen, selfconj in recs_by_item[s_]['pair_rows']:
            if ident_ is None:
                continue
            if kind == 'inv':
                name = typ[4:]
                fk = name in lists and lists[name]['F'] == lists[name]['K']
                part = partner.get(ident_, set())
                if fk and (len(part) != (0 if selfconj else 1)):
                    inj_bad.append(('partner', ident_))
                key = frozenset({ident_} | (part if fk else set()))
                cls = entry_class(ident_)
                if key in seen_inv or cls is None or int(cls[0]) != cen:
                    inj_bad.append(('inv', ident_))
                seen_inv.add(key)
            else:
                if ident_ in seen_odd:
                    inj_bad.append(('odd', ident_))
                seen_odd.add(ident_)
    c('I', 'Lemma 6.3(ii): rows with identifiers give distinct classes, with |C| = |Stab|; pairs saved by mass only are not extremal',
      not inj_bad and mass_only_ok, f'{len(seen_inv)} involution entries, {len(seen_odd)} odd rows with identifiers; '
      f'failures {inj_bad[:3]}; mass-only pairs {mass_only}')
    bad = sorted(k for k in set(expected) | set(row_mass) if expected.get(k) != row_mass.get(k, Q(0)))
    c('I', 'Lemma 6.3(ii): for every list, the rows of the 12,546 lattices of S_gen add up to its extremal mass less the trivial classes T_0 of S_(C_2^2)',
      len(sources) == 12546 and all(s_ in recs_by_item for s_ in sources) and not bad and t0_ok and len(expected) == 83 + 1 + 18,
      f'{len(expected)} lists; |T_0| = {len(t0)}; differing {bad[:5]}; T_0 trivial and unmet: {t0_ok}')

    # Proposition 4.6: sum f/|U| = M_2
    total = sum((Q(f, u) for _, u, f in gl['classes']), Q(0))
    c('I', 'Proposition 4.6: sum of f/|U| over the 76,742 Gaussian classes equals M_2',
      total == Q(gl['rootless_mass']) == value and len(gl['classes']) == len({k for k, _, _ in gl['classes']}) == 76742 == sub['list_classes']
      and all(f in (1, 2) for _, _, f in gl['classes']) and '$"76,742"$' in main_text, total)
    # The methods of the weights (proof of Proposition 4.6): the counts as stated, and every canonical-form class is a
    # discrete class (key k...) with |U| = 4 and f = 2, as the argument of the proof gives.
    wm = gl.get('weight_methods', [])
    by = Counter(wm)

    shown = re.search(r'For \$' + NUM + r'\$ lattices two methods were applied \(bliss and the search for \$' + NUM + r'\$, '
                      r'bliss and the canonical form for \$' + NUM + r'\$\), and they agree\. Of the others, \$' + NUM + r'\$ have '
                      r'their weight from the canonical form alone', main_text)
    shown2 = re.search(r'\$' + NUM + r'\$ have their weight from the backtrack search alone and \$' + NUM + r'\$ from bliss alone', main_text)
    stated = (tuple(map(num, shown.groups())) + tuple(map(num, shown2.groups()))) if shown and shown2 else None
    counted = (by['bliss+search'] + by['bliss+canon'], by['bliss+search'], by['bliss+canon'], by['canon'], by['search'], by['bliss'])
    canon_ok = all(k.startswith('k') == ('canon' in m.split('+')) and (('canon' not in m.split('+')) or (u, f) == (4, 2))
                   for (k, u, f), m in zip(gl['classes'], wm))
    c('I', 'Proposition 4.6: the methods of the weights, counted as in its proof; canonical-form classes are the k-keys, with |U| = 4, f = 2',
      len(wm) == len(gl['classes']) and all(wm) and set(by) <= {'canon', 'search', 'bliss', 'bliss+canon', 'bliss+search'}
      and stated == counted and canon_ok, f'stated {stated}; counted {counted}')

    failures = [x for x in c.results if not x['ok']]
    report = {'checks': c.results, 'failures': len(failures), 'verifier_sha256': sha256(Path(__file__))}
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(json.dumps(report, indent=1) + '\n')
    print(f'{len(c.results) - len(failures)} of {len(c.results)} checks passed')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
