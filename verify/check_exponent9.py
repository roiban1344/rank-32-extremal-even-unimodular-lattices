"""Finite check of the 9-torsion reduction (Appendix D.4.6) by integer modular arithmetic.

Usage: python3 verify/check_exponent9.py OUT.json

A = (Z/9)^2 with q(x) = 8 (x_1^2 + x_2^2) / 9 and B with the negative form.  For every
map m in GL_2(F_3) on 3A, the graph of m and its orthogonal complement are computed;
the quotient is elementary exactly for the 8 maps that are anti-isometries after
scaling.  The certificate certificates/exponent-9-reduction.json records the result.
"""
import itertools
import json
import sys
from pathlib import Path


def mv(m, x, modulus):
    return ((m[0] * x[0] + m[1] * x[1]) % modulus,
            (m[2] * x[0] + m[3] * x[1]) % modulus)


def main():
    rows = []
    vectors = list(itertools.product(range(9), repeat=4))
    for m in itertools.product(range(3), repeat=4):
        if (m[0] * m[3] - m[1] * m[2]) % 3 == 0:
            continue
        graph = {(3 * x, 3 * y, *(3 * c for c in mv(m, (x, y), 3)))
                 for x, y in itertools.product(range(3), repeat=2)}
        # Pairing against each of the two graph generators gives orthogonality.
        generators = [(3, 0, 3 * m[0], 3 * m[2]), (0, 3, 3 * m[1], 3 * m[3])]
        perp = [v for v in vectors if all(
            (8 * (v[0] * g[0] + v[1] * g[1] - v[2] * g[2] - v[3] * g[3])) % 9 == 0
            for g in generators)]
        elementary = all(tuple(3 * c % 9 for c in v) in graph for v in perp)
        orthogonal = all((sum(c * c for c in mv(m, x, 3)) - sum(c * c for c in x)) % 3 == 0
                         for x in itertools.product(range(3), repeat=2))
        assert elementary == orthogonal
        assert len(graph) == 9 and len(perp) == 729
        rows.append({'matrix': m, 'scaled_anti_isometry': orthogonal,
                     'quotient_elementary': elementary, 'quotient_size': len(perp) // len(graph)})
    assert len(rows) == 48 and sum(r['scaled_anti_isometry'] for r in rows) == 8
    report = {'GL2_3_maps': 48, 'admissible_graphs': 8, 'cases': rows}
    Path(sys.argv[1]).write_text(json.dumps(report, indent=2) + '\n')
    print('All 48 maps checked; exactly 8 give an elementary quotient.')


if __name__ == '__main__':
    main()
