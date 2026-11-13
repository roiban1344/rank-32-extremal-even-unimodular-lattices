"""Finite check of the element orbits of elementary quadratic 2-forms (Appendix D.4.6).

Usage: python3 verify/check_quadratic2_orbits.py OUT.json

For q = <1/2>^a + <3/2>^b on F_2^(a+b), 1 <= a + b <= 20, the counting of Appendix D.4.6
uses that two vectors are in one orbit of O(q) when they are both 0, both the characteristic
vector, or neither and q takes the same value on them.  A vector is described up to coordinate
permutations (which preserve q) by its support counts (i, j) on the two kinds of coordinates.
The program joins support counts by two kinds of actual isometries:
  * the transvection t_w(x) = x + (x.w) w for w of even weight with 2 q(w) = 2 mod 4, which
    preserves q; it maps support counts (i, j) to (i + r - 2u, j + s - 2v), where (r, s) are the
    support counts of w and (u, v) those of the overlap, for u + v odd;
  * for a, b >= 2, the map (I + J) P on two coordinates of each kind, J the all-ones matrix and
    P the exchange of the two pairs; its columns keep their q-values and are orthonormal for the
    dot pairing.
It then checks that the connected components are exactly the classes named above.  The
certificate certificates/quadratic-2-form-orbits.json records the result.  This proves the
finite orbit statement; it is not the classification of 2-elementary forms nor the mass formula.
"""
import json
import sys
from collections import defaultdict
from math import comb
from pathlib import Path


def check_form(a, b):
    patterns = [(i, j) for i in range(a + 1) for j in range(b + 1)]
    parent = {x: x for x in patterns}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def join(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[max(x, y)] = min(x, y)

    directions = [(r, s) for r in range(a + 1) for s in range(b + 1) if (r + 3 * s) % 4 == 2]
    for i, j in patterns:
        for r, s in directions:
            assert (r + s) % 2 == 0
            for u in range(max(0, r - a + i), min(i, r) + 1):
                for v in range(max(0, s - b + j), min(j, s) + 1):
                    if (u + v) % 2 == 0:
                        continue
                    target = (i + r - 2 * u, j + s - 2 * v)
                    assert target in parent
                    assert (target[0] + 3 * target[1] - i - 3 * j) % 4 == 0
                    join((i, j), target)
        if a >= 2 and b >= 2:
            for u in range(max(0, 2 - a + i), min(i, 2) + 1):
                for v in range(max(0, 2 - b + j), min(j, 2) + 1):
                    if (u + v) % 2:
                        target = (i - u + 2 - v, j - v + 2 - u)
                    else:
                        target = (i - u + v, j - v + u)
                    assert target in parent
                    assert (target[0] + 3 * target[1] - i - 3 * j) % 4 == 0
                    join((i, j), target)

    def label(x):
        if x == (0, 0):
            return 'zero'
        if x == (a, b):
            return 'characteristic'
        return 'norm:' + str((x[0] + 3 * x[1]) % 4)

    components = defaultdict(list)
    for x in patterns:
        components[find(x)].append(x)
    by_label = defaultdict(list)
    output = []
    for component in components.values():
        labels = {label(x) for x in component}
        assert len(labels) == 1, (a, b, component, labels)
        tag = next(iter(labels))
        by_label[tag].append(component)
        output.append({'label': tag, 'support_patterns': len(component),
                       'vectors': sum(comb(a, i) * comb(b, j) for i, j in component)})
    assert all(len(cs) == 1 for cs in by_label.values()), (a, b, by_label)
    assert sum(c['vectors'] for c in output) == 2 ** (a + b)
    return {'a': a, 'b': b, 'rho': a + b, 'sigma_mod8': (a - b) % 8,
            'directions_by_support': len(directions),
            'orbits': sorted(output, key=lambda x: x['label'])}


def main():
    cases = [check_form(a, rho - a) for rho in range(1, 21) for a in range(rho + 1)]
    Path(sys.argv[1]).write_text(json.dumps({'forms_checked': len(cases), 'max_rank': 20, 'cases': cases}) + '\n')


if __name__ == '__main__':
    main()
