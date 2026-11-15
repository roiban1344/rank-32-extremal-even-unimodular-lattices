"""The normal form of rank32.gluing with less memory and the same keys.

rank32.gluing is not changed; this module holds a subclass that does the same computation and
splits only the step that takes a lot of memory.

The memory went to the chain on the F side in Glue.__init__ (Glue._fixed_chain). At each stage of
the chain, Node.orbits_of is called on all 2^s points. A node of order <= SMALL keeps all its
elements, and orbits_of then builds an int64 array of size (order x number of points x s) at once:
2.64 GB for the pair of involution type (15, 17) of §4.1 of the first paper (F = sqrt(2) A_15^+ and K of rank 17; s = 15, order 720) and 8.6 GB for s = 16 and order 1024.

The value of orbits_of is determined point by point (the least point of the orbit of the point and
the size of the orbit). Computing it on slices of the points and joining the results gives the same
array. to_rep, child, the stabilizers and the normal form in C are not touched, so the keys (the
bytes of canon) and |Stab| are the same as for rank32.gluing.Glue. The tests check that the sliced
values equal the original ones and that the keys equal those of rank32.gluing.

rank32.gluing の正規形を、キーを変えずに少ないメモリで作る。

rank32.gluing は変えない。ここでは同じ計算を、メモリの大きい箇所だけ分けて行う部分クラスを置く。

メモリを食っていたのは Glue.__init__ の F 側の鎖（Glue._fixed_chain）である。鎖の各ステージで、全部の
点 2^s について Node.orbits_of を呼ぶ。位数 <= SMALL のノードは元を全部持っていて、orbits_of は
(位数 x 点の数 x s) の int64 の配列を一度に作る。第 1 の論文の §4.1 の対合の型 (15, 17) の組（F = sqrt(2) A_15^+、K は階数 17。s = 15、位数 720）で 2.64 GB、
s = 16 で位数 1024 なら 8.6 GB になる。

orbits_of の値は点ごとに独立に決まる（その点の軌道の最小の点と、軌道の大きさ）。だから点を区切って
計算してつなげれば、同じ配列になる。to_rep、child、安定化部分群、C の正規形には触れないので、
キー（canon の bytes）と |Stab| は rank32.gluing.Glue と同じである。試験は、区切った値と元の値の一致と、
キーが rank32.gluing と一致することを確かめる。
"""
import numpy as np

from . import gluing as glue

CHUNK_BYTES = 1 << 27     # largest (order x points x s) int64 array built in one slice (128 MB) / 1 回の区切りで作る (位数 x 点 x s) の int64 の上限（128 MB）


def orbits_of_chunked(node, y, d, chunk_bytes=CHUNK_BYTES):
    """The same values as Node.orbits_of(y, d), computed on slices of the points when the node keeps all its elements.

    Node.orbits_of(y, d) と同じ値を、元を全部持つノードでは点を区切って計算する。"""
    if node.elems is None:
        return node.orbits_of(y, d)
    y = np.asarray(y, dtype=np.int64)
    s = node.tree.s
    order = len(node.elems)
    step = max(1, chunk_bytes // (8 * order * s))
    reps = np.empty(len(y), dtype=np.int64)
    sizes = np.empty(len(y), dtype=np.int64)
    for a in range(0, len(y), step):
        r, z = node.orbits_of(y[a:a + step], d)
        reps[a:a + step] = r
        sizes[a:a + step] = z
    return reps, sizes


class Glue(glue.Glue):
    """The same keys as rank32.gluing.Glue; the chain on the F side is built in slices.

    rank32.gluing.Glue と同じキー。F 側の鎖を区切って作る。"""

    def _fixed_chain(self, T):
        levels = []
        node = T.root; pts = np.arange(1 << self.s)
        while node.order > 1:
            l, sz = orbits_of_chunked(node, pts, 0)
            reps = pts[(l == pts) & (sz > 1)]
            b = int(reps[np.argmin(sz[reps])])
            orb = pts[l == b]
            us = np.stack([glue.inverse(node.to_rep(int(c), 0)) for c in orb])
            levels.append((orb, us))
            node = node.child(b, 0)
        return levels
