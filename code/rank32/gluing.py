"""Gluing of involutions: double cosets of the anti-isometries phi between the discriminant forms of 2-elementary parts F and K.

Notation: D_F = F^#/F and D_K = K^#/K are vector spaces over F_2; an element is an integer mask
(bit i is the coefficient of basis vector i). A linear map is the list of "the mask of the image
of basis vector i" (a numpy array of length s). bar G_F and bar G_K are the images of Aut F and
Aut K in the discriminant groups. (a, b) in bar G_F x bar G_K sends phi to b phi a^{-1}.

Normal form canon(): the least image of Gamma = {(x, phi(x)) : x in S_F} (S_F a bar G_F-invariant
set spanning D_F) under the direct product bar G_F x bar G_K.
We go down the tree of stabilizer subgroups (one Node of the chain for each sequence of
representatives of points); at each stage we branch to every element of the orbit with the least
number of hits. At a leaf the group of (a, b) is trivial, so the number of leaves that reach the
least form is |Stab(phi)| (in the image). Leaves with the same form along the chain of F are merged
into one, and their weights are added.
The C part _gluing.c is compiled on first use (gcc -O2 -shared -fPIC) into the directory
$RANK32_BUILD (default $XDG_CACHE_HOME/rank32, that is ~/.cache/rank32). The keys depend on the constant SMALL and on the C code;
neither may be changed.

対合の貼り合わせ：2-初等的な部品 F, K の判別式形式の間の反等長 phi の両側剰余類を扱う道具。

記法：D_F = F^#/F、D_K = K^#/K は F_2 上のベクトル空間で、元は整数の mask（ビット i が基底 i の
係数）である。線形写像は「基底 i の像の mask」の列（長さ s の numpy 配列）で持つ。bar G_F、bar G_K は
Aut F、Aut K の判別式群への像である。(a, b) in bar G_F x bar G_K は phi を b phi a^{-1} に移す。

正規形 canon()：Gamma = {(x, phi(x)) : x in S_F}（S_F は bar G_F で不変で D_F を張る集合）の、
直積群 bar G_F x bar G_K の下の最小の像。
安定化部分群の鎖（その Node。点の代表の列ごとに 1 つ）を降りながら、各ステージで「当たる数が
最小の点の軌道」の元すべてに枝を分ける。葉で (a, b) の群は自明になるので、最小の形に到達した
葉の数が |Stab(phi)|（像の中の）である。F の鎖の形が同じ葉は、重みを足して 1 つにまとめる。
C の部分 _gluing.c は、初めて使うときに $RANK32_BUILD（既定は $XDG_CACHE_HOME/rank32、つまり ~/.cache/rank32）へコンパイルする
（gcc -O2 -shared -fPIC）。キーは定数 SMALL と C のコードに依るので、どちらも変えてはいけない。
"""
import ctypes
import hashlib
import os
from collections import deque
import random
import subprocess
from pathlib import Path

import numpy as np
from sage.all import libgap


_CORE = None


def _build_dir():
    """The directory for compiled helpers.

    コンパイルしたものを置くディレクトリ。"""
    d = os.environ.get("RANK32_BUILD")
    if d:
        return Path(d)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "rank32"


def core():
    """Compile _gluing.c if needed and load it.

    _gluing.c を必要ならコンパイルして読み込む。"""
    global _CORE
    if _CORE is None:
        src = Path(__file__).with_name("_gluing.c")
        build = _build_dir()
        so = build / ("_gluing-" + hashlib.sha256(src.read_bytes()).hexdigest()[:12] + ".so")
        if not so.exists():
            build.mkdir(parents=True, exist_ok=True)
            tmp = so.with_name(so.stem + f".{os.getpid()}.so")
            subprocess.check_call(["gcc", "-O2", "-shared", "-fPIC", "-o", str(tmp), str(src)])
            os.replace(tmp, so)
        L = ctypes.CDLL(str(so))
        L.tree_new.restype = ctypes.c_void_p
        L.node_add.restype = ctypes.c_int
        L.mat_enumerate.restype = ctypes.c_long
        for name, args in dict(
            tree_new=[ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p],
            node_add=[ctypes.c_void_p, ctypes.c_uint64, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p],
            node_set_tables=[ctypes.c_void_p, ctypes.c_int, ctypes.c_int] + [ctypes.c_void_p] * 4,
            child_set=[ctypes.c_void_p, ctypes.c_int, ctypes.c_uint32, ctypes.c_int],
            set_stages=[ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p],
            set_chain=[ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p],
            canon=[ctypes.c_void_p] + [ctypes.c_void_p] * 4 + [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p],
            stats=[ctypes.c_void_p, ctypes.c_void_p],
            mat_apply64=[ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p],
            mat_invert64=[ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int],
            orbit_tables=[ctypes.c_int, ctypes.c_int, ctypes.c_void_p] + [ctypes.c_void_p] * 4,
            mat_enumerate=[ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_long, ctypes.c_void_p],
            verify_order=[ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_uint64, ctypes.c_uint64, ctypes.c_int],
        ).items():
            getattr(L, name).argtypes = args
        _CORE = L
    return _CORE


def _ptrs(arrs):
    return (ctypes.c_void_p * len(arrs))(*[a.ctypes.data for a in arrs])


def apply(A, m):
    """Apply the linear map A to the array of masks m.

    線形写像 A を mask の配列 m に作用させる。"""
    A = np.ascontiguousarray(A, dtype=np.int64); m = np.asarray(m, dtype=np.int64, order="C")   # m may be 0-dimensional / m は 0 次元でもよい
    out = np.empty(m.shape, dtype=np.int64)
    core().mat_apply64(A.ctypes.data, len(A), m.ctypes.data, m.size, out.ctypes.data)
    return out


def compose(A, B):
    """A o B (B first).

    A o B（先に B）。"""
    return apply(A, B)


def table(A):
    """The table of A on all elements (length 2^s).

    A の全元の表（長さ 2^s）。"""
    s = len(A)
    t = np.zeros(1 << s, dtype=np.int64)
    for l in range(s):
        t[1 << l: 2 << l] = t[0: 1 << l] ^ A[l]
    return t


def inverse(A):
    """The inverse of A.

    A の逆。"""
    A = np.ascontiguousarray(A, dtype=np.int64)
    R = np.empty(len(A), dtype=np.int64)
    assert core().mat_invert64(A.ctypes.data, R.ctypes.data, len(A)) == 0, "singular matrix"
    return R


SMALL = 1024   # nodes of order <= SMALL keep all elements and build no tables; the keys depend on this value
               # これ以下の位数のノードは元を全部持ち、表を作らない。キーはこの値に依る
MAXTABLES = 6000   # most nodes that keep tables (7 bytes per point per stage) / 表を持つノードの上限（ステージあたり 7 B/点）
NRAND = 16    # number of generators of a stabilizer (made from nearly uniform elements); the order is then confirmed (verify_order in C, or GAP)
              # 安定化部分群の生成元の数（一様に近い元から作る）。位数は C の verify_order か GAP で確かめる


class Node:
    """One stage of the stabilizer chain of a group; its children fix a representative.

    Large group: for each stage (a set of points, domain) a table point -> representative of its
    orbit and a tree spanned from the representatives (built only for the stages that are needed).
    Small group (order <= SMALL): the array elems (order x s) of all elements, applied to arrays of
    points when needed.
    No GAP group is kept (a permutation group on 16384 points takes 5 MB per node; 785 nodes reached 4 GB).

    群の安定化部分群の 1 ステージ。ツリーの中で代表を固定する子を持つ。

    大きい群：ステージ（点の集合 domain）ごとに、点 -> 軌道の代表 の表と、代表から張ったツリーを持つ（必要に
    なったステージだけ作る）。
    小さい群（order <= SMALL）：元の配列 elems（order x s）を持ち、点の配列に対して都度計算する。
    GAP の群は持たない（16384 点の置換群を残すと 1 ノード 5 MB になり、785 ノードで 4 GB に達した）。"""
    __slots__ = ("tree", "gens", "order", "arr", "children", "elems", "pool", "cid", "keep")

    def __init__(self, tree, gens, order, elems=None):
        self.tree = tree; self.gens = gens; self.order = order
        self.children = {}; self.elems = elems; self.arr = {}; self.pool = None
        self.cid = None; self.keep = []

    def tables(self, d):
        """The tables of stage d: (label, size, parent, pgen, invs, trans), indexed by the position in the domain.

        ステージ d の表：(label, size, parent, pgen, invs, trans)。label などは domain の添字で引く。"""
        r = self.arr.get(d)
        if r is not None:
            return r
        tree = self.tree
        dom = tree.domains[d]; idx = tree.index[d]; n = len(dom)
        tabs = [np.ascontiguousarray(idx[apply(g, dom)], dtype=np.int32) for g in self.gens]
        label = np.empty(n, dtype=np.uint16); size = np.empty(n, dtype=np.uint16); parent = np.empty(n, dtype=np.uint16); pgen = np.empty(n, dtype=np.uint8)
        core().orbit_tables(n, len(tabs), _ptrs(tabs), label.ctypes.data, size.ctypes.data, parent.ctypes.data, pgen.ctypes.data)
        r = (label, size, parent, pgen, [inverse(g) for g in self.gens], {})
        self.arr[d] = r
        tree.note_tables(self)
        return r

    def orbits_of(self, y, d):
        """(representative point of the orbit, size of the orbit) for an array y of points of the domain of stage d.

        点の配列 y（ステージ d の domain の点）の (軌道の代表の点, 軌道の大きさ)。"""
        if self.elems is not None:
            bits = (y[None, :, None] >> np.arange(self.tree.s)[None, None, :]) & 1      # 1 x n x s
            im = np.bitwise_xor.reduce(bits * self.elems[:, None, :], axis=2)            # order x n
            o = np.sort(im, axis=0)
            return o[0], 1 + (o[1:] != o[:-1]).sum(axis=0)
        label, size, *_ = self.tables(d)
        ii = self.tree.index[d][y]
        l = label[ii]
        return self.tree.domains[d][l], size[l].astype(np.int64)

    def to_rep(self, p, d):
        """An element that sends the point p to the representative (least point) of its orbit.

        点 p を軌道の代表（軌道の最小の点）に移す元。"""
        p = int(p)
        if self.elems is not None:
            bits = np.nonzero((p >> np.arange(self.tree.s)) & 1)[0]
            im = np.bitwise_xor.reduce(self.elems[:, bits], axis=1)
            return self.elems[int(np.argmin(im))]
        label, size, parent, pgen, invs, trans = self.tables(d)
        r = trans.get(p)
        if r is not None:
            return r
        A = self.tree.identity
        q = int(self.tree.index[d][p])
        while parent[q] != q:
            A = compose(invs[pgen[q]], A)
            q = int(parent[q])
        trans[p] = A
        return A

    def random_element(self, rng):
        """A nearly uniform random element by product replacement.

        積の置き換え（product replacement）で一様に近い元を返す。"""
        if self.pool is None:
            g = list(self.gens)
            while len(g) < 10:
                g.append(g[len(g) % len(self.gens)])
            self.pool = g
            for _ in range(60):
                self._prstep(rng)
        return self._prstep(rng)

    def _prstep(self, rng):
        g = self.pool; i = rng.randrange(len(g)); j = rng.randrange(len(g) - 1); j += j >= i
        x = g[j] if rng.random() < .5 else inverse(g[j])
        g[i] = compose(g[i], x) if rng.random() < .5 else compose(x, g[i])
        return g[i]

    def child(self, rep, d):
        rep = int(rep)
        c = self.children.get(rep)
        if c is None:
            c = self.tree.stabilizer(self, rep, d); self.children[rep] = c
        return c


class Tree:
    """The tree of stabilizer subgroups, one for each sequence of representatives, of a group acting by linear maps (an image in the discriminant group).

    domains are the sets of points of the stages (arrays of masks), kept by the group (on the K side,
    the sets determined by the values of the quadratic form on the images of the stages of F).

    線形写像で作用する群（判別式の像）の、代表の列ごとの安定化部分群のツリー。

    domains はステージごとの点の集合（mask の配列）で、群はこれらを保つ（K 側では、F のステージの像の 2 次形式の
    値の条件で決まる集合）。"""

    def __init__(self, s, gens, order, domains, seed=1):
        self.s = s
        self.identity = 1 << np.arange(s, dtype=np.int64)
        self.domains = [np.asarray(d, dtype=np.int64) for d in domains]
        self.index = []
        for d in self.domains:
            ix = np.full(1 << s, -1, dtype=np.int32); ix[d] = np.arange(len(d)); self.index.append(ix)
        self.rng = random.Random(seed); self.faithful = None
        self.root = self._node([np.asarray(g, dtype=np.int64) for g in gens], order)
        self.nnodes = 1; self.nverified = 0
        self.c = None; self.byid = []; self.table_nodes = deque()

    def note_tables(self, node):
        """Bound the number of nodes that keep tables (drop the oldest; rebuild when needed); only once the tree is registered in C.

        表を持つノードの数を抑える（古いものから捨てる。要れば作り直す）。ツリーが C に登録済みのときだけ。"""
        if self.c is None:
            return
        self.table_nodes.append(node)
        while len(self.table_nodes) > MAXTABLES:
            old = self.table_nodes.popleft()
            if old is self.root or old is node:
                self.table_nodes.append(old)
                continue
            if old.cid is not None:
                for d in old.arr:
                    core().node_set_tables(self.c, old.cid, d, None, None, None, None)
            old.arr.clear()

    def _node(self, gens, order):
        if order == 1:
            return Node(self, [], 1)
        if order <= SMALL:
            return Node(self, gens, order, elems=self._enumerate(gens, order))
        return Node(self, gens, order)

    def _enumerate(self, gens, order):
        """List all elements (an order x s array); AssertionError if the generators are too few.

        元を全部並べる（order x s の配列）。生成元が足りなければ AssertionError。"""
        G = np.ascontiguousarray(np.stack(gens), dtype=np.int64)
        X = np.empty((order, self.s), dtype=np.int64)
        n = core().mat_enumerate(G.ctypes.data, len(G), self.s, order, X.ctypes.data)
        assert n == order, ("too few generators", n, order)
        return X

    def _generates(self, gens, order, last=False):
        """Whether the group of gens has order >= order. A yes from the randomised Schreier-Sims in C is certain (verify_order in _gluing.c).

        When it is not reached, the generators are usually really too few, so the caller adds
        generators and tries again. Only the last try is decided by GAP.

        gens の群の位数が order 以上か。C の乱択 Schreier-Sims は「はい」なら確実（_gluing.c の verify_order）。

        届かないときはたいてい生成元が本当に足りないので、呼び出し側が生成元を足して試し直す。
        最後の試行だけ GAP で決める。"""
        G = np.ascontiguousarray(np.stack(gens), dtype=np.int64)
        if core().verify_order(G.ctypes.data, len(G), self.s, order, self.rng.getrandbits(63), 40):
            return True
        return last and self._gap_order(gens) == order

    def _gap_order(self, gens):
        """The order of the group of gens by GAP, acting on a small invariant set faithful that spans the discriminant group (hence faithfully); on all points one call takes 0.5 seconds.

        元の位数を GAP で確かめる。点は判別式を張る小さい不変な集合 faithful（張れば忠実）。全体だと 1 回 0.5 秒かかる。"""
        if self.faithful is None:
            self.faithful = invariant_span(self.root, len(self.domains) - 1, self.s, self.domains[-1])
            self.findex = np.full(1 << self.s, -1, dtype=np.int32); self.findex[self.faithful] = np.arange(len(self.faithful))
        H = libgap.Group([libgap.PermList([int(x) + 1 for x in self.findex[apply(g, self.faithful)]]) for g in gens])
        return int(H.Size())

    def stabilizer(self, node, rep, d):
        """The child node: the stabilizer of the point rep of stage d in node.

        子のノード：node の中でステージ d の点 rep を固定する安定化部分群。"""
        self.nnodes += 1
        if node.elems is not None:
            bits = np.nonzero((rep >> np.arange(self.s)) & 1)[0]
            im = np.bitwise_xor.reduce(node.elems[:, bits], axis=1)
            keep = node.elems[im == rep]
            return Node(self, [], 1) if len(keep) == 1 else Node(self, [], len(keep), elems=keep)
        label, size, *_ = node.tables(d)
        order = node.order // int(size[self.index[d][rep]])
        if order == 1:
            return Node(self, [], 1)
        node.tables(d)[5].clear()
        # Schreier generators to_rep(h(rep)) o h from nearly uniform elements h
        # 一様に近い元 h から Schreier 生成元 to_rep(h(rep)) o h を作る
        gens = []; seen = set()
        tries = 0
        while True:
            tries += 1
            for _ in range(NRAND):
                h = node.random_element(self.rng)
                p = int(apply(h, np.array([rep], dtype=np.int64))[0])
                x = compose(node.to_rep(p, d), h)
                k = x.tobytes()
                if k not in seen and not (x == self.identity).all():
                    seen.add(k); gens.append(x)
            if order <= SMALL:
                try:
                    return Node(self, gens, order, elems=self._enumerate(gens, order))
                except AssertionError:
                    continue
            if self._generates(gens, order, last=tries >= 5):
                self.nverified += 1
                return Node(self, gens, order)
            assert tries < 5, "too few generators of the stabilizer"


def invariant_span(root, d, s, dom, split=False):
    """A set spanning the discriminant group: orbits of the group at node root, smallest first, added when they raise the rank; with split, the list of orbits (fixed points left out).

    ルートのノード root の群の軌道のうち小さいものから、階数が増えるものを足して判別式を張る集合を作る。split なら軌道ごとの列（固定点は除く）。"""
    label, size, *_ = root.tables(d)
    reps = [int(r) for r in np.nonzero((label == np.arange(len(dom))) & (dom > 0))[0]]
    reps.sort(key=lambda r: (size[r], r))
    basis = {}; chosen = []
    def reduce(x):
        while x:
            h = x.bit_length() - 1
            if h not in basis:
                return x
            x ^= basis[h]
        return 0
    for r in reps:
        mem = dom[label == r]
        added = False
        for m in mem:
            x = reduce(int(m))
            if x:
                basis[x.bit_length() - 1] = x; added = True
        if added and (size[r] > 1 or not split):
            chosen.append(mem)
        if len(basis) == s:
            break
    assert len(basis) == s
    return chosen if split else np.concatenate(chosen)


class Glue:
    """A pair of parts. F and K are dict(k, q2, mu2, B4, bits) (tables of the discriminant forms); gF and gK are generators acting on the discriminant groups.

    部品の組。F, K は dict(k, q2, mu2, B4, bits)（判別式の表）、gF, gK は判別式の上の生成元。"""

    def __init__(self, F, K, gF, gK, orderF, orderK):
        assert F["k"] == K["k"]
        self.F, self.K = F, K; self.s = F["k"]
        n = 1 << self.s
        self.gF = [np.asarray(g, dtype=np.int64) for g in gF]
        self.TF = Tree(self.s, gF, orderF, [np.arange(n)])
        self.stages = self._invariant_span(self.TF)
        # domains of the stages on the K side: points whose value 2q can be the image of a stage of F
        # K 側のステージの domain：F のステージの像が入る 2q の値の点
        doms = []
        for st in self.stages:
            vals = np.unique((-F["q2"][st]) % 4)
            doms.append(np.nonzero(np.isin(K["q2"], vals) & (np.arange(n) > 0))[0])
        doms.append(np.arange(n))
        self.TK = Tree(self.s, gK, orderK, doms)
        self.Q1F = np.nonzero(F["q2"] == 2)[0]
        self.chainF = self._fixed_chain(self.TF)

    def _fixed_chain(self, T):
        """The chain on the F side: at each stage, the base point is the representative of a smallest orbit among the moved points.

        F 側の鎖：各ステージで、動く点のうち軌道の小さいものの代表を基底にする。"""
        levels = []
        node = T.root; pts = np.arange(1 << self.s)
        while node.order > 1:
            l, sz = node.orbits_of(pts, 0)
            reps = pts[(l == pts) & (sz > 1)]
            b = int(reps[np.argmin(sz[reps])])
            orb = pts[l == b]
            us = np.stack([inverse(node.to_rep(int(c), 0)) for c in orb])
            levels.append((orb, us))
            node = node.child(b, 0)
        return levels

    def _invariant_span(self, T):
        """Span D by adding orbits, smallest first, that raise the rank; returns the list of stages (bar G-invariant sets of points).

        小さい軌道から順に、階数が増えるものを足して D を張る。ステージ（bar G 不変な点の集合）の列を返す。"""
        return invariant_span(T.root, 0, self.s, np.arange(1 << self.s), split=True)

    # --- basic operations on phi / phi の基本操作 ---------------------------------
    def is_antiisometry(self, phi):
        """Whether phi is an anti-isometry D_F -> D_K.

        phi が反等長 D_F -> D_K か。"""
        t = table(phi)
        return len(np.unique(t)) == 1 << self.s and bool((self.K["q2"][t] == (-self.F["q2"]) % 4).all())

    def is_extremal(self, phi):
        """Whether the glued lattice has minimum >= 4 (every nonzero class has 2 mu_F + 2 mu_K >= 8).

        貼り合わせた格子の最小ノルムが 4 以上か（0 でない類のすべてで 2 mu_F + 2 mu_K >= 8）。"""
        t = table(phi)
        return bool((self.F["mu2"][1:] + self.K["mu2"][t[1:]] >= 8).all())

    def reflection(self, v):
        """The reflection t_v(x) = x + 2b(x, v) v of D_F (2q(v) = 2 mod 4).

        D_F の反射 t_v(x) = x + 2b(x, v) v（2q(v) = 2 mod 4）。"""
        F = self.F; s = self.s
        vb = (v >> np.arange(s)) & 1
        e = ((F["B4"] @ vb) % 4) // 2      # 4b with basis vector i is 0 or 2 mod 4 / 基底 i との 4b mod 4 は 0 か 2
        return (1 << np.arange(s, dtype=np.int64)) ^ (e * v)

    def orbit_reflection(self, a, v):
        """The product of the reflections in the orbit of v under a in O(q_F), when the orbit is pairwise orthogonal; otherwise None.

        The product commutes with a, so (a, b) in Stab(phi) implies (a, b) in Stab(phi o product).
        Classes with a large |Stab| are found easily by this move.

        a in O(q_F) の下の v の軌道の反射の積（軌道が互いに直交するときだけ。それ以外は None）。

        この積は a と可換なので、(a, b) in Stab(phi) なら (a, b) in Stab(phi o 積)。|Stab| の大きい類は、
        この動きで探すと当たりやすい。"""
        F = self.F; s = self.s
        orb = [v]; x = int(apply(a, v))
        while x != v:
            orb.append(x); x = int(apply(a, x))
        bits = [(w >> np.arange(s)) & 1 for w in orb]
        for i in range(len(orb)):
            for j in range(i):
                if (bits[i] @ F["B4"] @ bits[j]) % 4:
                    return None
        t = self.reflection(orb[0])
        for w in orb[1:]:
            t = compose(t, self.reflection(w))
        return t

    def stab_prime_element(self, phi, rng):
        """A random element of prime order of Stab(phi) on the F side (None if Stab is trivial).

        Stab(phi) の F 側の素数位数の元を乱択で 1 つ（Stab が自明なら None）。"""
        _, n, els = self.canon(phi, want_stab=True)
        one = 1 << np.arange(self.s, dtype=np.int64)
        if n == 1 or not len(els):
            return None
        while True:
            a = els[rng.randrange(len(els))]
            for _ in range(rng.randrange(3)):
                a = compose(a, els[rng.randrange(len(els))])
            pows = [one]
            while True:
                x = compose(a, pows[-1])
                if np.array_equal(x, one):
                    break
                pows.append(x)
            m = len(pows)
            if m > 1:
                p = next(q for q in range(2, m + 1) if m % q == 0)
                return pows[m // p]

    # --- normal form / 正規形 ----------------------------------------------------
    def _setup_core(self):
        L = core(); T = self.TK
        self._dom32 = [d.astype(np.uint32) for d in T.domains]
        self._idx = T.index
        T.c = L.tree_new(self.s, len(T.domains), _ptrs(self._dom32), _ptrs(self._idx))
        self._st32 = [x.astype(np.uint32) for x in self.stages]
        self._stn = np.array([len(x) for x in self.stages], dtype=np.int32)
        L.set_stages(T.c, len(self.stages), _ptrs(self._st32), self._stn.ctypes.data)
        self._lorb = [o.astype(np.uint32) for o, _ in self.chainF]; self._lus = [u.astype(np.uint32) for _, u in self.chainF]
        self._ln = np.array([len(o) for o, _ in self.chainF], dtype=np.int32)
        L.set_chain(T.c, len(self.chainF), _ptrs(self._lorb), _ptrs(self._lus), self._ln.ctypes.data)
        self._register(T.root)
        self._need = np.zeros(4, dtype=np.uint32); self._form = np.zeros(self.s, dtype=np.uint32)
        self._stab = np.zeros(1, dtype=np.uint64); self._nel = np.zeros(1, dtype=np.int32)
        self._els = np.zeros((64, self.s), dtype=np.uint32)

    def _register(self, node):
        L = core(); T = self.TK
        if node.order == 1 or node.elems is not None:
            e = None
            if node.elems is not None:
                e = node.elems.astype(np.uint32); node.elems = e       # small nodes on the K side get no children (canon stops there), so no int64 copy is kept
                # K 側の小さいノードの子は作らない（canon がそこで止まる）ので int64 の複製は持たない
            node.cid = L.node_add(T.c, int(node.order), 1 if e is not None else 0, 0 if e is None else len(e), None if e is None else e.ctypes.data, None)
        else:
            inv = np.stack([inverse(g) for g in node.gens]).astype(np.uint32); node.keep.append(inv)
            node.cid = L.node_add(T.c, int(node.order), 0, 0, None, inv.ctypes.data)
        T.byid.append(node)
        assert T.byid[node.cid] is node

    def canon(self, phi, want_stab=False):
        """Return (bytes of the normal form psi, |Stab| in bar G_F x bar G_K [, elements of Stab on the F side (at most 64; they need not generate Stab)]).

        Go down the tree on the K side: at each stage keep only the branches whose table of hits per
        orbit of Y = psi(points of the stage) is least, and branch to every element of the orbit with
        the fewest hits. When the group of K is exhausted, keep the leaves with the least set Y and
        make psi least along the fixed chain on the F side. For (a, b) in Stab, b determines a, so the
        number of leaves reaching the least form is |Stab|.

        (正規形 psi の bytes, |Stab|（bar G_F x bar G_K の中の）[, Stab の F 側の元の列（最大 64 個。Stab 全体を生成するとは限らない）]) を返す。

        K 側のツリーを降りる：各ステージで Y = psi(ステージの点) の軌道ごとの当たる数の表が最小の枝だけ残し、当たる数が
        最小の軌道の元すべてに分かれる。K の群が尽きたら、葉を Y の集合の最小で絞り、F 側の固定の鎖で
        psi を最小にする。(a, b) in Stab では b が a を決めるので、最小の形に到達した葉の数が |Stab|。"""
        L = core(); T = self.TK
        if T.c is None:
            self._setup_core()
        phi = np.asarray(phi, dtype=np.int64)
        p32 = phi.astype(np.uint32)
        while True:
            rc = L.canon(T.c, p32.ctypes.data, self._form.ctypes.data, self._stab.ctypes.data, self._need.ctypes.data,
                         self._els.ctypes.data if want_stab else None, len(self._els), self._nel.ctypes.data)
            if rc == 0:
                break
            if rc == 1:
                node = T.byid[int(self._need[0])]; rep = int(self._need[1]); d = int(self._need[2])
                c = node.child(rep, d)
                if c.cid is None:
                    self._register(c)
                L.child_set(T.c, node.cid, rep, c.cid)
            elif rc == 2:
                node = T.byid[int(self._need[0])]; d = int(self._need[1])
                label, size, parent, pgen, *_ = node.tables(d)
                L.node_set_tables(T.c, node.cid, d, label.ctypes.data, size.ctypes.data, parent.ctypes.data, pgen.ctypes.data)
            else:
                raise AssertionError("the group of K does not act faithfully on the discriminant group")
        form = self._form.tobytes(); n = int(self._stab[0])
        if not want_stab:
            return form, n
        ne = int(self._nel[0]); bs = self._els[:ne].astype(np.int64)
        phinv = inverse(phi)
        els = [compose(phinv, compose(b, phi)) for b in bs]
        return form, n, els
