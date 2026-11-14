"""Automorphism group and isometry test of an extremal even unimodular lattice of rank 32 from its shortest vectors.

This is a second tool for |Aut L|, written independently of rank32.canonical (it builds a
different graph). It works on the 73440 lines (pairs +-x) of the 146880 vectors of norm 4.
  1. Profile: the lines are told apart by their classes b mod 2 (they must be distinct).
     For each line an exact integer invariant is computed by the C engine (moments of the
     numbers of pairs of lines with the same sum mod 2). The rarest value of the invariant
     is the key of the lattice; the lines with that value are the candidates, and the first
     candidate is the anchor.
  2. Local graph: the anchor x together with the 1240 norm-4 vectors y with (x, y) = 2 gives
     a local Gram matrix of size 1241. bliss returns its canonical labelling and the
     generators of its automorphism group.
  3. Stabilizer: a graph automorphism lifts to a rational isometry; the ones that keep the
     lattice form the stabilizer of the anchor vector x (with -1 they generate the stabilizer of
     its line). They are found as the set stabilizer
     of the glue in the discriminant group M^#/M of the sublattice M spanned by the local
     vectors (GAP), or, when M^#/M is not 2-elementary of rank <= 16, by a Schreier orbit
     of lattices.
  4. Orbit: every candidate line is compared with the anchor by the canonical local Gram
     matrix; a match that lifts to an integral isometry adds a generator. Then
     |Aut L| = 2 * |orbit of the anchor line| * |stabilizer of x| (the 2 is -1).
Every generator g is checked exactly: g^T G g = G (column convention).
Isometry: two lattices with different keys are not isometric; otherwise the anchor of the
first is compared with every candidate of the second, and a map T with T^T G_1 T = G_2 is
returned and checked.

The C engine _bliss_automorphisms.c is compiled on first use (gcc -O3 -fopenmp, linked to
bliss) into $RANK32_BUILD (default $XDG_CACHE_HOME/rank32, that is ~/.cache/rank32). bliss is looked up under
$BLISS_PREFIX, else $CONDA_PREFIX, else the prefix of the running Python. Temporary files
go to the directory OUT of this module (default: the system temporary directory); callers
may replace the attributes OUT and PROGRESS and the functions progress() and binary(). Progress is written to the file
PROGRESS when it is set (default: not written).

ランク 32 の極値偶ユニモジュラー格子の自己同型群と同型判定を、最短ベクトルから求める。

|Aut L| を求める二つ目の道具で、rank32.canonical とは独立に書かれている（作るグラフが違う）。
ノルム 4 のベクトル 146880 本の直線（+-x の組）73440 本の上で働く。
  1. プロファイル：直線を法 2 の類 b で区別する（類は互いに異なる必要がある）。
     直線ごとに、C のエンジンで整数の厳密な不変量を求める（法 2 の和が同じになる直線の組の
     数のモーメント）。最も少ない値を格子のキーとする。その値を持つ直線が候補で、最初の候補が
     アンカーである。
  2. 局所のグラフ：アンカー x と、(x, y) = 2 となるノルム 4 のベクトル y の 1240 本から、
     大きさ 1241 の局所の Gram 行列を作る。bliss がその標準ラベルと自己同型群の生成元を返す。
  3. 安定化群：グラフの自己同型は有理数の等長写像に持ち上がる。そのうち格子を保つものが、
     アンカーのベクトル x の安定化群である（-1 と合わせると、その直線の安定化群になる）。局所のベクトルが張る部分格子 M の判別群 M^#/M の中で、
     貼り合わせの集合の安定化群として求める（GAP）。M^#/M が階数 16 以下の 2-初等でない
     ときは、格子の Schreier 軌道で求める。
  4. 軌道：候補の直線をそれぞれ、標準の局所 Gram 行列でアンカーと比べる。一致して整数の
     等長写像に持ち上がれば、生成元が一つ増える。|Aut L| = 2 * |アンカーの直線の軌道| *
     |x の安定化群|（2 は -1 の分）。
生成元 g はすべて厳密に検算する：g^T G g = G（列の約束）。
同型：キーが違えば同型でない。同じなら、一つ目のアンカーを二つ目の候補の全部と比べ、
T^T G_1 T = G_2 となる T を返して検算する。

C のエンジン _bliss_automorphisms.c は、初めて使うときに $RANK32_BUILD（既定は $XDG_CACHE_HOME/rank32、つまり ~/.cache/rank32）へ
コンパイルする（gcc -O3 -fopenmp、bliss をリンクする）。bliss は $BLISS_PREFIX、なければ
$CONDA_PREFIX、なければ実行中の Python の prefix の下で探す。一時ファイルは、このモジュールの
OUT（既定はシステムの一時ディレクトリ）に置く。呼び出し側は、属性 OUT・PROGRESS と、
関数 progress()・binary() を差し替えてよい。進捗は PROGRESS が設定されているときだけ、そのファイルに書く（既定は書かない）。
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time

import numpy as np
from sage.all import ZZ, QQ, GF, matrix, pari, Permutation, PermutationGroup
from sage.libs.gap.libgap import libgap

SOURCE = Path(__file__).with_name('_bliss_automorphisms.c')
OUT = None          # directory for temporary files; None = the system default / 一時ファイルの置き場。None はシステムの既定
PROGRESS = None     # file for the progress summary; None = not written / 進捗の要約のファイル。None なら書かない
START = time.monotonic()
STATE = dict(stage='starting', status='running', processed=0, total=None,
             unit='lattice', failures=0, unresolved=0)
PROGRESS_LOCK = threading.Lock()


def _atomic(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False) + '\n')
    tmp.replace(path)


def progress(**kw):
    """Update the progress state, and write it to PROGRESS when that is set.

    進捗の状態を更新し、PROGRESS が設定されていれば、そのファイルに書く。"""
    with PROGRESS_LOCK:
        STATE.update(kw)
        if PROGRESS is not None:
            _atomic(Path(PROGRESS), dict(STATE,
                    updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                    elapsed_s=round(time.monotonic()-START, 1), heartbeat=True))


def _build_dir():
    d = os.environ.get('RANK32_BUILD')
    if d:
        return Path(d)
    return Path(os.environ.get('XDG_CACHE_HOME', Path.home() / '.cache')) / 'rank32'


def _bliss_prefix():
    return Path(os.environ.get('BLISS_PREFIX') or os.environ.get('CONDA_PREFIX') or sys.prefix)


def binary():
    """The compiled engine, built once per version of the C source.

    コンパイルしたエンジン。C のソースの版ごとに 1 回だけ作る。"""
    build = _build_dir()
    exe = build / ('bliss-automorphisms-' + hashlib.sha256(SOURCE.read_bytes()).hexdigest()[:12])
    if not exe.exists():
        env = _bliss_prefix()
        build.mkdir(parents=True, exist_ok=True)
        tmp = exe.with_name(exe.name + f'.{os.getpid()}.tmp')     # parallel callers do not collide / 並んで呼んでも取り合わない
        subprocess.run(['gcc', '-O3', '-fopenmp', str(SOURCE), '-o', str(tmp),
                        '-I'+str(env/'include'), '-L'+str(env/'lib'),
                        '-Wl,-rpath,'+str(env/'lib'), '-lbliss', '-lstdc++'], check=True)
        os.replace(tmp, exe)
    return exe


def rows(M):
    """The rows of a Sage matrix as lists of Python integers.

    Sage の行列の行を、Python の整数のリストにする。"""
    return [[int(x) for x in r] for r in M.rows()]


class Shell:
    """The lines of norm-4 vectors of a lattice, their classes mod 2 and the profile key.

    格子のノルム 4 のベクトルの直線と、その法 2 の類、プロファイルのキー。
    gram: the Gram matrix; tag: a label for the progress; threads: OpenMP threads of the engine;
    test=True allows other ranks (all shortest vectors must still have norm 4).
    gram は Gram 行列、tag は進捗のタグ、threads はエンジンの OpenMP のスレッド数。
    test=True なら 32 以外の階数も許す（最短ベクトルがすべてノルム 4 であることは要る）。
    """
    def __init__(self, gram, tag, threads=6, test=False):
        self.G = matrix(ZZ, gram)
        G = self.G
        if G != G.T or not G.is_positive_definite():
            raise ValueError('symmetric positive definite Gram required')
        if not test and (G.nrows()!=32 or G.det()!=1 or any(x%2 for x in G.diagonal())):
            raise ValueError('32-dimensional even unimodular Gram required')
        progress(stage='short_vectors', current=tag)
        r = pari(G).qfminim(4)
        X = np.array(rows(matrix(ZZ, r[2]).T), dtype=np.int64)
        self.Gn = np.array(rows(G), dtype=np.int64)
        # Guard integer overflow before any numpy dot product. / numpy の積の前に整数のあふれを防ぐ
        if max(abs(int(x)) for x in X.flat)**2*max(abs(int(x)) for x in self.Gn.flat)*G.nrows()**2 >= 2**62:
            raise ValueError('basis coefficients too large; reduce the Gram first')
        norms = np.sum((X @ self.Gn)*X, axis=1)
        if not test and (int(r[0])!=146880 or len(X)!=73440 or not np.all(norms==4)):
            raise ValueError('input is not extremal: shell must consist of 146880 norm-4 vectors')
        if G.nrows()>32 or not np.all(norms==4):
            raise ValueError('norm-4 shell required')
        self.X = X
        self.b = ((X&1).astype(np.uint64)<<np.arange(G.nrows(),dtype=np.uint64)).sum(1).astype(np.uint32)
        if len(np.unique(self.b))!=len(X):
            raise ValueError('minimal lines must have distinct mod-2 classes')
        progress(stage='line_profiles')
        with tempfile.TemporaryDirectory(dir=OUT) as td:
            p = Path(td)/'profile'
            with p.open('wb') as f:
                f.write(np.uint32(len(X)).tobytes());f.write(self.b.tobytes())
            subprocess.run([str(binary()),'profile',str(p),str(p)+'.out'], check=True,
                           env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
            self.sig = np.fromfile(str(p)+'.out',dtype=np.uint64).reshape(-1,3)
        unique, counts = np.unique(self.sig,axis=0,return_counts=True)
        self.key = tuple(unique[np.lexsort((unique[:,2],unique[:,1],unique[:,0],counts))[0]])
        self.candidates = np.flatnonzero(np.all(self.sig==self.key,axis=1)).tolist()
        self.anchor = self.candidates[0]
        self.sort = np.argsort(self.b)

    def local(self, index):
        """The local Gram matrix at the line `index`, its canonical order and the bliss generators.

        直線 index での局所の Gram 行列、その標準の順序、bliss の生成元。"""
        x = self.X[index]
        ips = self.X @ (self.Gn @ x)
        idx = np.flatnonzero(np.abs(ips)==2)
        V = np.vstack([x, self.X[idx]*np.sign(ips[idx,None])])
        if len(V)!=1241 and self.G.nrows()==32:
            raise ValueError('unexpected neighborhood size')
        Q = V @ self.Gn @ V.T
        if np.max(np.abs(Q))>4:
            raise ValueError('Gram overflow')
        progress(stage='local_graph', anchor=int(index), local_vertices=len(V))
        with tempfile.TemporaryDirectory(dir=OUT) as td:
            p = Path(td)/'graph'
            with p.open('wb') as f:
                f.write(np.uint32(len(V)).tobytes());f.write(Q.astype(np.int8).tobytes())
            result = subprocess.run([str(binary()),'graph',str(p),str(p)],check=True,capture_output=True,text=True)
            stats = json.loads(result.stdout)
            lab = np.fromfile(str(p)+'.lab',dtype=np.uint32)
            gens = np.fromfile(str(p)+'.gens',dtype=np.uint32).reshape(-1,len(V))
        canonical = np.argsort(lab)
        # Complete Gram in canonical core order is the exact comparison, not a hash.
        # 標準の順序に並べた Gram 行列の全体で比べる（ハッシュではない）
        return dict(V=V, Q=Q, canonical=canonical, gens=gens, stats=stats)

    def permutation(self, R):
        """The permutation of the lines induced by the row action R.

        行の作用 R が直線に引き起こす置換。"""
        # R acts on row coordinates. Compute its action modulo 2 without large dots.
        # R は行の座標に作用する。大きな積を使わずに、法 2 での作用を求める
        bits = np.array(rows(R),dtype=np.int64)&1
        masks = (bits.astype(np.uint64)<<np.arange(bits.shape[1],dtype=np.uint64)).sum(1).astype(np.uint32)
        transformed = np.zeros(len(self.b),dtype=np.uint32)
        for a in range(4):
            table = np.zeros(256,dtype=np.uint32)
            for c in range(256):
                for j in range(8):
                    if c>>j&1 and 8*a+j<len(masks):table[c]^=masks[8*a+j]
            transformed ^= table[(self.b>>(8*a))&255]
        positions = np.searchsorted(self.b[self.sort],transformed)
        if np.any(positions>=len(self.b)) or not np.array_equal(self.b[self.sort][positions],transformed):
            raise ValueError('isometry does not permute shortest lines')
        return self.sort[positions]


def lift(A, B, G, H):
    """Return the exact rational row map R with R H R^t=G.

    R H R^t = G となる有理数の行の写像 R を厳密に返す。"""
    A = matrix(ZZ,A.tolist());B = matrix(ZZ,B.tolist())
    indices = A.T.pivots()
    if len(indices)!=G.nrows():
        raise ValueError('local vectors do not span rational space')
    R = matrix(QQ,A[list(indices)]).inverse()*B[list(indices)]
    if A*R!=B or R*H*R.T!=G:
        raise ValueError('graph map did not lift to rational isometry')
    return R


def integral(R):
    """R as an integer matrix if it is integral with determinant +-1, else None.

    R が整数で行列式が +-1 なら整数の行列として返し、そうでなければ None。"""
    if any(x.denominator()!=1 for x in R.list()):return None
    R=matrix(ZZ,R)
    return R if abs(R.det())==1 else None


def lattice_key(R):
    """Exact canonical row lattice, including denominators (not saturation).

    行が張る格子の厳密な標準形（分母を含む。飽和ではない）。"""
    d=R.denominator()
    H=matrix(ZZ,d*R).hermite_form()
    return tuple(QQ(x)/d for x in H.list())


def integral_stabilizer_fallback(a,G):
    """Schreier orbit of Z^n under the finite rational local graph group.

    A graph symmetry can preserve the index-two neighbor sublattice but not L.
    Stabilizing the lattice, rather than discarding bad generators, is essential.

    有限の有理数の局所グラフの群による Z^n の Schreier 軌道。
    グラフの対称性は、指数 2 の近傍の部分格子を保っても L を保つとは限らない。
    悪い生成元を捨てるのではなく、格子の安定化群を取ることが肝心である。
    """
    rational=[lift(a['V'],a['V'][p],G,G) for p in a['gens']]
    K=PermutationGroup([Permutation((p+1).tolist()) for p in a['gens']]) if len(a['gens']) else None
    identity=matrix.identity(QQ,G.nrows())
    representatives={lattice_key(identity):identity};todo=[identity]
    generators=[];seen=set()
    while todo:
        A=todo.pop()
        for R in rational:
            B=A*R;key=lattice_key(B)
            if key not in representatives:
                if len(representatives)>=1024:
                    raise ValueError('non-elementary discriminant fallback orbit exceeds 1024')
                representatives[key]=B;todo.append(B)
                progress(stage='lattice_orbit',lattice_orbit=len(representatives))
            S=integral(B*representatives[key].inverse())
            if S is None:raise ValueError('Schreier lattice stabilizer failure')
            k=tuple(S.list())
            if S!=identity and k not in seen:
                generators.append(S);seen.add(k)
    order=ZZ(K.order() if K else 1)
    if order%len(representatives):raise ValueError('nonintegral stabilizer order')
    return generators,order//len(representatives),representatives


def integral_stabilizer(a,G):
    """Set stabilizer of the glue in the finite discriminant group M*/M.

    For a 2-elementary M, encode M*/M by at most 16 binary coordinates;
    GAP set stabilizers avoid storing tens of thousands of rational lattices.

    有限の判別群 M*/M の中の、貼り合わせの集合の安定化群。
    M が 2-初等なら、M*/M を 16 個以下の 2 進の座標で表す。GAP の集合の安定化群を使えば、
    数万個の有理数の格子を持たずに済む。
    """
    V=matrix(ZZ,a['V'].tolist())
    M=matrix(ZZ,pari(V.T).mathnf()).T
    Mi=M.inverse();Q=M*G*M.T;D=2*Q.inverse()
    if any(x.denominator()!=1 for x in D.list()) or any(x.denominator()!=1 for x in (2*Mi).list()):
        return integral_stabilizer_fallback(a,G)
    B=matrix(GF(2),D).row_space().basis_matrix()
    rank=B.nrows();pivots=list(B.pivots());size=1<<rank;n=len(a['V'])
    if rank>16:return integral_stabilizer_fallback(a,G)
    vectors=((np.arange(size)[:,None]>>np.arange(rank))&1).astype(np.uint32)
    powers=1<<np.arange(rank,dtype=np.uint32)
    rational=[lift(a['V'],a['V'][p],G,G) for p in a['gens']]
    combined=[]
    for p,R in zip(a['gens'],rational):
        A=matrix(ZZ,M*R*Mi)
        T=(B*matrix(GF(2),A))[:,pivots]
        action=((vectors@np.array(rows(T),dtype=np.uint32))&1)@powers
        combined.append(Permutation(np.concatenate([p+1,n+action+1]).astype(int).tolist()))
    K=PermutationGroup(combined) if combined else PermutationGroup([Permutation(list(range(1,n+size+1)))])
    def subset(R):
        C=2*R*Mi
        if any(x.denominator()!=1 for x in C.list()):return None
        C=matrix(GF(2),C)
        # Membership in M*/M, not just projection to pivot columns.
        # M*/M に属することを確かめる（ピボットの列への射影だけでなく）
        if C[:,pivots]*B!=C:return None
        basis=C[:,pivots].row_space().basis_matrix()
        masks=[sum(int(row[j])<<j for j in range(rank)) for row in basis.rows()]
        values=[0]
        for mask in masks:values += [v^mask for v in values]
        return sorted(n+v+1 for v in values)
    base=subset(matrix.identity(QQ,G.nrows()))
    progress(stage='glue_stabilizer',discriminant_size=size,glue_size=len(base))
    group=K.gap();base=libgap(base)
    H=libgap.Stabilizer(group,base,libgap.OnSets)
    generators=[]
    for p in libgap.GeneratorsOfGroup(H):
        perm=np.array([int(x)-1 for x in libgap.ListPerm(p,n)],dtype=int)[:n]
        R=integral(lift(a['V'],a['V'][perm],G,G))
        if R is None:raise ValueError('discriminant stabilizer not integral')
        generators.append(R)
    context=dict(group=group,base=base,subset=subset,local=a,G=G,
                 index=int(libgap.Size(group))//int(libgap.Size(H)))
    return generators,ZZ(libgap.Size(H)),context


def integral_coset(R,representatives):
    """Correct the rational map R by the stabilizer data to an integral isometry, or None.

    有理数の写像 R を、安定化群のデータで補正して整数の等長写像にする。できなければ None。"""
    if R is None:return None
    if 'group' in representatives:
        context=representatives;target=context['subset'](R.inverse())
        if target is None:return None
        p=libgap.RepresentativeAction(context['group'],context['base'],libgap(target),libgap.OnSets)
        if p==libgap.fail:return None
        a=context['local'];n=len(a['V'])
        perm=np.array([int(x)-1 for x in libgap.ListPerm(p,n)],dtype=int)[:n]
        A=lift(a['V'],a['V'][perm],context['G'],context['G'])
        S=integral(A*R)
        if S is None:raise ValueError('glue coset failed integer check')
        return S
    A=representatives.get(lattice_key(R.inverse()))
    return integral(A*R) if A is not None else None


def local_map(a,b,G,H):
    """The rational map between two local graphs with equal canonical Gram matrices, or None.

    標準の Gram 行列が等しい二つの局所グラフの間の有理数の写像。等しくなければ None。"""
    ca,cb=a['canonical'],b['canonical']
    if not np.array_equal(a['Q'][np.ix_(ca,ca)],b['Q'][np.ix_(cb,cb)]):return None
    return lift(a['V'][ca],b['V'][cb],G,H)


def automorphisms(shell):
    """|Aut L| and generators (column convention g^T G g = G) of the lattice of `shell`.

    shell の格子の |Aut L| と生成元（列の約束 g^T G g = G）。"""
    a = shell.local(shell.anchor)
    generators,stabilizer,representatives=integral_stabilizer(a,shell.G)
    perms=[shell.permutation(R) for R in generators]
    orbit={shell.anchor}; rejected=set();attempts=0
    def close_orbit():
        todo=list(orbit)
        while todo:
            i=todo.pop()
            for p in perms:
                j=int(p[i])
                if j not in orbit:orbit.add(j);todo.append(j)
    for target in shell.candidates:
        if target in orbit:continue
        attempts+=1
        b=shell.local(target)
        R=integral_coset(local_map(a,b,shell.G,shell.G),representatives)
        if R is None:rejected.add(target);continue
        generators.append(R);perms.append(shell.permutation(R));close_orbit()
        progress(stage='orbit', orbit_lines=len(orbit), attempted=attempts,
                 candidate_lines=len(shell.candidates))
    if orbit & rejected:raise ValueError('inconsistent orbit classification')
    generators.append(-matrix.identity(ZZ,shell.G.nrows()))
    result=dict(order=str(2*len(orbit)*stabilizer),generators=[rows(R.T) for R in generators],
                method='mod2-profile/local-Gram-bliss/orbit-stabilizer',check=True,
                stabilizer_order=str(stabilizer),orbit_lines=len(orbit),candidate_lines=len(shell.candidates),
                local_graph_attempts=attempts+1,local_stats=a['stats'],lattice_orbit=representatives.get('index',len(representatives)))
    # g is returned in the column convention (g^T G g = G), whatever the row action used inside.
    # g は列の約束（g^T G g = G）で返す。内部で使う行の作用には依らない。
    for g in result['generators']:
        g=matrix(ZZ,g)
        if g.T*shell.G*g!=shell.G:raise ValueError('generator check failed')
    return result


def isometry(first, second):
    """Whether the lattices of two shells are isometric, with T^T G_first T = G_second when they are.

    二つの shell の格子が同型か。同型なら T^T G_first T = G_second となる T も返す。"""
    if first.key!=second.key:return dict(isometric=False,check=True,reason='profile mismatch')
    a=first.local(first.anchor)
    _,_,representatives=integral_stabilizer(a,first.G)
    for target in second.candidates:
        b=second.local(target)
        R=integral_coset(local_map(a,b,first.G,second.G),representatives)
        if R is not None:
            T=R.inverse().T
            if T.T*first.G*T!=second.G:raise ValueError('isometry check failed')
            return dict(isometric=True,transform=rows(T),check=True)
    return dict(isometric=False,check=True,reason='all anchor candidates exhausted')
