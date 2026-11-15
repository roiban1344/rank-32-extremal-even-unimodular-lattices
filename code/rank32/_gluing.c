/* The fast part of the normal form for the gluing of involutions; the explanation is in rank32/gluing.py.
   A matrix is the list uint32[s] of "the mask of the image of basis vector i". The nodes of the tree are
   created on the Python side; when one is missing, the return code asks for it (and the call is repeated).
   Compilation: rank32/gluing.py does it on first use (gcc -O2 -shared -fPIC).

   対合の貼り合わせの正規形の速い部分。説明は rank32/gluing.py にある。
   行列は「基底 i の像の mask」の列 uint32[s] で持つ。ツリーのノードは Python 側で作り、足りないときは
   戻り値で要求する（呼び出しをやり直して続ける）。
   コンパイル：rank32/gluing.py が初めて使うときに行う（gcc -O2 -shared -fPIC）。 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define MAXS 20
#define MAXD 6
typedef uint32_t u32; typedef uint64_t u64;

typedef struct {
    int small; u64 order; int nelem; const u32 *elems;          /* small group: nelem x s / 小さい群：nelem x s */
    const uint16_t *label[MAXD], *size[MAXD], *parent[MAXD]; const uint8_t *pgen[MAXD]; /* large group: indexed by the position in the domain / 大きい群：domain の添字で引く */
    const u32 *invs;                                             /* large group: ngens x s / 大きい群：ngens x s */
} node_t;


typedef struct {
    int s, ndom;
    const u32 *dom[MAXD]; const int32_t *idx[MAXD];
    node_t *nodes; int nnodes, cap;
    u64 *hk; int32_t *hv; u64 hcap, hn;                          /* (node, rep) -> child / (node, rep) -> 子 */
    /* the stages and the chain of F / F のステージと鎖 */
    int nst; const u32 *st[MAXD]; int stn[MAXD]; u32 *kx; int kn;
    int nlev; const u32 *lorb[64]; const u32 *lus[64]; int ln[64];
    /* work space / 作業領域 */
    u32 *fr, *fr2; size_t frcap;                                 /* leaf: psi[s], b[s], weight u64 / 葉：psi[s], b[s], 重み u64 */
    u64 calls, leaves, misses;
} tree_t;

tree_t *tree_new(int s, int ndom, const u32 **dom, const int32_t **idx) {
    tree_t *t = calloc(1, sizeof(tree_t));
    t->s = s; t->ndom = ndom;
    for (int i = 0; i < ndom; i++) { t->dom[i] = dom[i]; t->idx[i] = idx[i]; }
    t->cap = 1024; t->nodes = calloc(t->cap, sizeof(node_t));
    t->hcap = 1 << 12; t->hk = calloc(t->hcap, 8); t->hv = calloc(t->hcap, 4);
    t->frcap = 1024; t->fr = malloc(t->frcap * (2 * s + 2) * 4); t->fr2 = malloc(t->frcap * (2 * s + 2) * 4);
    return t;
}

int node_add(tree_t *t, u64 order, int small, int nelem, const u32 *elems, const u32 *invs) {
    if (t->nnodes == t->cap) { t->nodes = realloc(t->nodes, 2 * t->cap * sizeof(node_t)); memset(t->nodes + t->cap, 0, t->cap * sizeof(node_t)); t->cap *= 2; }
    node_t *n = &t->nodes[t->nnodes];
    n->order = order; n->small = small; n->nelem = nelem; n->elems = elems; n->invs = invs;
    return t->nnodes++;
}

void node_set_tables(tree_t *t, int id, int d, const uint16_t *label, const uint16_t *size, const uint16_t *parent, const uint8_t *pgen) {
    node_t *n = &t->nodes[id]; n->label[d] = label; n->size[d] = size; n->parent[d] = parent; n->pgen[d] = pgen;
}

static inline u64 hsh(u64 k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k; }

static int32_t child_get(tree_t *t, int id, u32 rep) {
    u64 k = ((u64)(id + 1) << 32) | rep, m = t->hcap - 1, i = hsh(k) & m;
    while (t->hk[i]) { if (t->hk[i] == k) return t->hv[i]; i = (i + 1) & m; }
    return -1;
}

void child_set(tree_t *t, int id, u32 rep, int child) {
    if (2 * (t->hn + 1) > t->hcap) {
        u64 oc = t->hcap; u64 *ok = t->hk; int32_t *ov = t->hv;
        t->hcap *= 2; t->hk = calloc(t->hcap, 8); t->hv = calloc(t->hcap, 4);
        for (u64 j = 0; j < oc; j++) if (ok[j]) { u64 i = hsh(ok[j]) & (t->hcap - 1); while (t->hk[i]) i = (i + 1) & (t->hcap - 1); t->hk[i] = ok[j]; t->hv[i] = ov[j]; }
        free(ok); free(ov);
    }
    u64 k = ((u64)(id + 1) << 32) | rep, m = t->hcap - 1, i = hsh(k) & m;
    while (t->hk[i]) i = (i + 1) & m;
    t->hk[i] = k; t->hv[i] = child; t->hn++;
}

void set_stages(tree_t *t, int nst, const u32 **st, const int *stn) {
    t->nst = nst; for (int i = 0; i < nst; i++) { t->st[i] = st[i]; t->stn[i] = stn[i]; }
    /* points used to compare leaves: from the first stage until there are at least 64 (a union of invariant sets, hence invariant)
       葉の比較に使う点：先頭のステージから合計 64 点以上になるまで（不変な集合の和なので不変） */
    int k = 0, tot = 0; while (k < nst && tot < 64) tot += stn[k++];
    t->kx = malloc(tot * 4); t->kn = tot; tot = 0;
    for (int i = 0; i < k; i++) { memcpy(t->kx + tot, st[i], stn[i] * 4); tot += stn[i]; }
}
void set_chain(tree_t *t, int nlev, const u32 **lorb, const u32 **lus, const int *ln) { t->nlev = nlev; for (int i = 0; i < nlev; i++) { t->lorb[i] = lorb[i]; t->lus[i] = lus[i]; t->ln[i] = ln[i]; } }

/* tables of the matrix A by 8 bits (2 tables if s <= 16, else 3) / 行列 A の 8 ビットずつの表（s <= 16 なら 2 枚。それ以上は 3 枚） */
typedef struct { u32 t[3][256]; } ltab;
static inline void mk(ltab *L, const u32 *A, int s) {
    for (int p = 0; p < 3; p++) {
        int lo = 8 * p; L->t[p][0] = 0;
        for (int v = 1; v < 256; v++) {
            int b = __builtin_ctz(v); int i = lo + b;
            L->t[p][v] = L->t[p][v & (v - 1)] ^ (i < s ? A[i] : 0);
        }
    }
}
static inline u32 ap(const ltab *L, u32 x) { return L->t[0][x & 255] ^ L->t[1][(x >> 8) & 255] ^ L->t[2][(x >> 16) & 255]; }
static inline u32 ap_slow(const u32 *A, u32 x) { u32 r = 0; while (x) { r ^= A[__builtin_ctz(x)]; x &= x - 1; } return r; }

typedef struct { u32 rep, cnt, size; } hent;
static int cmp_u32(const void *a, const void *b) { u32 x = *(const u32 *)a, y = *(const u32 *)b; return x < y ? -1 : x > y; }
static int cmp64(const void *a, const void *b) { u64 x = *(const u64 *)a, z = *(const u64 *)b; return x < z ? -1 : x > z; }
static int cmp_h(const void *a, const void *b) {
    const hent *x = a, *y = b;
    if (x->cnt != y->cnt) return x->cnt < y->cnt ? -1 : 1;
    if (x->size != y->size) return x->size < y->size ? -1 : 1;
    return x->rep < y->rep ? -1 : x->rep > y->rep;
}

/* bring psi to its least form along the chain of F / psi を F の鎖で最小の形へ */
static void chain_min(tree_t *t, u32 *psi) {
    int s = t->s; ltab L; u32 nw[MAXS];
    for (int l = 0; l < t->nlev; l++) {
        mk(&L, psi, s);
        const u32 *o = t->lorb[l]; int n = t->ln[l]; u32 mn = ~0u; int j = 0;
        for (int i = 0; i < n; i++) { u32 r = ap(&L, o[i]); if (r < mn) { mn = r; j = i; } }
        const u32 *u = t->lus[l] + (size_t)j * s;
        for (int i = 0; i < s; i++) nw[i] = ap(&L, u[i]);
        memcpy(psi, nw, s * 4);
    }
}

/* leaves with the same form along the chain of F have the same subtree, so they are merged and their weights added
   (if psi1 = psi2 a, then b2^{-1} b1 is in Stab)
   F の鎖の形が同じ葉は、その先のツリーが同じなので、1 つにまとめて重みを足す（psi1 = psi2 a のとき b2^{-1} b1 は Stab の元） */
static __thread const u32 *g_forms; static __thread int g_s;
static int cmp_form(const void *a, const void *b) {
    u32 x = *(const u32 *)a, z = *(const u32 *)b;
    int c = memcmp(g_forms + (size_t)x * g_s, g_forms + (size_t)z * g_s, g_s * 4);
    return c ? c : (x < z ? -1 : x > z);
}

static void invert(const u32 *A, u32 *R, int s) {
    u32 M[MAXS], I[MAXS];
    /* eliminate on the columns of A (images of the basis) as rows; R[i] is the image of basis vector i under the map with R(A[j]) = e_j
       A の列（基底の像）を行とみて掃き出す。R[i] は R(A[j]) = e_j となる線形写像の基底 i の像 */
    for (int i = 0; i < s; i++) { M[i] = A[i]; I[i] = 1u << i; }
    for (int c = 0; c < s; c++) {
        int p = c; while (!((M[p] >> c) & 1)) p++;
        u32 t = M[p]; M[p] = M[c]; M[c] = t; t = I[p]; I[p] = I[c]; I[c] = t;
        for (int r = 0; r < s; r++) if (r != c && ((M[r] >> c) & 1)) { M[r] ^= M[c]; I[r] ^= I[c]; }
    }
    /* now M[c] = e_c and I[c] is the combination of the basis giving e_c: R(e_c) = I[c]
       いま M[c] = e_c、I[c] は e_c を与える基底の組合せ：R(e_c) = I[c] */
    for (int c = 0; c < s; c++) R[c] = I[c];
}

static void add_elem(const u32 *b0, const u32 *b1, u32 *elems, int maxel, int *ne, int s) {
    if (!elems || *ne >= maxel) return;
    u32 R[MAXS]; invert(b0, R, s);
    u32 *o = elems + (size_t)(*ne) * s;
    for (int i = 0; i < s; i++) o[i] = ap_slow(R, b1[i]);
    (*ne)++;
}

static inline u64 getw(const u32 *it, int s) { return ((u64)it[2 * s + 1] << 32) | it[2 * s]; }
static inline void setw(u32 *it, int s, u64 v) { it[2 * s] = (u32)v; it[2 * s + 1] = (u32)(v >> 32); }

/* merge the front fr[0..nf); returns the new nf / フロンティア fr[0..nf) をまとめる。戻り値は新しい nf */
static size_t merge(tree_t *t, u32 *fr, size_t nf, u32 *elems, int maxel, int *ne) {
    if (nf < 2) return nf;
    int s = t->s; size_t w = 2 * s + 2;
    /* the buffers are sized by the number of entries nf * s, not by nf alone: s can differ from one call to the next
       in the same thread
       領域の大きさは要素数 nf * s で決める（件数 nf だけでは足りない）。同じスレッドの中でも、呼び出しごとに s が違いうる */
    static __thread u32 *forms = NULL, *ord = NULL; static __thread size_t fcap = 0, ocap = 0;
    if (fcap < nf * (size_t)s) { fcap = 2 * nf * (size_t)s; forms = realloc(forms, fcap * 4); }
    if (ocap < nf) { ocap = 2 * nf; ord = realloc(ord, ocap * 4); }
    for (size_t f = 0; f < nf; f++) {
        memcpy(forms + f * s, fr + f * w, s * 4);
        chain_min(t, forms + f * s);
        ord[f] = f;
    }
    g_forms = forms; g_s = s;
    qsort(ord, nf, 4, cmp_form);
    size_t out = 0;
    static __thread u32 *tmp = NULL; static __thread size_t tc = 0;      /* for the same reason, sized by nf * w / 同じ理由で、要素数（nf * w）で比べる */
    if (tc < nf * (size_t)w) { tc = 2 * nf * (size_t)w; tmp = realloc(tmp, tc * 4); }
    for (size_t i = 0; i < nf; ) {
        size_t j = i + 1; u32 fi = ord[i]; u64 wt = getw(fr + (size_t)fi * w, s);
        while (j < nf && !memcmp(forms + (size_t)ord[j] * s, forms + (size_t)fi * s, s * 4)) {
            const u32 *it = fr + (size_t)ord[j] * w;
            wt += getw(it, s);
            add_elem(fr + (size_t)fi * w + s, it + s, elems, maxel, ne, s);
            j++;
        }
        u32 *o = tmp + out * w; memcpy(o, fr + (size_t)fi * w, w * 4); setw(o, s, wt); out++;
        i = j;
    }
    memcpy(fr, tmp, out * w * 4);
    return out;
}

/* Return value: 0 done, 1 a child is needed (need = node, rep, d), 2 tables are needed (need = node, d).
   The normal form goes to form[s] and |Stab| to *stab. If elems != NULL, up to maxel elements b of Stab on the K side
   are written to elems (their number to nel). They arise while merging, so they need not generate all of Stab
   (they are used only to thin out neighbours; too few is safe).
   戻り値：0 完了、1 子が要る (need = node, rep, d)、2 表が要る (need = node, d)。
   form[s] に正規形、*stab に |Stab|。elems != NULL なら Stab の K 側の元 b を maxel 個まで elems に入れる（nel に個数）。
   まとめの過程で出た元なので、Stab 全体を生成するとは限らない（使い道は近傍の間引きだけで、足りなくても安全）。 */
int canon(tree_t *t, const u32 *phi, u32 *form, u64 *stab, u32 *need, u32 *elems, int maxel, int *nel) {
    int s = t->s; size_t w = 2 * s + 2; ltab L;
    static __thread u32 *y = NULL, *tmpy = NULL, *keys = NULL, *bestk = NULL; static __thread size_t ycap = 0;
    size_t maxst = t->kn; for (int i = 0; i < t->nst; i++) if ((size_t)t->stn[i] > maxst) maxst = t->stn[i];
    if (ycap < maxst + 1) { ycap = maxst + 1; y = realloc(y, ycap * 4); tmpy = realloc(tmpy, ycap * 4); keys = realloc(keys, ycap * 4); bestk = realloc(bestk, ycap * 4); }
    static __thread hent *h = NULL, *hb = NULL; static __thread size_t hcap = 0;
    if (hcap < maxst + 1) { hcap = maxst + 1; h = realloc(h, hcap * sizeof(hent)); hb = realloc(hb, hcap * sizeof(hent)); }
    size_t nf = 1; int ne = 0;
    u32 *fr = t->fr, *fr2 = t->fr2;
    for (int i = 0; i < s; i++) { fr[i] = phi[i]; fr[s + i] = 1u << i; }
    setw(fr, s, 1);
    int id = 0; t->calls++;
    static __thread u32 *sel = NULL; static __thread size_t selcap = 0;
    while (id >= 0 && t->nodes[id].order > 1) {
        const node_t *n = &t->nodes[id];
        if (n->small) {
            /* small group: try every element; keep only the (leaf, element) with the least function (minimum and sum of squares) of the set Y
               小さい群：元を全部試す。Y の集合の関数（最小値と平方和）が最小の (葉, 元) だけ残す */
            const u32 *st0 = t->kx; int n0 = t->kn; u64 bk = ~0ULL; size_t nf2 = 0; ltab Lh;
            for (size_t f = 0; f < nf; f++) {
                const u32 *psi = fr + f * w;
                mk(&L, psi, s);
                for (int i = 0; i < n0; i++) y[i] = ap(&L, st0[i]);
                for (int e = 0; e < n->nelem; e++) {
                    const u32 *hh = n->elems + (size_t)e * s;
                    mk(&Lh, hh, s);
                    u32 mn = ~0u; u32 sum = 0;
                    for (int i = 0; i < n0; i++) { u32 v = ap(&Lh, y[i]); if (v < mn) mn = v; sum += v * v; }
                    u64 k = ((u64)mn << 32) | sum;
                    if (k > bk) continue;
                    if (k < bk) { bk = k; nf2 = 0; }
                    if (nf2 == t->frcap) {
                        t->frcap *= 2; t->fr = realloc(t->fr, t->frcap * w * 4); t->fr2 = realloc(t->fr2, t->frcap * w * 4);
                        fr = t->fr; fr2 = t->fr2; psi = fr + f * w;
                    }
                    u32 *o = fr2 + nf2 * w;
                    for (int a = 0; a < 2 * s; a++) o[a] = ap(&Lh, psi[a]);
                    o[2 * s] = psi[2 * s]; o[2 * s + 1] = psi[2 * s + 1];
                    nf2++;
                }
            }
            u32 *tp = fr; fr = fr2; fr2 = tp; t->fr = fr; t->fr2 = fr2; nf = nf2;
            t->leaves += nf;
            id = -1; break;
        }
        int bestd = -1; size_t besth = 0; u32 bestrep = 0;
        size_t nsel = 0, used = 0;
        if (selcap < 16) { selcap = 1024; sel = realloc(sel, selcap * 4); }
        for (size_t f = 0; f < nf; f++) {
            u32 *psi = fr + f * w;
            mk(&L, psi, s);
            int d = -1; size_t m = 0; size_t hn = 0;
            for (int dd = 0; dd < t->nst; dd++) {
                int ns = t->stn[dd];
                const u32 *st = t->st[dd];
                if (!n->label[dd]) { need[0] = id; need[1] = dd; return 2; }
                m = 0;
                for (int i = 0; i < ns; i++) {
                    u32 yy = ap(&L, st[i]);
                    int32_t ii = t->idx[dd][yy]; int32_t l = n->label[dd][ii]; u32 size = n->size[dd][l];
                    if (size > 1) { keys[m] = t->dom[dd][l]; tmpy[m] = size; y[m] = yy; m++; }
                }
                if (m) { d = dd; break; }
            }
            if (d < 0) { need[0] = id; return -1; }   /* the group of K acts trivially on the points but has order > 1 / K の群が点の上で自明なのに位数 > 1 */
            {
                static __thread u64 *o64 = NULL; static __thread size_t oc = 0;
                if (oc < m + 1) { oc = 2 * (m + 1); o64 = realloc(o64, oc * 8); }
                for (size_t i = 0; i < m; i++) o64[i] = ((u64)keys[i] << 32) | tmpy[i];
                qsort(o64, m, 8, cmp64);
                for (size_t i = 0; i < m; ) {
                    size_t j = i; while (j < m && (o64[j] >> 32) == (o64[i] >> 32)) j++;
                    h[hn].rep = o64[i] >> 32; h[hn].size = (u32)o64[i]; h[hn].cnt = j - i; hn++; i = j;
                }
            }
            qsort(h, hn, sizeof(hent), cmp_h);
            int c;
            if (bestd >= 0) {
                if (d != bestd) c = d < bestd ? -1 : 1;
                else {
                    c = 0;
                    if (hn != besth) c = hn < besth ? -1 : 1;
                    for (size_t i = 0; !c && i < hn; i++) c = cmp_h(&h[i], &hb[i]);
                }
            } else c = -1;
            if (c < 0) { bestd = d; besth = hn; memcpy(hb, h, hn * sizeof(hent)); bestrep = h[0].rep; nsel = 0; used = 0; }
            if (c <= 0) {
                u32 need_sz = 2 + h[0].cnt;
                if (used + need_sz > selcap) { selcap = 2 * (used + need_sz); sel = realloc(sel, selcap * 4); }
                sel[used] = f; sel[used + 1] = h[0].cnt; size_t q = used + 2;
                for (size_t i = 0; i < m; i++) if (keys[i] == h[0].rep) sel[q++] = y[i];
                used = q; nsel++;
            }
        }
        int32_t ch = -1;
        {
            u64 co = n->order / (u64)n->size[bestd][n->label[bestd][t->idx[bestd][bestrep]]];
            if (co > 1) {
                ch = child_get(t, id, bestrep);
                if (ch < 0) { need[0] = id; need[1] = bestrep; need[2] = bestd; return 1; }
            }
        }
        size_t tot = 0, q = 0;
        for (size_t i = 0; i < nsel; i++) { tot += sel[q + 1]; q += 2 + sel[q + 1]; }
        if (tot > t->frcap) {
            t->frcap = 2 * tot; t->fr = realloc(t->fr, t->frcap * w * 4); t->fr2 = realloc(t->fr2, t->frcap * w * 4); fr = t->fr; fr2 = t->fr2;
        }
        size_t nf2 = 0; q = 0;
        u32 B[MAXS];
        const uint16_t *par = n->parent[bestd]; const uint8_t *pg = n->pgen[bestd];
        for (size_t i = 0; i < nsel; i++) {
            u32 fi = sel[q], cnt = sel[q + 1]; const u32 *ys = sel + q + 2; q += 2 + cnt;
            const u32 *psi = fr + (size_t)fi * w;
            for (u32 j = 0; j < cnt; j++) {
                /* walk up the tree: point c satisfies tabs[pgen[c]][parent[c]] = c, so multiply by the inverses in turn
                   ツリーを親へ辿る：点 c は tabs[pgen[c]][parent[c]] = c なので、逆元を順に掛ける */
                for (int a = 0; a < s; a++) B[a] = 1u << a;
                int32_t c = t->idx[bestd][ys[j]];
                while (par[c] != c) {
                    const u32 *inv = n->invs + (size_t)pg[c] * s;
                    for (int a = 0; a < s; a++) B[a] = ap_slow(inv, B[a]);
                    c = par[c];
                }
                mk(&L, B, s);
                u32 *o = fr2 + nf2 * w;
                for (int a = 0; a < 2 * s; a++) o[a] = ap(&L, psi[a]);
                o[2 * s] = psi[2 * s]; o[2 * s + 1] = psi[2 * s + 1];
                nf2++;
            }
        }
        u32 *tp = fr; fr = fr2; fr2 = tp; t->fr = fr; t->fr2 = fr2; nf = nf2; id = ch;
        nf = merge(t, fr, nf, elems, maxel, &ne);
        t->leaves += nf;
    }
    /* leaves: keep those with the least set Y, then the least form along the chain of F / 葉：Y の集合の最小で絞り、F の鎖で最小の形へ */
    const u32 *st0 = t->kx; int n0 = t->kn;
    int have = 0; u32 best[MAXS]; u64 cnt = 0; size_t b0 = 0;
    size_t kept = 0;
    for (size_t f = 0; f < nf; f++) {
        u32 *psi = fr + f * w;
        mk(&L, psi, s);
        for (int i = 0; i < n0; i++) keys[i] = ap(&L, st0[i]);
        qsort(keys, n0, 4, cmp_u32);
        int c = kept ? memcmp(keys, bestk, n0 * 4) : -1;
        if (c < 0) { memcpy(bestk, keys, n0 * 4); kept = 1; memmove(fr, psi, w * 4); }
        else if (c == 0) { memmove(fr + kept * w, psi, w * 4); kept++; }
    }
    for (size_t f = 0; f < kept; f++) {
        u32 *psi = fr + f * w; u32 tmp[MAXS]; memcpy(tmp, psi, s * 4);
        chain_min(t, tmp);
        int c = have ? memcmp(tmp, best, s * 4) : -1;
        if (c < 0) { memcpy(best, tmp, s * 4); have = 1; cnt = getw(psi, s); b0 = f; }
        else if (c == 0) { cnt += getw(psi, s); add_elem(fr + b0 * w + s, psi + s, elems, maxel, &ne, s); }
    }
    memcpy(form, best, s * 4); *stab = cnt; if (nel) *nel = ne;
    return 0;
}

/* helpers that speed up building the tree on the Python side; a matrix is int64[s] (mask of the image of basis vector i)
   Python 側のツリーの作成を速くするヘルパー。行列は int64[s]（基底 i の像の mask） */
void mat_apply64(const int64_t *A, int s, const int64_t *m, size_t n, int64_t *out) {
    if (n <= 16) {
        for (size_t k = 0; k < n; k++) { u64 x = (u64)m[k]; int64_t r = 0; while (x) { r ^= A[__builtin_ctzll(x)]; x &= x - 1; } out[k] = r; }
        return;
    }
    static __thread int64_t T[3][256];
    for (int p = 0; p < 3; p++) {
        int lo = 8 * p; T[p][0] = 0;
        for (int v = 1; v < 256; v++) { int i = lo + __builtin_ctz(v); T[p][v] = T[p][v & (v - 1)] ^ (i < s ? A[i] : 0); }
    }
    for (size_t k = 0; k < n; k++) { u64 x = (u64)m[k]; out[k] = T[0][x & 255] ^ T[1][(x >> 8) & 255] ^ T[2][(x >> 16) & 255]; }
}

/* inverse matrix; -1 if singular / 逆行列。特異なら -1 */
int mat_invert64(const int64_t *A, int64_t *R, int s) {
    u32 M[MAXS], I[MAXS];
    for (int i = 0; i < s; i++) { M[i] = (u32)A[i]; I[i] = 1u << i; }
    for (int c = 0; c < s; c++) {
        int p = c; while (p < s && !((M[p] >> c) & 1)) p++;
        if (p == s) return -1;
        u32 t = M[p]; M[p] = M[c]; M[c] = t; t = I[p]; I[p] = I[c]; I[c] = t;
        for (int r = 0; r < s; r++) if (r != c && ((M[r] >> c) & 1)) { M[r] ^= M[c]; I[r] ^= I[c]; }
    }
    for (int c = 0; c < s; c++) R[c] = I[c];
    return 0;
}

/* orbits of the permutations tabs: label is the least index of the orbit, size[label] its size, and the tree tabs[pgen[c]][parent[c]] = c
   (the representative has parent[c] = c)
   置換の組 tabs の軌道。label は軌道の最小の添字、size[label] は軌道の大きさ、tabs[pgen[c]][parent[c]] = c のツリー（代表は parent[c] = c） */
void orbit_tables(int n, int ngens, const int32_t **tabs, uint16_t *label, uint16_t *size, uint16_t *parent, uint8_t *pgen) {
    int32_t *q = malloc((size_t)n * 4); uint8_t *seen = calloc(n, 1);
    memset(size, 0, (size_t)n * 2); memset(pgen, 0, n);
    for (int i = 0; i < n; i++) {
        if (seen[i]) continue;
        int h = 0, e = 0; q[e++] = i; seen[i] = 1; parent[i] = i; label[i] = i;
        while (h < e) {
            int c = q[h++];
            for (int g = 0; g < ngens; g++) {
                int im = tabs[g][c];
                if (!seen[im]) { seen[im] = 1; parent[im] = c; pgen[im] = g; label[im] = i; q[e++] = im; }
            }
        }
        size[i] = (uint16_t)e;
    }
    free(q); free(seen);
}

/* list all elements; gens is ngens x s, out is maxn x s; returns the number of elements (-1 if it exceeds maxn)
   元を全部並べる。gens は ngens x s、out は maxn x s。戻り値は元の数（maxn を超えたら -1） */
long mat_enumerate(const int64_t *gens, int ngens, int s, long maxn, int64_t *out) {
    u64 cap = 4; while (cap < 2 * (u64)maxn + 2) cap *= 2;
    int32_t *h = calloc(cap, 4);
    long cnt = 1;
    for (int i = 0; i < s; i++) out[i] = (int64_t)1 << i;
    {
        u64 v = 14695981039346656037ULL; for (int i = 0; i < s; i++) { v ^= (u64)out[i]; v *= 1099511628211ULL; }
        u64 k = hsh(v) & (cap - 1); h[k] = 1;
    }
    int64_t z[MAXS];
    for (long e = 0; e < cnt; e++) {
        for (int g = 0; g < ngens; g++) {
            const int64_t *G = gens + (size_t)g * s; const int64_t *E = out + (size_t)e * s;
            u64 v = 14695981039346656037ULL;
            for (int i = 0; i < s; i++) {
                u64 x = (u64)E[i]; int64_t r = 0; while (x) { r ^= G[__builtin_ctzll(x)]; x &= x - 1; }
                z[i] = r; v ^= (u64)r; v *= 1099511628211ULL;
            }
            u64 k = hsh(v) & (cap - 1); int found = 0;
            while (h[k]) { if (!memcmp(out + (size_t)(h[k] - 1) * s, z, s * 8)) { found = 1; break; } k = (k + 1) & (cap - 1); }
            if (found) continue;
            if (cnt == maxn) { free(h); return -1; }
            memcpy(out + (size_t)cnt * s, z, s * 8); h[k] = (int32_t)(cnt + 1); cnt++;
        }
    }
    free(h);
    return cnt;
}

/* Randomised Schreier-Sims to confirm that the group H of the generators has order >= order.
   Every point of every basic orbit is reached by an element of H, so the product of the orbit sizes is a lower bound
   for |H|; when it reaches order, |H| >= order holds regardless of chance. If the caller knows that H lies in a group
   of order order, H equals that group. Return 1 = confirmed (the product equals order), 0 = not confirmed (the product
   stays below order: too few generators, or bad luck; or it exceeds order, which contradicts the caller's group).
   乱択 Schreier-Sims で「生成元の群 H の位数が order 以上」を確かめる。
   各ステージの軌道の元はすべて H の元で実現されているので、軌道の大きさの積は |H| の下界であり、積が order に達すれば
   確率によらず |H| >= order。呼び出し側で H が位数 order の群に含まれていれば、H はその群に等しい。
   戻り値 1 = 確かめた（積が order に等しい）、0 = 確かめられなかった（積が order に届かない：生成元が足りないか、乱択の運。
   または積が order を超えた：呼び出し側の群と矛盾する）。 */
typedef struct { u32 base; int norb, cap; u32 *pts; u32 *tr; int ng, gcap; u32 *gens; int32_t *pos; } rss_t;
static inline void comp(const u32 *A, const u32 *B, u32 *o, int s) { for (int i = 0; i < s; i++) o[i] = ap_slow(A, B[i]); }
static void rss_addpt(rss_t *l, u32 pt, const u32 *u, int s) {
    if (l->norb == l->cap) { l->cap *= 2; l->pts = realloc(l->pts, (size_t)l->cap * 4); l->tr = realloc(l->tr, (size_t)l->cap * s * 4); }
    l->pts[l->norb] = pt; memcpy(l->tr + (size_t)l->norb * s, u, s * 4); l->pos[pt] = l->norb; l->norb++;
}
static void rss_addgen(rss_t *l, const u32 *g, int s) {
    if (l->ng == l->gcap) { l->gcap *= 2; l->gens = realloc(l->gens, (size_t)l->gcap * s * 4); }
    memcpy(l->gens + (size_t)l->ng * s, g, s * 4); int newg = l->ng++; u32 z[MAXS];
    /* apply the new generator to all current points, and all generators to new points / 新しい生成元を今の全部の点に、新しい点には全部の生成元を */
    int done_old = l->norb;
    for (int p = 0; p < l->norb; p++) {
        int g0 = p < done_old ? newg : 0;
        for (int j = g0; j < l->ng; j++) {
            const u32 *G = l->gens + (size_t)j * s; u32 y = ap_slow(G, l->pts[p]);
            if (l->pos[y] < 0) { comp(G, l->tr + (size_t)p * s, z, s); rss_addpt(l, y, z, s); }
        }
    }
}
int verify_order(const int64_t *gens, int ngens, int s, u64 order, u64 seed, int maxfail) {
    if (order == 1) return 1;
    u64 rs = seed * 2 + 1;
    #define RND() (rs ^= rs << 13, rs ^= rs >> 7, rs ^= rs << 17, rs)
    int np = ngens < 10 ? 10 : ngens; u32 *pool = malloc((size_t)np * s * 4);
    for (int i = 0; i < np; i++) for (int a = 0; a < s; a++) pool[(size_t)i * s + a] = (u32)gens[(size_t)(i % ngens) * s + a];
    rss_t L[MAXS + 1]; int nl = 0; u64 prod = 1; int fail = 0, ok = 0; u32 g[MAXS], R[MAXS], z[MAXS];
    size_t npt = (size_t)1 << s;
    for (int it = 0; fail < maxfail; it++) {
        /* product replacement; the first rounds only mix / 積の置き換え。最初の数回は混ぜるだけ */
        for (int rep = 0; rep < (it == 0 ? 30 : 1); rep++) {
            int i = RND() % np, j = RND() % (np - 1); j += j >= i;
            u32 *Pi = pool + (size_t)i * s; const u32 *Pj = pool + (size_t)j * s;
            if (RND() & 1) comp(Pi, Pj, z, s); else comp(Pj, Pi, z, s);
            memcpy(Pi, z, s * 4); memcpy(g, z, s * 4);
        }
        /* sift / ふるう */
        int l = 0, added = 0;
        for (; l < nl; l++) {
            u32 y = ap_slow(g, L[l].base);
            if (L[l].pos[y] < 0) {
                u64 before = L[l].norb; rss_addgen(&L[l], g, s);
                prod = prod / before * (u64)L[l].norb; added = 1; break;
            }
            invert((const u32 *)(L[l].tr + (size_t)L[l].pos[y] * s), R, s);
            comp(R, g, z, s); memcpy(g, z, s * 4);
        }
        if (!added) {
            int mv = -1; for (int i = 0; i < s; i++) if (g[i] != (1u << i)) { mv = i; break; }
            if (mv >= 0) {
                rss_t *n = &L[nl++]; n->base = 1u << mv; n->norb = 0; n->cap = 16; n->pts = malloc(16 * 4); n->tr = malloc((size_t)16 * s * 4);
                n->ng = 0; n->gcap = 8; n->gens = malloc((size_t)8 * s * 4); n->pos = malloc(npt * 4); memset(n->pos, 0xff, npt * 4);
                u32 id[MAXS]; for (int i = 0; i < s; i++) id[i] = 1u << i;
                rss_addpt(n, n->base, id, s); rss_addgen(n, g, s);
                prod *= (u64)n->norb; added = 1;
            }
        }
        if (added) { fail = 0; if (prod >= order) { ok = prod == order; break; } } else fail++;
    }
    for (int i = 0; i < nl; i++) { free(L[i].pts); free(L[i].tr); free(L[i].gens); free(L[i].pos); }
    free(pool);
    return ok;
}

void stats(tree_t *t, u64 *o) { o[0] = t->calls; o[1] = t->leaves; o[2] = t->misses; o[3] = t->nnodes; }
