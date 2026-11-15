/* Automorphisms and canonical form of a directed graph with coloured edges (bliss); the core of rank32/canonical.py.
   Input (binary, int32 values): n, m, then m edges (i, j, c). c is a non-zero colour number (the edge points from i to j).
   Directed coloured edges are encoded in layers: for each colour c, every vertex v gets copies O_c(v), I_c(v)
   joined to v, and the edge (i, j, c) becomes O_c(i) - I_c(j). The automorphisms found by bliss restrict
   faithfully to permutations of the vertices 0..n-1.
   Output: <out>.gens (the generators as permutations, int32 x n per generator), <out>.lab (the canonical
   labelling, int32 x n), and JSON on standard output.
   hash is a hash of the set of edges in canonical labelling that does not depend on their order (equal for
   isometric lattices; when two hashes agree the caller builds the map and checks it).

   色つき有向グラフの自己同型と標準形（bliss）。rank32/canonical.py の本体。
   入力（バイナリ、int32 の並び）：n, m, 続けて辺 m 本の (i, j, c)。c は 0 でない色の番号（i -> j の向きつき）。
   有向・色つきの辺は層で表す：色 c ごとに各頂点 v の複製 O_c(v), I_c(v) を作り、v - O_c(v), v - I_c(v) を
   結び、辺 (i, j, c) は O_c(i) - I_c(j) とする。bliss の自己同型は頂点 0..n-1 の置換に忠実に落ちる。
   出力：<out>.gens（生成元の置換、生成元ごとに int32 x n）、<out>.lab（標準形のラベル、int32 x n）、標準出力に JSON。
   hash は、標準形での辺の集合の、順序によらないハッシュ（同型な格子なら一致する。一致したときの判定は、
   呼び出し側で写像を検算する）。 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bliss/bliss_C.h>

static FILE *genout;
static uint32_t core;
static unsigned long ngens;

static void hook(void *p, unsigned n, const unsigned *a) {
  (void)p; (void)n;
  if (fwrite(a, 4, core, genout) != core) exit(4);
  ngens++;
}

static uint64_t mix(uint64_t x) {
  x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27; x *= 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static int cmp_int(const void *a, const void *b) {
  int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
  return (x > y) - (x < y);
}

int main(int argc, char **argv) {
  if (argc != 3) return 2;
  FILE *f = fopen(argv[1], "rb");
  if (!f) return 3;
  int32_t head[2];
  if (fread(head, 4, 2, f) != 2) return 3;
  uint32_t n = (uint32_t)head[0];
  uint64_t m = (uint64_t)(uint32_t)head[1];
  int32_t *e = malloc(12 * m);
  if (!e || fread(e, 12, m, f) != m) return 3;
  fclose(f);
  core = n;
  /* the list of colours, ascending; the layer number follows the colour value, so it does not depend on the lattice
     色の一覧（昇順）。層の番号は色の値の順なので、格子によらない */
  int32_t *cs = malloc(4 * (m + 1));
  for (uint64_t k = 0; k < m; k++) cs[k] = e[3 * k + 2];
  qsort(cs, m, 4, cmp_int);
  uint64_t nc = 0;
  for (uint64_t k = 0; k < m; k++) if (k == 0 || cs[k] != cs[k - 1]) cs[nc++] = cs[k];
  BlissGraph *g = bliss_new(0);
  for (uint32_t i = 0; i < n; i++) bliss_add_vertex(g, 0);
  uint32_t *base = malloc(4 * (nc + 1));
  uint32_t count = n;
  for (uint64_t c = 0; c < nc; c++) {
    base[c] = count;
    /* vertex colours come from the colour value: layers of the same colour get the same vertex colour, O and I differ
       頂点の色は色の値から決める：同じ色の層は同じ頂点色。O と I で別の色 */
    for (uint32_t i = 0; i < n; i++) {
      unsigned vo = bliss_add_vertex(g, (unsigned)(2 * c + 1));
      unsigned vi = bliss_add_vertex(g, (unsigned)(2 * c + 2));
      bliss_add_edge(g, i, vo);
      bliss_add_edge(g, i, vi);
    }
    count += 2 * n;
  }
  for (uint64_t k = 0; k < m; k++) {
    int32_t key = e[3 * k + 2];
    int32_t *pos = bsearch(&key, cs, nc, 4, cmp_int);
    uint64_t c = (uint64_t)(pos - cs);
    bliss_add_edge(g, base[c] + 2 * (uint32_t)e[3 * k], base[c] + 2 * (uint32_t)e[3 * k + 1] + 1);
  }
  char path[4096];
  snprintf(path, sizeof path, "%s.gens", argv[2]);
  genout = fopen(path, "wb");
  BlissStats stats;
  const unsigned *lab = bliss_find_canonical_labeling(g, hook, NULL, &stats);
  fclose(genout);
  snprintf(path, sizeof path, "%s.lab", argv[2]);
  FILE *o = fopen(path, "wb");
  if (fwrite(lab, 4, n, o) != n) return 4;
  fclose(o);
  /* hash of the canonical form: the sum over edges of a mix of (lab(i), lab(j), c), colour value included
     標準形のハッシュ：色の値を含めて (lab(i), lab(j), c) の和 */
  uint64_t h1 = 0, h2 = 0;
  for (uint64_t k = 0; k < m; k++) {
    uint64_t a = lab[e[3 * k]], b = lab[e[3 * k + 1]], c = (uint64_t)(uint32_t)e[3 * k + 2];
    uint64_t x = mix(a * 0x100000001b3ULL + b * 0x9e3779b97f4a7c15ULL + c * 0xc2b2ae3d27d4eb4fULL + 1);
    h1 += x;
    h2 += mix(x ^ 0x2545f4914f6cdd1dULL);
  }
  printf("{\"group_size\":\"%.0Lf\",\"generators\":%lu,\"nodes\":%lu,\"vertices\":%u,\"hash\":\"%016llx%016llx\"}\n",
         (long double)stats.group_size_approx, ngens, stats.nof_nodes, count,
         (unsigned long long)h1, (unsigned long long)h2);
  bliss_release(g);
  free(e); free(cs); free(base);
  return 0;
}
