/* Canonical form of a coloured Gram graph (bliss); the core of rank32/involution_8_24_classes.py.
   The layers of the graph are built as in rank32/_bliss_automorphisms.c (the graph of the local Gram
   matrix), but no vertex is fixed and the vertices get the colours col.
   canon(n, Gram, col, lab, gens, maxgens) returns the number of generators (-1 if there are more than
   maxgens); lab and gens refer to the first n vertices.

   彩色 Gram グラフの標準形（bliss）。rank32/involution_8_24_classes.py の本体。
   グラフの層の作り方は rank32/_bliss_automorphisms.c（局所 Gram 行列のグラフ）と同じだが、頂点を
   固定せず、頂点に色 col を与える。
   canon(n, Gram, col, lab, gens, maxgens) は生成元の個数を返す（maxgens を超えたら -1）。lab と gens は
   先頭の n 頂点についてのもの。 */
#include <stdint.h>
#include <stdlib.h>
#include <bliss/bliss_C.h>
static uint32_t core, *out_gens, count_gens, max_gens;
static void hook(void *p, unsigned n, const unsigned *a) {
  (void)p; (void)n;
  if (count_gens < max_gens)
    for (uint32_t i = 0; i < core; i++) out_gens[(size_t)count_gens * core + i] = a[i];
  count_gens++;
}
int canon(uint32_t n, const int8_t *q, const uint32_t *col, uint32_t *lab, uint32_t *gens, uint32_t maxgens) {
  int present[256] = {0}, layer[256];
  for (uint32_t i = 0; i < n; i++)
    for (uint32_t j = i + 1; j < n; j++)
      if (q[(size_t)i * n + j]) present[q[(size_t)i * n + j] + 128] = 1;
  BlissGraph *g = bliss_new(0);
  for (uint32_t i = 0; i < n; i++) bliss_add_vertex(g, col[i]);
  uint32_t cnt = n;
  for (int v = 0; v < 256; v++)
    if (present[v]) {
      layer[v] = cnt;
      for (uint32_t i = 0; i < n; i++) { bliss_add_vertex(g, v + 16); bliss_add_edge(g, i, cnt++); }
    }
  for (uint32_t i = 0; i < n; i++)
    for (uint32_t j = i + 1; j < n; j++) {
      int v = q[(size_t)i * n + j];
      if (v) bliss_add_edge(g, layer[v + 128] + i, layer[v + 128] + j);
    }
  core = n; out_gens = gens; count_gens = 0; max_gens = maxgens;
  BlissStats st;
  const unsigned *l = bliss_find_canonical_labeling(g, hook, NULL, &st);
  for (uint32_t i = 0; i < n; i++) lab[i] = l[i];
  bliss_release(g);
  return count_gens > maxgens ? -1 : (int)count_gens;
}
