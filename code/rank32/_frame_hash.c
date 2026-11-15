/* The frame distribution a_k of an extremal even unimodular lattice of rank 32, counted exactly by XOR on L/2L.
   Called by rank32/frame_fingerprint.py.
   Input (binary): int32 n, uint32 b[n] (the classes mod 2 of the n lines of minimal vectors, pairwise distinct),
                   uint8 qbyte[4][256], uint8 bbyte[6][256][256] (tables of the quadratic form q on classes:
                   q(c) = sum_a qbyte[a][c_a] + sum_{a<b} bbyte[ab][c_a][c_b] mod 2)
   Usage: frame_hash <in.bin> <top> <threads>
   For every pair of lines {i, j} (i < j) the class c = b_i ^ b_j is counted (all n(n-1)/2 pairs). If (x_i, x_j) is even
   then q(c) = 0, otherwise 1. A class with q = 0 is the class of a minimal vector (one of the b) when the pair has inner
   product +-2, and a class containing vectors of norm 8 when the inner product is 0.
   The lines are split into groups by the top PL bits; round p counts only the pairs of groups g and g^p, so each pair
   is touched once.
   Output (one line of JSON): hist0[m] = number of classes with q = 0, not a class of a minimal vector, hit by m pairs;
   hist1 (q = 1); mincls[m] (the m of the classes of minimal vectors); top (for each group, up to top classes with q = 0,
   not minimal, with the largest m, as (class, m)); times.

   階数 32 の極値偶ユニモジュラー格子の「枠の分布」a_k を、L/2L の XOR で厳密に数える。
   rank32/frame_fingerprint.py から呼ぶ。
   入力（バイナリ）：int32 n、uint32 b[n]（最短ベクトルの直線 n 本の mod 2 の類。全部相異なる）、
                    uint8 qbyte[4][256]、uint8 bbyte[6][256][256]（類の 2 次形式 q の表。
                    q(c) = sum_a qbyte[a][c_a] + sum_{a<b} bbyte[ab][c_a][c_b] mod 2）
   使い方：frame_hash <in.bin> <top> <threads>
   直線の対 {i, j}（i < j）ごとに類 c = b_i ^ b_j を数える（全 n(n-1)/2 対）。(x_i, x_j) が偶数なら q(c) = 0、奇数なら 1。
   q = 0 の類は、内積 +-2 の対なら最短ベクトルの類（b のどれか）、内積 0 の対ならノルム 8 を含む類である。
   最上位 PL ビットで直線を組に分け、組 g と g^p の対だけを第 p 回で数えるので、各対は 1 回だけ触る。
   出力（JSON 1 行）：hist0[m] = 「q = 0 で最短ベクトルの類でなく、m 対が入った類」の個数、hist1（q = 1）、
   mincls[m]（最短ベクトルの類の m）、top（組ごとに、q = 0 の非最短の類で m の大きい順に top 個の (類, m)）、時間。 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#ifndef PL
#define PL 4      /* split by the top PL bits; the table has 2^(32-PL) uint16 (512 MB per thread for PL=4) / 最上位 PL ビットで分ける */
#endif
#define S (32 - PL)
#define TSIZE (1u << S)
#define MAXM 65535

static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static int n; static uint32_t *b; static uint8_t qb[4][256]; static uint8_t bb[6][256][256];
static inline int qform(uint32_t c){
  uint8_t c0=c, c1=c>>8, c2=c>>16, c3=c>>24;
  return (qb[0][c0]^qb[1][c1]^qb[2][c2]^qb[3][c3]^bb[0][c0][c1]^bb[1][c0][c2]^bb[2][c0][c3]^bb[3][c1][c2]^bb[4][c1][c3]^bb[5][c2][c3])&1;
}
static int cmpu(const void*x,const void*y){uint32_t a=*(const uint32_t*)x,c=*(const uint32_t*)y;return (a>c)-(a<c);}
typedef struct {uint32_t c; uint32_t m;} cm;

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: frame-hash32 in.bin top threads\n");return 2;}
  int top=atoi(argv[2]), nth=atoi(argv[3]);
  FILE*f=fopen(argv[1],"rb"); if(!f){perror("open");return 1;}
  if(fread(&n,4,1,f)!=1) return 1;
  b=malloc(4L*n); if(fread(b,4,n,f)!=(size_t)n) return 1;
  if(fread(qb,1,sizeof qb,f)!=sizeof qb) return 1;
  if(fread(bb,1,sizeof bb,f)!=sizeof bb) return 1;
  fclose(f);
  double t0=now();
  /* sort the lines by group / 組ごとに並べる */
  int G=1<<PL; int *cnt=calloc(G+1,sizeof(int));
  for(int i=0;i<n;i++) cnt[(b[i]>>S)+1]++;
  for(int g=0;g<G;g++) cnt[g+1]+=cnt[g];
  uint32_t *sb=malloc(4L*n); int *fill=malloc(sizeof(int)*G); memcpy(fill,cnt,sizeof(int)*G);
  for(int i=0;i<n;i++) sb[fill[b[i]>>S]++]=b[i];
  uint32_t *bsorted=malloc(4L*n); memcpy(bsorted,b,4L*n); qsort(bsorted,n,4,cmpu);
  long *hist0=calloc(MAXM+1,sizeof(long)), *hist1=calloc(MAXM+1,sizeof(long)), *minh=calloc(MAXM+1,sizeof(long));
  cm *tops=calloc((size_t)top*G,sizeof(cm)); int *ntops=calloc(G,sizeof(int));
  long saturated=0, min_q1=0;
  #pragma omp parallel num_threads(nth) reduction(+:saturated,min_q1)
  {
    uint16_t *T=malloc(2L*TSIZE);
    long *h0=calloc(MAXM+1,sizeof(long)), *h1=calloc(MAXM+1,sizeof(long)), *hm=calloc(MAXM+1,sizeof(long));
    #pragma omp for schedule(dynamic,1)
    for(int p=0;p<G;p++){
      memset(T,0,2L*TSIZE);
      for(int g=0;g<G;g++){int h=g^p; if(h<g) continue;
        for(int i=cnt[g];i<cnt[g+1];i++){uint32_t bi=sb[i]; int j0=(h==g)?i+1:cnt[h];
          for(int j=j0;j<cnt[h+1];j++){uint32_t c=(bi^sb[j])&(TSIZE-1); if(T[c]<MAXM) T[c]++; else saturated++;}}}
      /* take out the classes of minimal vectors / 最短ベクトルの類を抜き出す */
      for(int i=cnt[p];i<cnt[p+1];i++){uint32_t c=sb[i]&(TSIZE-1); hm[T[c]]++; if(qform(sb[i])) min_q1++; T[c]=0;}
      /* (i) q(c) is computed only when m = 36: a class with q = 1 contains a vector +-y of norm 6 and has m = 36
         (the number of decompositions y = a + b, forced by the design); with q = 0, m = 36 happens only for the
         classes of frames with k = 9. (ii) the top list is updated only when m exceeds its current minimum.
         (i) q(c) は m = 36 のときだけ計算する。q = 1 の類はノルム 6 のベクトル +-y を含み、m = 36（デザインで強制）。
         q = 0 で m = 36 になるのは k = 9 の枠の類だけ。(ii) 上位 top 個の記録は、今の最小より大きい m のときだけ入れ替える。 */
      cm *tp=tops+(size_t)p*top; int nt=0; int wmin=0;
      for(uint32_t c=0;c<TSIZE;c++){uint16_t m=T[c]; if(!m) continue; uint32_t full=((uint32_t)p<<S)|c;
        if(m==36 && qform(full)){h1[m]++; continue;}
        h0[m]++;
        if(nt<top){tp[nt].c=full;tp[nt].m=m;nt++; if(nt==top){wmin=0; for(int k=1;k<top;k++) if(tp[k].m<tp[wmin].m) wmin=k;}}
        else if(m>tp[wmin].m){tp[wmin].c=full;tp[wmin].m=m; wmin=0; for(int k=1;k<top;k++) if(tp[k].m<tp[wmin].m) wmin=k;}}
      ntops[p]=nt;
    }
    #pragma omp critical
    {for(int m=0;m<=MAXM;m++){hist0[m]+=h0[m];hist1[m]+=h1[m];minh[m]+=hm[m];}}
    free(T);free(h0);free(h1);free(hm);
  }
  double t1=now();
  printf("{\"n\":%d,\"seconds\":%.2f,\"threads\":%d,\"saturated\":%ld,\"min_q1\":%ld",n,t1-t0,nth,saturated,min_q1);
  const char*names[3]={"hist0","hist1","mincls"}; long*hs[3]={hist0,hist1,minh};
  for(int k=0;k<3;k++){printf(",\"%s\":[",names[k]); int first=1;
    for(int m=0;m<=MAXM;m++) if(hs[k][m]){printf("%s[%d,%ld]",first?"":",",m,hs[k][m]); first=0;} printf("]");}
  printf(",\"top\":["); int first=1;
  for(int p=0;p<G;p++) for(int k=0;k<ntops[p];k++){cm*e=tops+(size_t)p*top+k; printf("%s[%u,%u]",first?"":",",e->c,e->m); first=0;}
  printf("]}\n");
  return 0;
}
