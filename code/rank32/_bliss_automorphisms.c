/* Engine of rank32/bliss_automorphisms.py: a bounded-memory invariant of the lines of
   shortest vectors, and the coloured graph of a local Gram matrix for bliss.
   profile: input uint32 n, uint32 b[n] (the classes mod 2 of the lines)
            -> output uint64 moments[n][3].
   graph:   input uint32 n, int8 Gram[n*n]
            -> canonical labels of the n core vertices (file OUT.lab) and ALL
               automorphism generators found by bliss (file OUT.gens); counts on stdout.
   Each nonzero off-diagonal Gram value has its own layer of vertices.
   Compiled by bliss_automorphisms.py with gcc -O3 -fopenmp and linked to bliss.

   rank32/bliss_automorphisms.py のエンジン。最短ベクトルの直線の、メモリを抑えた不変量と、
   局所の Gram 行列の色付きグラフ（bliss に渡す）を作る。
   profile：入力 uint32 n、uint32 b[n]（直線の法 2 の類）-> 出力 uint64 moments[n][3]。
   graph：  入力 uint32 n、int8 Gram[n*n] -> 中心の n 頂点の標準ラベル（ファイル OUT.lab）と、
            bliss が見つけた自己同型の生成元の全部（ファイル OUT.gens）。数は標準出力に書く。
   Gram の非対角の 0 でない値ごとに、頂点の層を 1 つずつ置く。
   bliss_automorphisms.py が gcc -O3 -fopenmp でコンパイルし、bliss をリンクする。 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bliss/bliss_C.h>
#include <omp.h>
static FILE *out; static uint32_t core;
static void hook(void *p,unsigned n,const unsigned *a){
  (void)p;(void)n; if(fwrite(a,4,core,out)!=core) exit(4);
}
int main(int argc,char **argv){
  if(argc!=4) return 2;
  FILE *f=fopen(argv[2],"rb"); if(!f) return 3;
  uint32_t n; if(fread(&n,4,1,f)!=1) return 3;
  if(!strcmp(argv[1],"profile")){
    uint32_t *b=malloc(4L*n); if(fread(b,4,n,f)!=n) return 3; fclose(f);
    uint64_t *sig=calloc(3L*n,8);
    /* 4096 buckets, 2 MB/thread. No hash collisions: exact integer moments.
       4096 個のバケット、スレッドあたり 2 MB。ハッシュの衝突は無い（厳密な整数のモーメント）。 */
    #pragma omp parallel
    {
      uint16_t *t=malloc(2L*(1<<20)); uint64_t *s=calloc(3L*n,8);
      #pragma omp for schedule(dynamic)
      for(unsigned p=0;p<4096;p++){
        memset(t,0,2L*(1<<20));
        /* Bucket index lists avoid rescanning all pairs 4096 times.
           バケットごとの添字の列で、全部の組を 4096 回走査し直すのを避ける。 */
        unsigned cnt[4097]={0};
        for(unsigned i=0;i<n;i++) cnt[(b[i]>>20)+1]++;
        for(int g=0;g<4096;g++) cnt[g+1]+=cnt[g];
        unsigned *idx=malloc(4L*n), fill[4096]; memcpy(fill,cnt,sizeof fill);
        for(unsigned i=0;i<n;i++) idx[fill[b[i]>>20]++]=i;
        for(unsigned g=0;g<4096;g++) {unsigned h=g^p;if(h<g)continue;
          for(unsigned a=cnt[g];a<cnt[g+1];a++){
            unsigned i=idx[a],start=g==h?a+1:cnt[h];
            for(unsigned c=start;c<cnt[h+1];c++){
              unsigned v=(b[i]^b[idx[c]])&((1<<20)-1);
              if(t[v]==65535) exit(5); t[v]++;
            }
          }
        }
        for(unsigned i=0;i<n;i++){
          unsigned h=(b[i]>>20)^p;
          for(unsigned c=cnt[h];c<cnt[h+1];c++){
            unsigned j=idx[c]; if(i==j) continue;
            uint64_t m=t[(b[i]^b[j])&((1<<20)-1)];
            s[3L*i]+=m; s[3L*i+1]+=m*m; s[3L*i+2]+=m*m*m;
          }
        }
        free(idx);
      }
      #pragma omp critical
      {for(unsigned long i=0;i<3L*n;i++)sig[i]+=s[i];}
      free(t);free(s);
    }
    out=fopen(argv[3],"wb"); fwrite(sig,8,3L*n,out); fclose(out);free(sig);free(b);
  }else if(!strcmp(argv[1],"graph")){
    core=n; int8_t *q=malloc((size_t)n*n);
    if(fread(q,1,(size_t)n*n,f)!=(size_t)n*n) return 3; fclose(f);
    int present[9]={0},layer[9];
    for(unsigned i=0;i<n;i++)for(unsigned j=i+1;j<n;j++)if(q[(size_t)i*n+j])present[q[(size_t)i*n+j]+4]=1;
    BlissGraph *g=bliss_new(0);
    for(unsigned i=0;i<n;i++) bliss_add_vertex(g,i==0?1:0);
    unsigned count=n;
    for(int v=0;v<9;v++)if(present[v]){
      layer[v]=count;
      for(unsigned i=0;i<n;i++){bliss_add_vertex(g,v+2);bliss_add_edge(g,i,count++);}
    }
    for(unsigned i=0;i<n;i++)for(unsigned j=i+1;j<n;j++){
      int v=q[(size_t)i*n+j];if(v)bliss_add_edge(g,layer[v+4]+i,layer[v+4]+j);
    }
    char path[4096];snprintf(path,sizeof path,"%s.gens",argv[3]);out=fopen(path,"wb");
    BlissStats stats; const unsigned *lab=bliss_find_canonical_labeling(g,hook,NULL,&stats);fclose(out);
    snprintf(path,sizeof path,"%s.lab",argv[3]);out=fopen(path,"wb");fwrite(lab,4,n,out);fclose(out);
    printf("{\"nodes\":%lu,\"generators\":%lu,\"vertices\":%u}\n",stats.nof_nodes,stats.nof_generators,count);
    bliss_release(g);free(q);
  }else return 2;
  return 0;
}
