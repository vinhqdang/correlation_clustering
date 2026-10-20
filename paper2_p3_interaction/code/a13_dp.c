// a13_dp.c  exact min over partitions of vertex set: counted disagreements, with hard constraints from T.
// input: n; then n lines of n ints: G adjacency; n lines: counted (0/1); n lines: T code (0 none,1 T-edge must be joined,2 T-nonedge must be separated)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(){int n;scanf("%d",&n);int A[24][24],C[24][24],T[24][24];
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)scanf("%d",&A[i][j]);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)scanf("%d",&C[i][j]);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)scanf("%d",&T[i][j]);
 unsigned cm[24]={0},ce[24]={0},te[24]={0},tn[24]={0};
 long Etot=0;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++){if(i==j)continue;
  if(C[i][j])cm[i]|=1u<<j; if(C[i][j]&&A[i][j])ce[i]|=1u<<j;
  if(T[i][j]==1)te[i]|=1u<<j; if(T[i][j]==2)tn[i]|=1u<<j;
  if(i<j&&C[i][j]&&A[i][j])Etot++;}
 unsigned full=(n==32)?0xffffffffu:((1u<<n)-1); size_t M=(size_t)1<<n;
 int *g=malloc(M*sizeof(int)); unsigned *teu=malloc(M*4); unsigned char *viol=malloc(M);
 int *pin=malloc(M*sizeof(int)),*ein=malloc(M*sizeof(int));
 int *f=malloc(M*sizeof(int));
 pin[0]=ein[0]=0;teu[0]=0;viol[0]=0;
 const int INF=1<<28;
 g[0]=0;
 for(size_t S=1;S<M;S++){int low=__builtin_ctz(S);unsigned rest=S&(S-1);
  pin[S]=pin[rest]+__builtin_popcount(cm[low]&rest);
  ein[S]=ein[rest]+__builtin_popcount(ce[low]&rest);
  teu[S]=teu[rest]|te[low];
  viol[S]=viol[rest]|((tn[low]&S)!=0);
  if(viol[S]||(teu[S]&~S))g[S]=INF; else g[S]=pin[S]-2*ein[S];}
 f[0]=0;
 for(size_t S=1;S<M;S++){unsigned low=S&-S;unsigned rest=S^low;int best=INF;
  // enumerate sub subsets of rest
  unsigned sub=rest;
  while(1){unsigned Cc=sub|low; if(g[Cc]<INF){int v=g[Cc]+f[S^Cc]; if(v<best)best=v;}
   if(sub==0)break; sub=(sub-1)&rest;}
  f[S]=best;}
 printf("%ld\n",Etot+f[full]);return 0;}
