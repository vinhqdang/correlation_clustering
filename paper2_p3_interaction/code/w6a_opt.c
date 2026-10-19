// exact cluster editing OPT by subset DP, n<=22.  stdin: n then n adjacency bitmasks (decimal).  stdout: OPT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(){int n;scanf("%d",&n);unsigned adj[32];for(int i=0;i<n;i++)scanf("%u",&adj[i]);
 unsigned N=1u<<n; int *e=malloc(N*sizeof(int)); int *best=malloc(N*sizeof(int)); int *f=malloc(N*sizeof(int));
 e[0]=0;for(unsigned m=1;m<N;m++){int v=__builtin_ctz(m);unsigned r=m&(m-1);e[m]=e[r]+__builtin_popcount(adj[v]&r);}
 for(unsigned m=1;m<N;m++){int s=__builtin_popcount(m);f[m]=s*(s-1)/2-2*e[m];}
 best[0]=0;
 for(unsigned m=1;m<N;m++){unsigned low=m&-m,rest=m^low,sub=rest;int bv=1<<29;
  while(1){unsigned c=sub|low;int v=f[c]+best[m^c];if(v<bv)bv=v;if(sub==0)break;sub=(sub-1)&rest;}
  best[m]=bv;}
 printf("%d\n",e[N-1]+best[N-1]);return 0;}
