// exact cluster editing by subset DP, n<=20.  adj: array of n bitmasks.
#include <stdlib.h>
#include <string.h>
static int N; static unsigned adjm[32]; static int *f; static long long *cnt;
static int cost_cluster(unsigned C, unsigned S){
    int inside=0,out=0; unsigned rest=S&~C;
    for(int v=0;v<N;v++) if(C>>v&1){ inside+=__builtin_popcount(adjm[v]&C); out+=__builtin_popcount(adjm[v]&rest);}
    int c=__builtin_popcount(C); return c*(c-1)/2-inside/2+out;
}
int solve(int n, unsigned *adj, int *labels_out, int maxsol, int *nsol_out){
    N=n; for(int i=0;i<n;i++) adjm[i]=adj[i];
    unsigned full=(n==32)?0xffffffffu:((1u<<n)-1);
    f=(int*)malloc(sizeof(int)*(full+1)); cnt=(long long*)malloc(sizeof(long long)*(full+1));
    f[0]=0; cnt[0]=1;
    for(unsigned S=1;S<=full;S++){
        unsigned low=S&-S, rest=S^low, sub=rest; int best=1<<30; long long c=0;
        while(1){
            unsigned C=sub|low; int v=cost_cluster(C,S)+f[S&~C];
            if(v<best){best=v;c=cnt[S&~C];} else if(v==best) c+=cnt[S&~C];
            if(sub==0)break; sub=(sub-1)&rest;
        }
        f[S]=best; cnt[S]=c;
    }
    // enumerate up to maxsol optimal partitions (iterative DFS with explicit stack)
    int nsol=0;
    // recursive via stack arrays
    unsigned stackS[40], stackC[40], stackSub[40]; int depth=0;
    unsigned cur[40]; 
    // simple recursion emulate
    struct {unsigned S,sub; int started;} st[40];
    st[0].S=full; st[0].started=0; depth=0;
    while(depth>=0 && nsol<maxsol){
        if(st[depth].S==0){
            // record
            for(int d=0;d<depth;d++){ for(int v=0;v<n;v++) if(cur[d]>>v&1) labels_out[nsol*n+v]=d; }
            nsol++; depth--; continue;
        }
        unsigned S=st[depth].S, low=S&-S, rest=S^low;
        if(!st[depth].started){ st[depth].sub=rest; st[depth].started=1; }
        else { if(st[depth].sub==0){ depth--; continue; } st[depth].sub=(st[depth].sub-1)&rest; }
        // find next valid sub
        unsigned sub=st[depth].sub; int found=0;
        while(1){
            unsigned C=sub|low;
            if(cost_cluster(C,S)+f[S&~C]==f[S]){found=1; st[depth].sub=sub; cur[depth]=C; break;}
            if(sub==0)break; sub=(sub-1)&rest;
        }
        if(!found){ depth--; continue; }
        st[depth+1].S=S&~cur[depth]; st[depth+1].started=0; depth++;
    }
    *nsol_out=nsol; int res=f[full]; free(f); free(cnt); return res;
}
