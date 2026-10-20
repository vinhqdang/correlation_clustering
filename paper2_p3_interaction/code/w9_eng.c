// Wave 9 engine: soft-triangle-free CSP tables, forest DP, weighted region packing.
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAXN 64
#define MAXM 20
#define INF 100000
static const int VAL[5]={1,2,4,3,7};   // L R M O X as edit masks (bit0=ab, bit1=bc, bit2=ac)
static const int COSTX[5]={0,0,0,1,2}; // extra cost
static int n,m;
static int H[MAXM][3];
static unsigned char A[MAXN][MAXN];
static int own[MAXN][MAXN]; // member*3+k or -1
static unsigned char un[MAXM][5];
static unsigned char bin[MAXM][MAXM][5][5]; // valid only if adjacent (shared vertex)
static int adjm[MAXM][MAXM];
static int pairk(int mem,int u,int v){ // pair index within member: 0=ab 1=bc 2=ac
    int a=H[mem][0],b=H[mem][1],c=H[mem][2];
    if((u==a&&v==b)||(u==b&&v==a)) return 0;
    if((u==b&&v==c)||(u==c&&v==b)) return 1;
    return 2;
}
// returns 0 ok, -1 static violation, -2 soft triangle
static int build(void){
    memset(own,-1,sizeof own);
    for(int i=0;i<m;i++){
        int a=H[i][0],b=H[i][1],c=H[i][2];
        own[a][b]=own[b][a]=i*3; own[b][c]=own[c][b]=i*3+1; own[a][c]=own[c][a]=i*3+2;
    }
    memset(un,1,sizeof un); memset(bin,1,sizeof bin); memset(adjm,0,sizeof adjm);
    for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){
        int sh=0;
        for(int x=0;x<3;x++)for(int y=0;y<3;y++) if(H[i][x]==H[j][y]) sh++;
        if(sh) adjm[i][j]=adjm[j][i]=1;
    }
    for(int p=0;p<n;p++)for(int q=p+1;q<n;q++)for(int r=q+1;r<n;r++){
        int u[3]={p,q,r};
        int pr[3][2]={{p,q},{q,r},{p,r}};
        int base[3], ow[3];
        int owners[3],no=0;
        for(int k=0;k<3;k++){
            base[k]=A[pr[k][0]][pr[k][1]];
            ow[k]=own[pr[k][0]][pr[k][1]];
            if(ow[k]>=0){ int mm=ow[k]/3; int f=0; for(int z=0;z<no;z++) if(owners[z]==mm) f=1; if(!f) owners[no++]=mm; }
        }
        if(no==0){ if(base[0]+base[1]+base[2]==2) return -1; continue; }
        if(no==3) return -2;
        if(no==1){
            int mm=owners[0];
            for(int s=0;s<5;s++){
                int mk=VAL[s]; int sum=0;
                for(int k=0;k<3;k++){ int a=base[k]; if(ow[k]>=0){ int pk=ow[k]%3; if(mk>>pk&1) a^=1; } sum+=a; }
                if(sum==2) un[mm][s]=0;
            }
        } else {
            int mA=owners[0], mB=owners[1];
            for(int s=0;s<5;s++)for(int t=0;t<5;t++){
                int mkA=VAL[s], mkB=VAL[t]; int sum=0;
                for(int k=0;k<3;k++){ int a=base[k]; if(ow[k]>=0){ int mm=ow[k]/3, pk=ow[k]%3; int mk=(mm==mA)?mkA:mkB; if(mk>>pk&1) a^=1; } sum+=a; }
                if(sum==2){ bin[mA][mB][s][t]=0; bin[mB][mA][t][s]=0; }
            }
        }
    }
    return 0;
}

static int ee[1<<MAXM];
static int dpv[MAXM][5];
static int subset_cost(int R){
    // forest DP over members in R. returns INF if infeasible.
    int visited=0, order[MAXM], par[MAXM], tot=0;
    for(int root=0;root<m;root++){
        if(!(R>>root&1) || (visited>>root&1)) continue;
        int cnt=0; order[cnt]=root; par[root]=-1; visited|=1<<root; int head=0; cnt=1;
        while(head<cnt){ int v=order[head++]; for(int w=0;w<m;w++) if((R>>w&1)&&adjm[v][w]&&!(visited>>w&1)){ visited|=1<<w; par[w]=v; order[cnt++]=w; } }
        for(int idx=cnt-1;idx>=0;idx--){
            int v=order[idx];
            for(int s=0;s<5;s++){
                if(!un[v][s]){ dpv[v][s]=INF; continue; }
                dpv[v][s]=COSTX[s];
            }
            // children already folded into dpv[v] when processed (we fold child into parent below)
        }
        // second pass: fold children bottom-up properly
        int acc[MAXM][5];
        for(int idx=0;idx<cnt;idx++){ int v=order[idx]; for(int s=0;s<5;s++) acc[v][s]=un[v][s]?COSTX[s]:INF; }
        for(int idx=cnt-1;idx>0;idx--){
            int v=order[idx], p=par[v];
            for(int s=0;s<5;s++){
                if(acc[p][s]>=INF) continue;
                int best=INF;
                for(int t=0;t<5;t++) if(bin[p][v][s][t] && acc[v][t]<best) best=acc[v][t];
                if(best>=INF) acc[p][s]=INF; else acc[p][s]+=best;
            }
        }
        int best=INF; for(int s=0;s<5;s++) if(acc[root][s]<best) best=acc[root][s];
        if(best>=INF) return INF;
        tot+=best;
    }
    return tot;
}
static int fpack[1<<MAXM];
int eval(int n_,const unsigned char *adj,int m_,const int *Hin,int *out){
    n=n_; m=m_;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++) A[i][j]=adj[i*n+j];
    for(int i=0;i<m;i++)for(int k=0;k<3;k++) H[i][k]=Hin[i*3+k];
    int rc=build(); if(rc) return rc;
    int full=(1<<m)-1;
    for(int R=1;R<=full;R++) ee[R]=subset_cost(R);
    fpack[0]=0;
    for(int av=1;av<=full;av++){
        int low=av&-av; int best=fpack[av&~low];
        for(int sub=av;sub;sub=(sub-1)&av){
            if(!(sub&low)) continue;
            int v=ee[sub]; if(v>=INF) v=INF; 
            if(v>0){ int c=v+fpack[av&~sub]; if(c>best) best=c; }
        }
        fpack[av]=best;
    }
    out[0]=ee[full]; out[1]=fpack[full];
    return 0;
}
int get_e(int R){ return ee[R]; }

static int subset_cost_relaxed(int R,int relax){
    unsigned char save[MAXM][5];
    memcpy(save,un,sizeof un);
    for(int i=0;i<m;i++) if(relax>>i&1) for(int s=0;s<5;s++) un[i][s]=1;
    int v=subset_cost(R);
    memcpy(un,save,sizeof un);
    return v;
}
int get_er(int R,int relax){ return subset_cost_relaxed(R,relax); }

// light: build tables only; then cost of full set with relax masks
int build_only(int n_,const unsigned char *adj,int m_,const int *Hin){
    n=n_; m=m_;
    for(int i=0;i<n;i++)for(int j=0;j<n;j++) A[i][j]=adj[i*n+j];
    for(int i=0;i<m;i++)for(int k=0;k<3;k++) H[i][k]=Hin[i*3+k];
    return build();
}
int cost_full_relaxed(int relax){ return subset_cost_relaxed((1<<m)-1,relax); }
// number of single-member relaxations that reduce the full cost; fills out[0..m-1] with reduced cost
int single_relax(int *out){
    int base=subset_cost((1<<m)-1);
    int cnt=0;
    for(int a=0;a<m;a++){ int v=subset_cost_relaxed((1<<m)-1,1<<a); out[a]=v; if(v<base) cnt++; }
    return base;
}
