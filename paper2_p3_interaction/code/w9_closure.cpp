// finite closure of rooted-subtree profiles (soft-triangle-free forest CSP, states = partitions of the 3 vertices of a member)
#include <bits/stdc++.h>
using namespace std;
int CAP=3; const int INF=9;
// states: 0=L 1=R 2=M 3=O 4=X ; y=(ab,bc,ac)
int Y[5][3]={{0,1,0},{1,0,0},{1,1,1},{0,0,0},{0,0,1}};
int XC[5]={0,0,0,1,2};
int viewIdx(int s,int role){ // role 0=a:(ab,ac) 1=b:(ab,bc) 2=c:(bc,ac); index = y1*2+y2 matching VIEWS order (0,0),(0,1),(1,0),(1,1)
    int a,b;
    if(role==0){a=Y[s][0];b=Y[s][2];}
    else if(role==1){a=Y[s][0];b=Y[s][1];}
    else {a=Y[s][1];b=Y[s][2];}
    return a*2+b;
}
// lists as bitmasks over states
vector<int> LISTS={0b11111,0b01001 /*L,O*/,0b11000 /*O,X*/,0b01010 /*R,O*/,0b01000,0b00001,0b00010,0b10000,0b00100,0};
bool compat(int vu,int vv,int eps){ // vu,vv view idx; eps bits: eps[i][j] at bit (i*2+j)
    for(int i=0;i<2;i++)for(int j=0;j<2;j++){
        int ui=(vu>>(1-i))&1, vj=(vv>>(1-j))&1;
        int e=(eps>>(i*2+j))&1;
        if(e){ if(ui!=vj) return false;} else { if(ui&&vj) return false; }
    }
    return true;
}
typedef array<int,4> F4;
int addc(int a,int b){ if(a>=INF||b>=INF) return INF; int s=a+b; return s<=CAP?s:CAP+1; }
F4 normF(F4 p){ int m=INF; for(int x:p) m=min(m,x); if(m>=INF) return p; F4 r; for(int i=0;i<4;i++) r[i]=p[i]>=INF?INF:p[i]-m; return r; }
F4 induced(const F4&pr,int eps){ F4 out; for(int vu=0;vu<4;vu++){ int best=INF; for(int vv=0;vv<4;vv++) if(pr[vv]<INF&&compat(vu,vv,eps)) best=min(best,pr[vv]); out[vu]=best;} return normF(out); }
int main(int argc,char**argv){
    if(argc>1) CAP=atoi(argv[1]);
    set<F4> IND, PROFS; 
    F4 noch={0,0,0,0}; IND.insert(noch);
    int round=0;
    while(true){
        round++;
        bool added=false;
        vector<F4> ind(IND.begin(),IND.end());
        for(int role=0;role<3;role++){
            int o0=(role==0)?1:0, o1=(role==2)?1:2;
            for(int lst:LISTS)
                for(auto&f1:ind) for(auto&f2:ind){
                    F4 prof={INF,INF,INF,INF};
                    for(int s=0;s<5;s++){ if(!(lst>>s&1)) continue;
                        int c=XC[s];
                        c=addc(c,f1[viewIdx(s,o0)]); c=addc(c,f2[viewIdx(s,o1)]);
                        if(c>=INF) continue;
                        int w=viewIdx(s,role); if(c<prof[w]) prof[w]=c;
                    }
                    F4 n=normF(prof);
                    if(!PROFS.count(n)){PROFS.insert(n);added=true;}
                }
        }
        bool nw=false;
        for(auto&pr:PROFS) for(int eps=0;eps<16;eps++){ F4 f=induced(pr,eps); if(!IND.count(f)){IND.insert(f);nw=true;} }
        fprintf(stderr,"round %d profiles %zu induced %zu\n",round,PROFS.size(),IND.size());
        if(!added&&!nw) break;
    }
    printf("CAP %d profiles %zu induced %zu\n",CAP,PROFS.size(),IND.size());
    map<pair<int,int>,int> lv;
    for(auto&p:PROFS){ int mx=-1,ninf=0; for(int x:p){ if(x>=INF) ninf++; else mx=max(mx,x);} lv[{mx,ninf}]++; }
    for(auto&kv:lv) printf("maxfinite %d ninf %d : %d\n",kv.first.first,kv.first.second,kv.second);
    if(argc>2){ for(auto&p:PROFS){ printf("%d %d %d %d\n",p[0],p[1],p[2],p[3]); } }
}
