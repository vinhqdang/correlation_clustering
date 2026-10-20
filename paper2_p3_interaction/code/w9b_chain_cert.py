"""K3-chain with pendant guards: unbounded minimal certificate (Core Bound counterexample at l=0, c=2, tw(I(H))=2)."""
import itertools, sys, time
from w9b_lib import *
from w9b_leaf import leaf_cost
from w9b_k3chain import A0,B0,C0,HARD0

def build(N):
    """returns inst, T (list of far pairs forming the certificate), info"""
    vid={}
    def V(j,l):
        if l==2 and j>0: return V(j-1,4)
        if (j,l) not in vid: vid[(j,l)]=len(vid)
        return vid[(j,l)]
    H=[];E=set();soft=set()
    for j in range(N):
        for t in (A0,B0,C0): H.append(tuple(V(j,x) for x in t))
    for (a,b,c) in H: E.add(key(a,b));E.add(key(b,c)); soft|={key(a,b),key(b,c),key(a,c)}
    for j in range(N):
        for (x,y),s in HARD0.items():
            u,v=V(j,x),V(j,y)
            if key(u,v) in soft: continue
            if s: E.add(key(u,v))
    nv=len(vid)
    T=[]
    pend=[]
    def add_pendant(adj_to, nonadj_to):
        nonlocal nv
        w=nv; nv+=1
        E.add(key(w,adj_to))
        T.append(key(w,adj_to)); T.append(key(w,nonadj_to))
        pend.append(w)
    # guard {L,R} on C_j=(p,c,r): pair (p,r) apart : w adjacent to p only; far pairs (w,p),(w,r)
    for j in range(N): add_pendant(V(j,0),V(j,3))
    # A_1 = R : pair (p,a) apart in A=(q,p,a)... guard x_bc=0 on A_1: pair (p_1,a_1): w adjacent to a_1 only, far pairs (w,a),(w,p)
    add_pendant(V(0,2),V(0,0))
    # B_N = R : pair (r,b) apart in B=(q,r,b): w adjacent to b_N only, far pairs (w,b),(w,r)
    add_pendant(V(N-1,4),V(N-1,3))
    I=Inst(nv,E,H)
    return I,T

if __name__=='__main__':
    for N in range(1,int(sys.argv[1])+1):
        I,T=build(N)
        assert I.valid(), 'invalid'
        far=set(I.far_pairs())
        assert all(p in far for p in T)
        t0=time.time()
        base=len(I.H)
        c_full,w=leaf_cost(I,T)
        c_empty,_=leaf_cost(I,[])
        res=[]
        for e in T:
            c,_=leaf_cost(I,[p for p in T if p!=e]); res.append(c)
        # independent oracle
        o_full=opt_csp_cpsat(I,T); o_empty=opt_csp_cpsat(I,[])
        o_res=[opt_csp_cpsat(I,[p for p in T if p!=e]) for e in T]
        print('N',N,'members',len(I.H),'n',I.n,'c',I.cmax(),'pair-maximal',I.pair_maximal(),'soft triangles',len(I.soft_triangles()),'|T|',len(T),
              '| bucket: F(T)=',c_full,'F(empty)=',c_empty,'F(T-e) values',sorted(set(res)),
              '| cpsat: F(T)=',o_full,'F(empty)=',o_empty,'F(T-e)',sorted(set(o_res)),round(time.time()-t0,1),flush=True)
