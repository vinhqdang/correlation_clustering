"""Audit 14: independent re-measurement on the THICK instance of the reduction (built with w10_red.build_gt):
validity of the packing (induced P3, pair-disjoint), pair-maximality (no all-hard induced P3), c, u_max (soft pairs inside a super-vertex),
soft triangles, components of I(H), and a pathwidth UPPER BOUND of I(H) from an explicit linear order (vertex separation), for several n."""
import sys, os, random, itertools, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w10_red import build_gt, layout
from w10_lib import Inst

def k2(a,b): return (a,b) if a<b else (b,a)

def cellpos(name, m):
    """position (cell index) of a vertex name"""
    # chain vertex: (tag, ('s',k), i) or (tag, ('e',j), i);  guests: ('z',tag,u), ('y',tag,u), ('x',tag,u), ('c',c1,c2,alpha)
    if name[0] in ('z','y','x'): return 2*name[2]+1
    if name[0]=='c': return None
    kind,idx=name[1]
    return idx//2 if kind=='s' else idx
def analyse(k,n,p,seed=1,dens=0.6):
    rng=random.Random(seed)
    S=[[set() for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for j in range(k):
            for a in range(1,n+1):
                for b in range(1,n+1):
                    if rng.random()<dens: S[i][j].add((a,b))
            if not S[i][j]: S[i][j].add((1,1))
    t0=time.time()
    W,meta=build_gt(k,n,S,p,thin_mode=False)
    I=Inst(W); N=I.n; E=I.E; mem=I.members
    inv={v:nm for nm,v in I.idx.items()}
    soft=set(); ok=True
    for (a,b,c) in mem:
        ok &= (k2(a,b) in E and k2(b,c) in E and k2(a,c) not in E)
        for q in (k2(a,b),k2(b,c),k2(a,c)):
            if q in soft: ok=False
            soft.add(q)
    hardnb=collections.defaultdict(set)
    for q in E:
        if q not in soft: hardnb[q[0]].add(q[1]); hardnb[q[1]].add(q[0])
    allhard=0
    for v,nb in hardnb.items():
        nb=sorted(nb)
        for i in range(len(nb)):
            for j in range(i+1,len(nb)):
                q=k2(nb[i],nb[j])
                if q not in E and q not in soft: allhard+=1
    cnt=collections.Counter()
    for t in mem:
        for v in t: cnt[v]+=1
    c=max(cnt.values())
    # u_max: soft pairs inside a super-vertex (same (tag,thinvertex))
    def sv(v):
        nm=inv[v]
        return (nm[0],nm[1]) if len(nm)==3 and isinstance(nm[1],tuple) else None
    inside=collections.Counter()
    for q in soft:
        a,b=sv(q[0]),sv(q[1])
        if a is not None and a==b: inside[q[0]]+=1; inside[q[1]]+=1
    umax=max(inside.values()) if inside else 0
    # I(H) components & vertex separation
    byv=collections.defaultdict(list)
    for h,t in enumerate(mem):
        for v in t: byv[v].append(h)
    adj=[set() for _ in mem]
    for v,L in byv.items():
        for a in L:
            adj[a].update(L)
    for a in range(len(mem)): adj[a].discard(a)
    seen=[False]*len(mem); comps=[]
    for s in range(len(mem)):
        if not seen[s]:
            st=[s]; seen[s]=True; size=0
            while st:
                x=st.pop(); size+=1
                for y in adj[x]:
                    if not seen[y]: seen[y]=True; st.append(y)
            comps.append(size)
    # order of members by key = max cell position of its vertices (comparator guests: cell of separator)
    def key(h):
        ps=[]
        for v in mem[h]:
            nm=inv[v]; cp=cellpos(nm,0) if nm[0]!='c' else None
            if nm[0]=='c':
                cp=None
            if cp is not None: ps.append(cp)
        return max(ps) if ps else 0
    U,m,sg,sep=layout(n)
    order=sorted(range(len(mem)),key=lambda h:(key(h),h))
    rank={h:i for i,h in enumerate(order)}
    # vertex separation of I(H) under this order
    last={h:max([rank[y] for y in adj[h]]+[rank[h]]) for h in range(len(mem))}
    delta=collections.Counter()
    for h in range(len(mem)):
        if last[h]>rank[h]:
            delta[rank[h]]+=1; delta[last[h]]-=1
    cur=0; best=0
    for i in range(len(mem)):
        cur+=delta[i]; best=max(best,cur)
    print(f'k={k} n={n} p={p}: |V|={N} |H|={len(mem)} valid={ok} pair-maximal(all-hard P3 count)={allhard} c={c} 4p+2={4*p+2} u_max={umax} '
          f'components(I(H))={len(comps)} sizes={sorted(comps)[:4]}.. vertex-sep UB on pw(I(H))={best} (22k^2p^2={22*k*k*p*p}) {round(time.time()-t0,1)}s',flush=True)
if __name__=='__main__':
    k,n,p=map(int,sys.argv[1:4]); analyse(k,n,p)
