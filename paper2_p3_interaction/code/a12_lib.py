"""Audit-12 library, written from scratch (itertools/random only)."""
import itertools, random

def K(a,b): return (a,b) if a<b else (b,a)

def opt_dp(n,E):
    """exact cluster-editing OPT by subset DP: cost = |E| + sum_clusters (nonedges(S)-edges(S))."""
    E=set(K(*e) for e in E)
    adj=[0]*n
    for a,b in E: adj[a]|=1<<b; adj[b]|=1<<a
    N=1<<n
    w=[0]*N   # nonedges - edges inside S
    for S in range(1,N):
        v=(S&-S).bit_length()-1; R=S&(S-1)
        # pairs (v,u) u in R
        ne=bin(adj[v]&R).count('1'); nn=bin(R).count('1')-ne
        w[S]=w[R]+nn-ne
    INF=10**9; f=[INF]*N; f[0]=0
    for S in range(1,N):
        low=S&-S; rest=S^low; sub=rest; best=INF
        while True:
            C=sub|low
            c=f[S^C]+w[C]
            if c<best: best=c
            if sub==0: break
            sub=(sub-1)&rest
        f[S]=best
    return len(E)+f[N-1]

def is_cluster(n,E):
    adj=[set() for _ in range(n)]
    for a,b in E: adj[a].add(b); adj[b].add(a)
    for b in range(n):
        nb=sorted(adj[b])
        for i in range(len(nb)):
            for j in range(i+1,len(nb)):
                if nb[j] not in adj[nb[i]]: return False
    return True

def p3s(n,E):
    E=set(E); out=[]
    for b in range(n):
        nb=[v for v in range(n) if K(v,b) in E]
        for i in range(len(nb)):
            for j in range(i+1,len(nb)):
                if K(nb[i],nb[j]) not in E: out.append((nb[i],b,nb[j]))
    return out

def prs(t):
    a,b,c=t; return [K(a,b),K(b,c),K(a,c)]

def rand_graph(n,p,rng):
    return set(K(a,b) for a in range(n) for b in range(a+1,n) if rng.random()<p)

def rand_packing(n,E,rng,maximal=True,drop=0.0):
    P=p3s(n,E); rng.shuffle(P); used=set(); H=[]
    for t in P:
        if all(p not in used for p in prs(t)):
            H.append(t); used.update(prs(t))
    if not maximal:
        H=[t for t in H if rng.random()>drop] or H[:1]
    return H

def soft_set(H): return {p for t in H for p in prs(t)}

def has_soft_triangle(H):
    own={}
    for i,t in enumerate(H):
        for p in prs(t): own[p]=i
    V=sorted({v for t in H for v in t})
    for a,b,c in itertools.combinations(range(max(V)+1),3):
        o=[own.get(K(a,b)),own.get(K(b,c)),own.get(K(a,c))]
        if None not in o and len(set(o))==3: return True
    return False

def cross_pairs(H):
    """hard-or-not pairs {p,r}: p in V(A)\\q, r in V(B)\\q for distinct members A,B sharing q"""
    S=set()
    for i in range(len(H)):
        for j in range(i+1,len(H)):
            sh=set(H[i])&set(H[j])
            for q in sh:
                for p in set(H[i])-{q}:
                    for r in set(H[j])-{q}:
                        if p!=r: S.add(K(p,r))
    return S

# states on member (a,b,c): final adjacency (y_ab,y_bc,y_ac); edit set is y xor (1,1,0)
ST={'L':(0,1,0),'R':(1,0,0),'M':(1,1,1),'O':(0,0,0),'X':(0,0,1)}
def edits_of(y):
    base=(1,1,0); return [i for i in range(3) if y[i]!=base[i]]
def apply_state(E,t,y):
    E=set(E); P=prs(t)
    for i in range(3):
        if (K(*P[i]) in E)!=bool(y[i]):
            E.symmetric_difference_update({P[i]})
    return E
