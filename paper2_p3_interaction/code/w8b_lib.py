"""Wave 8 (branching on far pairs): independent library. numpy + networkx only.
Instance: n, edge set, members (a,b,c) centre b.  Exact OPT by enumeration of set partitions (n<=11).
Local pairs: via a tree decomposition of I(H); X_v = bags containing all members through v.
Pair (a,b) is LOCAL iff X_a & X_b nonempty (a,b in V(H)); every pair touching V\\V(H) is non-local unless soft.
"""
import itertools, numpy as np, networkx as nx
from networkx.algorithms.approximation import treewidth_min_fill_in

def partitions(n):
    # restricted growth strings
    out=[]
    def rec(i,cur,m):
        if i==n: out.append(tuple(cur)); return
        for k in range(m+1):
            cur.append(k); rec(i+1,cur,max(m,k+1)); cur.pop()
    rec(0,[],0)
    return np.array(out,dtype=np.int8)

_PCACHE={}
def parts(n):
    if n not in _PCACHE: _PCACHE[n]=partitions(n)
    return _PCACHE[n]

def pairs_of(n): return [(i,j) for i in range(n) for j in range(i+1,n)]

class Inst:
    def __init__(self,n,edges,members):
        self.n=n; self.E=set(tuple(sorted(e)) for e in edges); self.H=[tuple(m) for m in members]
        self.pairs=pairs_of(n); self.pidx={p:i for i,p in enumerate(self.pairs)}
        self.soft=set()
        for (a,b,c) in self.H:
            for p in [(a,b),(b,c),(a,c)]: self.soft.add(tuple(sorted(p)))
    def valid(self):
        seen=set()
        for (a,b,c) in self.H:
            if not(self.adj(a,b) and self.adj(b,c) and not self.adj(a,c)): return False
            for p in [(a,b),(b,c),(a,c)]:
                p=tuple(sorted(p))
                if p in seen: return False
                seen.add(p)
        return True
    def adj(self,a,b): return tuple(sorted((a,b))) in self.E
    def pair_maximal(self):
        n=self.n
        for a,b,c in itertools.combinations(range(n),3):
            ps=[(a,b),(b,c),(a,c)]
            if any(p in self.soft for p in ps): continue
            if sum(self.adj(*p) for p in ps)==2: return False
        return True
    def IH(self):
        g=nx.Graph(); g.add_nodes_from(range(len(self.H)))
        for i,j in itertools.combinations(range(len(self.H)),2):
            if set(self.H[i])&set(self.H[j]): g.add_edge(i,j)
        return g
    def decomposition(self):
        g=self.IH()
        if len(self.H)==0: return 0,[frozenset()]
        if g.number_of_edges()==0: return 0,[frozenset([i]) for i in g.nodes]
        w,T=treewidth_min_fill_in(g)
        return w,list(T.nodes)
    def local_matrix(self):
        w,bags=self.decomposition()
        mem={v:frozenset(i for i,m in enumerate(self.H) if v in m) for v in range(self.n)}
        X={}
        for v in range(self.n):
            if not mem[v]: X[v]=frozenset()
            else: X[v]=frozenset(k for k,b in enumerate(bags) if mem[v]<=b)
        loc={}
        for (a,b) in self.pairs: loc[(a,b)]=bool(X[a]&X[b])
        self.width=w; return loc

class Oracle:
    """disagreement matrix for all partitions; F(T) = min cost over partitions respecting T-pairs, non-local non-T pairs free."""
    def __init__(self,inst):
        self.I=inst; n=inst.n; self.P=parts(n)
        pr=inst.pairs
        a=np.array([p[0] for p in pr]); b=np.array([p[1] for p in pr])
        same=(self.P[:,a]==self.P[:,b])
        self.same=same
        self.Ev=np.array([p in inst.E for p in pr])
        self.loc=inst.local_matrix()
        self.counted=np.array([(p in inst.soft) or self.loc[p] for p in pr])
        self.N=[i for i,p in enumerate(pr) if not self.counted[i]]
        self.flip=np.zeros(len(pr),bool)
    def dis(self):
        E=self.Ev^self.flip
        return same_ne(self.same,E)
def same_ne(same,E): return same!=E[None,:]
