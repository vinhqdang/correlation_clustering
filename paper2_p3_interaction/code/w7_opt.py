"""Exact cluster editing by subset DP (n <= 14): f(S)=min over clusters C containing lowest vertex of S of cost(C)+f(S\\C),
cost(C)=nonedges inside C + edges between C and S\\C counted once (edges leaving to the rest of S).  Enumerates all optimal partitions."""
import itertools, numpy as np
from functools import lru_cache
def solve(n, edges, enum=True, limit=100000):
    adj=[0]*n
    for u,v in edges: adj[u]|=1<<v; adj[v]|=1<<u
    pc=lambda x: bin(x).count('1')
    full=(1<<n)-1
    @lru_cache(None)
    def f(S):
        if S==0: return 0
        low=S&-S; rest=S^low
        best=10**9
        sub=rest
        while True:
            C=sub|low
            # edits: nonedges inside C + edges from C to S\C
            inside=0; out=0
            for v in range(n):
                if C>>v&1:
                    inside+=pc(adj[v]&C); out+=pc(adj[v]&(S&~C))
            ne=(pc(C)*(pc(C)-1)//2)-inside//2
            c=ne+out+f(S&~C)
            if c<best: best=c
            if sub==0: break
            sub=(sub-1)&rest
        return best
    opt=f(full)
    sols=[]
    if enum:
        def rec(S,acc):
            if len(sols)>=limit: return
            if S==0: sols.append(tuple(acc)); return
            low=S&-S; rest=S^low; sub=rest
            while True:
                C=sub|low
                inside=0; out=0
                for v in range(n):
                    if C>>v&1:
                        inside+=pc(adj[v]&C); out+=pc(adj[v]&(S&~C))
                ne=(pc(C)*(pc(C)-1)//2)-inside//2
                if ne+out+f(S&~C)==f(S): rec(S&~C,acc+[C])
                if sub==0: break
                sub=(sub-1)&rest
        rec(full,[])
    return opt, sols
def cl(C): return [i for i in range(64) if C>>i&1]
