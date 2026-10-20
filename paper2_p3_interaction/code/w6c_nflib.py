import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6c_lib import *
import networkx as nx, itertools, numpy as np
from networkx.generators.atlas import graph_atlas_g
def packings_of_size(A, k):
    P3=p3_list(A); ps=[pairs_of(t) for t in P3]
    out=[]
    def rec(i, chosen, used):
        if len(chosen)==k: out.append(list(chosen)); return
        if len(P3)-i < k-len(chosen): return
        for j in range(i,len(P3)):
            if not (ps[j]&used):
                chosen.append(j); rec(j+1,chosen,used|ps[j]); chosen.pop()
    rec(0,[],set())
    return [[P3[j] for j in ch] for ch in out]
def minimal_nonfaces(singfam, anch):
    # singfam: list of frozensets (singleton sets of optimal partitions) restricted to anch
    fam=[s&frozenset(anch) for s in singfam]
    def isface(S): return any(S<=f for f in fam)
    res=[]
    for r in range(1,len(anch)+1):
        for S in itertools.combinations(anch,r):
            S=frozenset(S)
            if any(T<=S for T in res): continue
            if not isface(S): res.append(S)
    return res
