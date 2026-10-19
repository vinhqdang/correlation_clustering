"""Inspect the 5-member flip gadget found by w6c_flip4.py (n=8, anchor x=7)."""
import itertools, numpy as np
from w6c_lib import *
names='LRM'
n=8
members=[(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)]
x=7
lab,same,P=part_table(n); pidx={p:i for i,p in enumerate(P)}
soft=set(); fe=set()
for m in members: soft|=pairs_of(m); fe|={key(m[0],m[1]),key(m[1],m[2])}
hard=[p for p in P if p not in soft and x not in p]
j=10
E=set(fe)
for i,p in enumerate(hard):
    if (j>>i)&1: E.add(p)
print('hard edges (non-member):',sorted(p for p in E if p not in fe))
print('member edges:',sorted(fe))
A=adj_from_edges(n,E)
print('valid packing:',check_packing(A,members),'pair-maximal:',is_pair_maximal(A,members))
opt,labs,sm,cost,PP=brute_opt(A)
print('OPT',opt,'|H|',len(members))
def sig(m,l):
    a,b,c=m; s=(l[a]==l[b],l[b]==l[c],l[a]==l[c])
    return {(0,1,0):'L',(1,0,0):'R',(1,1,1):'M'}.get(tuple(int(t) for t in s),'?')
for l,c in zip(labs,cost):
    if c==opt:
        cl={}
        for v,t in enumerate(l): cl.setdefault(t,[]).append(v)
        print('clusters',[ ''.join(map(str,b)) for b in cl.values()],' states',[sig(m,l) for m in members],' x singleton:',len(cl[l[x]])==1)
