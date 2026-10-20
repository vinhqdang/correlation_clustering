"""Large-scale instance check of Theorem A (soft-triangle-free) and of the contrast (soft triangles present)."""
import random, sys, itertools
import networkx as nx
from w7_lib import *
from w7_cores import cores
from w7_rel7 import closed_all_majorities
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2])
cnt={'tf':[0,0,0,0],'tri':[0,0,0,0]}   # instances, closed under all majorities, 2-decomposable, instances with a core of size>=3
forest=0
ntry=0
while sum(c[0] for c in cnt.values())<N:
    ntry+=1
    m=rng.randint(3,7); n=rng.randint(int(1.4*m)+2,3*m+1)
    I=gen_planted(n,m,rng,cmax=rng.choice([2,2,3,4]),maxcl=rng.choice([3,4,5,6]),tries=1500)
    if I is None: continue
    sg=I.sigma()
    if len(sg)<2 or len(sg)>80: continue
    k='tri' if I.soft_triangles() else 'tf'
    if k=='tf' and nx.is_forest(I.Ig()): forest+=1
    cnt[k][0]+=1
    cnt[k][1]+=closed_all_majorities(sg)
    cnt[k][2]+=two_decomposable(sg)
    an,co=cores(I,sg)
    cnt[k][3]+= bool(co and max(len(c) for c in co)>=3)
print('tries',ntry,'forest instances among tf:',forest)
for k,v in cnt.items(): print(k,'instances',v[0],'closed-all-majorities',v[1],'2-decomposable',v[2],'with core>=3',v[3])
