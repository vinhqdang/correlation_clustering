import itertools, random, sys, io, contextlib
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from w7_lib import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2])
# prism: members 0,1,2 form triangle (shared vertices t01,t02,t12), members 3,4,5 triangle (t34,t35,t45), matching vertices m03,m14,m25
names=['t01','t02','t12','t34','t35','t45','m03','m14','m25']; vid={x:i for i,x in enumerate(names)}
mem=[['t01','t02','m03'],['t01','t12','m14'],['t02','t12','m25'],['t34','t35','m03'],['t34','t45','m14'],['t35','t45','m25']]
mem=[[vid[x] for x in m] for m in mem]
soft=set()
for m in mem:
    for a,b in itertools.combinations(m,2): soft.add(key(a,b))
hard=[p for p in itertools.combinations(range(9),2) if p not in soft]
print('prism hard pairs',len(hard))
tot=0; nd=0; nclose=0; ex=None
tries=0
while tot<N:
    tries+=1
    orders=[tuple(rng.sample(m,3)) for m in mem]
    A=np.zeros((9,9),bool)
    for (a,b,c) in orders: A[a,b]=A[b,a]=A[b,c]=A[c,b]=True
    # sparse random hard edges (dense random graphs are rarely tight)
    pr=rng.choice([0.05,0.1,0.2,0.3])
    for p in hard:
        if rng.random()<pr: A[p[0],p[1]]=A[p[1],p[0]]=True
    I=Inst(A,orders)
    if not I.valid(): continue
    sg=I.sigma()
    if len(sg)<2: continue
    tot+=1
    if not two_decomposable(sg): nd+=1; ex=(orders,A.astype(int).tolist(),sg)
print('prism: tight instances with |Sigma|>=2:',tot,'(tries %d)'%tries,'non-2-decomposable:',nd)
if ex: print(ex)
