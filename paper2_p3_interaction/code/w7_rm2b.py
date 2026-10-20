import random, itertools, sys
from w7_rm2 import tree_paths, depth, nu, solve
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); fails=0; tot=0; gain=0; fails2=0
for _ in range(N):
    n=rng.randint(10,22)
    parent=[-1]+[rng.randrange(max(0,i-rng.choice([1,2,4,8])),i) for i in range(1,n)]
    k=rng.randint(6,12)
    confs=[]
    while len(confs)<k:
        s,t=rng.sample(range(n),2); confs.append((s,t))
    ends=sorted({v for c in confs for v in c})
    tun=set()
    for _ in range(rng.randint(5,11)):
        a,b=rng.sample(ends if rng.random()<0.85 else range(n),2)
        if sum(1 for t in tun if a in t)<3 and sum(1 for t in tun if b in t)<3: tun.add(tuple(sorted((a,b))))
    tun=sorted(tun)
    if len(tun)>11: tun=tun[:11]
    full=solve(parent,confs,tun)
    pk=nu(parent,confs); pe={v for c in pk for v in c}
    restr={t for t in tun if t[0] in pe or t[1] in pe}
    r1=solve(parent,confs,tun,restr)
    restr2={t for t in tun if t[0] in pe and t[1] in pe}
    r2=solve(parent,confs,tun,restr2)
    tot+=1; gain+= full<len(pk)
    if r1!=full: fails+=1; ex=(parent,confs,tun,full,r1)
    if r2!=full: fails2+=1
print('instances',tot,'tunnel gain',gain,'restricted(any endpoint in packing) fails',fails,' restricted(both endpoints in packing) fails',fails2)
if fails: print(ex)
