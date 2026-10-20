import random, itertools, sys, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from w7_rm2 import tree_paths, depth, nu, solve
def all_max_packings(parent,confs):
    P=[tree_paths(parent,s,t)[0] for s,t in confs]
    k=len(confs); best=len(nu(parent,confs)); res=[]
    for sub in itertools.combinations(range(k),best):
        ok=True; used=set()
        for i in sub:
            if P[i]&used: ok=False;break
            used|=P[i]
        if ok: res.append([confs[i] for i in sub])
    return res
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); fails=0; tot=0; gain=0; npk=0
for _ in range(N):
    n=rng.randint(8,18)
    parent=[-1]+[rng.randrange(max(0,i-rng.choice([1,2,4,8])),i) for i in range(1,n)]
    k=rng.randint(5,10)
    confs=[]
    while len(confs)<k:
        s,t=rng.sample(range(n),2); confs.append((s,t))
    ends=sorted({v for c in confs for v in c})
    tun=set()
    for _ in range(rng.randint(5,10)):
        a,b=rng.sample(ends,2)
        if sum(1 for t in tun if a in t)<3 and sum(1 for t in tun if b in t)<3: tun.add(tuple(sorted((a,b))))
    tun=sorted(tun)
    full=solve(parent,confs,tun)
    pks=all_max_packings(parent,confs)
    tot+=1; gain+= full<len(nu(parent,confs))
    for pk in pks:
        npk+=1
        pe={v for c in pk for v in c}
        restr={t for t in tun if t[0] in pe or t[1] in pe}
        if solve(parent,confs,tun,restr)!=full: fails+=1; ex=(parent,confs,tun,pk,full)
print('instances',tot,'with gain',gain,'(instance,max packing) pairs',npk,'failures',fails)
if fails: print(ex)
