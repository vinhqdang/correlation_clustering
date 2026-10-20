import random, sys, itertools
from w7_lib import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2])
stat={'tf':[0,0,0],'tri':[0,0,0]}  # count instances (|Sigma|>=3), closed M-maj, 2-dec
perms=list(itertools.permutations(ST))
def rule_maj(rule):
    def m(a,b,c):
        out=[]
        for x,y,z in zip(a,b,c):
            if x==y or x==z: out.append(x)
            elif y==z: out.append(y)
            else: out.append(rule[(x,y,z)])
        return ''.join(out)
    return m
def closed_rule(vecs,m):
    vs=set(vecs)
    for a in vecs:
        for b in vecs:
            for c in vecs:
                if m(a,b,c) not in vs: return False
    return True
n_try=0; ex_fail=None
while n_try<N:
    n=rng.randint(7,12); m=rng.randint(3,6)
    I=gen_planted(n,m,rng,cmax=rng.choice([2,3]),maxcl=rng.choice([3,4,5]))
    if I is None: continue
    sg=I.sigma()
    if len(sg)<3 or len(sg)>60: continue
    n_try+=1
    tri=bool(I.soft_triangles())
    k='tri' if tri else 'tf'
    stat[k][0]+=1
    rule={p:rng.choice(ST) for p in perms}
    ok=closed_mmaj_fast(sg) and closed_rule(sg,rule_maj(rule))
    stat[k][1]+=ok
    stat[k][2]+=two_decomposable(sg)
    if (not tri) and not ok and ex_fail is None: ex_fail=(I.A.astype(int).tolist(),I.H,sg)
print(stat,'FAIL example (soft-triangle-free, not closed):',ex_fail)
