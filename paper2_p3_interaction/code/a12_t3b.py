import sys,random
from a12_t3 import *
rng=random.Random(int(sys.argv[1])); N=int(sys.argv[2]); l=int(sys.argv[3])
c={'tot':0,'yes0':0,'yes_need_far':0,'no':0,'mism':0,'cyc':0}
for _ in range(N):
    n=rng.randint(7,11); E,H=gen(rng,n,5,rng.choice([.2,.3,.4]))
    ins=Inst(n,E,H); st={}
    a=solveU(ins,frozenset(),l,st); b=brute(ins,l); b0=brute(ins,0)
    c['tot']+=1; c['yes0']+=b0; c['yes_need_far']+=(b and not b0); c['no']+=(not b); c['mism']+=(a!=b)
    ig=I_edges=[(i,j) for i in range(len(H)) for j in range(i+1,len(H)) if set(H[i])&set(H[j])]
    if len(ig)>=len(H): c['cyc']+=1
    if a!=b: print('MISMATCH',n,sorted(E),H)
print(c)
