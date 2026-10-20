"""compare the CSP leaf oracle with the exact CP-SAT oracle (far pairs respected) on random instances."""
import sys, random, itertools, time
from w9b_lib import *
from w9b_leaf import *
from w9b_rand2 import gen as gen_planted
from w9b_rand import gen as gen_rand
def run(seed,N,m,nv,cmax,mode):
    rng=random.Random(seed); agree=dis=skip=0; widths=[]; nfar=[]
    for _ in range(N):
        I=(gen_planted(rng,m,cmax,nv,rng.choice([0.05,0.1,0.2,0.3])) if mode=='planted' else gen_rand(rng,m,cmax,nv,rng.choice([0.2,0.4])))
        if I is None: continue
        far=I.far_pairs()
        r1=opt_cpsat(I,respect=far)
        c,w=leaf_cost(I)
        e1=None if r1 is None else r1[0]
        if e1==c: agree+=1; widths.append(w); nfar.append(len(far))
        else: dis+=1; print('DISAGREE',e1,c,I.n,sorted(I.E),I.H)
    return agree,dis,max(widths) if widths else None,sum(nfar)/max(1,len(nfar))
if __name__=='__main__':
    t0=time.time()
    for (seed,N,m,nv,cmax,mode) in [(1,60,4,9,3,'planted'),(2,60,5,10,3,'planted'),(3,60,4,9,2,'rand'),(4,60,6,13,2,'planted'),(5,40,5,11,3,'rand')]:
        print(seed,N,m,nv,cmax,mode,run(seed,N,m,nv,cmax,mode),round(time.time()-t0),flush=True)
