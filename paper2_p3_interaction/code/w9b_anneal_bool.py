import sys, random, itertools, time, math
import numpy as np
from w9b_lib import *
from w9b_octtree import enum_sigma
from w9b_boolproj import nonbij_triples
import w9b_anneal as A
def score(I):
    if not I.valid() or I.cmax()>A.CMAX: return -1,0
    S=enum_sigma(I,limit=1500)
    if len(S)<3 or len(S)>=1500: return -1,0
    Sg=np.array(S,dtype=np.int8)
    nb=nonbij_triples(Sg,len(I.H))
    return min(len(nb),20)+min(len(S),30)*0.01, len(nb)
A.score=lambda I: (lambda r:(r[0],r[1]))(score(I))
if __name__=='__main__':
    seed=int(sys.argv[1]); n=int(sys.argv[2]); m=int(sys.argv[3]); steps=int(sys.argv[4]); A.CMAX=int(sys.argv[5])
    t0=time.time()
    best,k=A.run(seed,n,m,steps)
    print('seed',seed,'n',n,'m',m,'c<=',A.CMAX,'nonbij boolean triples',k,'time',round(time.time()-t0),flush=True)
    if k and k>0: print(best[1],sorted(best[2]))
