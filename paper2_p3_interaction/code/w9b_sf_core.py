"""soft-triangle-free instances: distribution of the maximum minimal inconsistent list family (k2) of Sigma_0 (far pairs dropped)."""
import sys, random, time
import numpy as np
from collections import Counter
from w9b_lib import *
from w9b_octtree import enum_sigma
from w9b_sigma import max_core
from w9b_locality import rand_tree_instance
seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); cmax=int(sys.argv[4]); mode=sys.argv[5]
rng=random.Random(seed); d1=Counter(); d2=Counter(); t0=time.time(); tot=0
for _ in range(N):
    I=rand_tree_instance(rng,m,cmax,rng.choice([0.0,0.0,0.3]))
    if I is None or not I.valid(): continue
    tri=len(I.soft_triangles())>0
    if (mode=='sf' and tri) or (mode=='tri' and not tri): continue
    S=enum_sigma(I,limit=4000,drop_far=True)
    if len(S)<2 or len(S)>=4000: continue
    Sg=np.array(S,dtype=np.int8)
    k1=max_core(Sg,singleton_only=True,maxk=7); k2=max_core(Sg,maxk=7)
    d1[k1]+=1; d2[k2]+=1; tot+=1
print(mode,'m',m,'c<=',cmax,'instances',tot,'k1',sorted(d1.items()),'k2',sorted(d2.items()),round(time.time()-t0),'s')
