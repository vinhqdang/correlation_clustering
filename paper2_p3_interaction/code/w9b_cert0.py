"""l=0 certificates (csp variant): minimal far-pair sets T with F'(T) > |H| ; soft-triangle-free vs with soft triangles."""
import sys, random, time
from collections import Counter
from w9b_lib import *
from w9b_leaf import leaf_cost, min_certificate
from w9b_locality import rand_tree_instance
def run(seed,N,m,cmax,mode):
    rng=random.Random(seed); dist=Counter(); tried=0; best=None
    t0=time.time()
    for _ in range(N):
        I=rand_tree_instance(rng,m,cmax,rng.choice([0.0,0.0,0.3]))
        if I is None or not I.valid(): continue
        tri=len(I.soft_triangles())>0
        if (mode=='sf' and tri) or (mode=='tri' and not tri): continue
        c0,_=leaf_cost(I,T=[])
        if c0!=len(I.H): continue
        cN,_=leaf_cost(I)
        tried+=1
        if cN==len(I.H): dist[0]+=1; continue
        T=min_certificate(I,len(I.H))
        dist[len(T)]+=1
        if best is None or len(T)>best[0]: best=(len(T),I.n,sorted(I.E),I.H,T)
        if time.time()-t0>float(sys.argv[5]): break
    return tried,sorted(dist.items()),best
if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); m=int(sys.argv[3]); cmax=int(sys.argv[4])
    for mode in ('sf','tri'):
        tried,d,best=run(seed,N,m,cmax,mode)
        print(mode,'m',m,'c<=',cmax,'tight instances tried',tried,'certificate size distribution',d,flush=True)
        if best and best[0]>=4: print('  largest',best)
