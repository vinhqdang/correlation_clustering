"""Triangular-lattice strips: vertices (i,j) i in Z_L, j in [0..Nrows]; members = up-triangles; translation-invariant hard statuses by displacement."""
import sys, itertools, random, time
import numpy as np
from collections import Counter
from w9b_lib import *
from w9b_sigma import max_core
from w9b_octtree import enum_sigma

def strip(L,R,orient,status,wrap=True):
    """R = number of up-triangle rows; vertices rows 0..R.  status: dict (di,dj)->0/1 for displacement classes (dj>=0; if dj==0 then di>0 by canonical)."""
    vid={}
    def V(i,j):
        i%=L
        if (i,j) not in vid: vid[(i,j)]=len(vid)
        return vid[(i,j)]
    H=[]
    for j in range(R):
        for i in range(L):
            t=(V(i,j),V(i+1,j),V(i,j+1))
            H.append(tuple(itertools.permutations(t))[orient])
    n=len(vid)
    inv={v:k for k,v in vid.items()}
    E=set()
    soft=set()
    for (a,b,c) in H:
        E.add(key(a,b)); E.add(key(b,c)); soft|={key(a,b),key(b,c),key(a,c)}
    for u,v in itertools.combinations(range(n),2):
        if key(u,v) in soft: continue
        (i1,j1),(i2,j2)=inv[u],inv[v]
        di=(i2-i1); dj=j2-j1
        if dj<0 or (dj==0 and di<0): di,dj=-di,-dj
        # canonical modulo L: choose representative with min |di|
        di%=L
        if di>L//2: di-=L
        k=(di,dj)
        if k not in status: continue   # far displacement: default non-edge (0)
        if status[k]: E.add(key(u,v))
    return Inst(n,E,H)

def classes(L,maxd=2):
    cl=set()
    for dj in range(0,maxd+1):
        for di in range(-maxd,maxd+1):
            if dj==0 and di<=0: continue
            d=di%L
            if d>L//2: d-=L
            cl.add((d,dj))
    return sorted(cl)

if __name__=='__main__':
    L=int(sys.argv[1]); R=int(sys.argv[2]); maxd=int(sys.argv[3])
    cls=classes(L,maxd)
    print('L',L,'rows',R,'classes',len(cls),cls,flush=True)
    res=Counter(); good=[]
    t0=time.time()
    for orient in range(6):
        for bits in itertools.product((0,1),repeat=len(cls)):
            status=dict(zip(cls,bits))
            I=strip(L,R,orient,status)
            if not I.valid(): continue
            S=enum_sigma(I,limit=2000)
            if len(S)<2 or len(S)>=2000: continue
            Sg=np.array(S,dtype=np.int8)
            k1=max_core(Sg,singleton_only=True,maxk=7)
            res[k1]+=1
            if k1>=3: good.append((orient,bits,len(S),k1))
    print('dist',sorted(res.items()),'time',round(time.time()-t0))
    print('k1>=3 templates:',len(good))
    import pickle; pickle.dump(good,open('/tmp/w9b_strip_%d_%d_%d.pkl'%(L,R,maxd),'wb'))
    for g in good[:10]: print(g)
