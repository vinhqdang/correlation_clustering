"""annealing over unit cells: chain of N cells glued out-vertex -> in-vertex; objective = singleton list-core of the N=3 chain (real locks are singleton)."""
import sys, random, itertools, math, time
import numpy as np
from w9b_lib import *
from w9b_octtree import enum_sigma
from w9b_sigma import max_core
def build_chain(cell, N, crossbits):
    """cell: (v, members(list of triples), hard_edges(set), vin, vout, crosspairs list of (x,y)), crossbits: bits per crosspair"""
    v,mem,hard,vin,vout,cp=cell
    vid={}
    def V(j,l):
        if j>0 and l==vin: return V(j-1,vout)
        if (j,l) not in vid: vid[(j,l)]=len(vid)
        return vid[(j,l)]
    H=[];E=set()
    for j in range(N):
        for t in mem: H.append(tuple(V(j,x) for x in t))
    for (a,b,c) in H: E.add(key(a,b));E.add(key(b,c))
    for j in range(N):
        for (x,y) in hard: E.add(key(V(j,x),V(j,y)))
    for j in range(N-1):
        for (x,y),bit in zip(cp,crossbits):
            if bit: E.add(key(V(j,x),V(j+1,y)))
    return Inst(len(vid),E,H)
def score(cell,bits,CMAX):
    v,mem,hard,vin,vout,cp=cell
    cnt={}
    for t in mem:
        for x in t: cnt[x]=cnt.get(x,0)+1
    if cnt.get(vin,0)!=1 or cnt.get(vout,0)!=1 or vin==vout: return -1,None
    best=None
    res=[]
    for N in (2,3):
        I=build_chain(cell,N,bits)
        if not I.valid() or I.cmax()>CMAX: return -1,None
        S=enum_sigma(I,limit=3000)
        if len(S)<2 or len(S)>=3000: return -1,None
        Sg=np.array(S,dtype=np.int8)
        res.append(max_core(Sg,singleton_only=True,maxk=8))
    return res[1]*10+res[0]+0.001, tuple(res)
def rand_cell(rng,v,m,CMAX):
    for _ in range(500):
        mem=[];used=set();cnt=[0]*v;tries=0
        while len(mem)<m and tries<300:
            tries+=1
            t=tuple(rng.sample(range(v),3)); ps=[key(t[0],t[1]),key(t[1],t[2]),key(t[0],t[2])]
            if any(p in used for p in ps) or any(cnt[x]>=CMAX for x in t): continue
            mem.append(t); used|=set(ps)
            for x in t: cnt[x]+=1
        if len(mem)<m: continue
        priv=[x for x in range(v) if cnt[x]==1]
        if len(priv)<2: continue
        vin,vout=rng.sample(priv,2)
        soft=used
        pool=[p for p in itertools.combinations(range(v),2) if p not in soft]
        hard=set(p for p in pool if rng.random()<0.15)
        cp=[(x,y) for x in range(v) for y in range(v)]
        cp=rng.sample(cp,4)
        return (v,mem,hard,vin,vout,cp)
    return None
def run(seed,v,m,steps,CMAX):
    rng=random.Random(seed)
    for _ in range(200):
        cell=rand_cell(rng,v,m,CMAX)
        if cell is None: continue
        bits=[0]*4
        s,info=score(cell,bits,CMAX)
        if s>=0: break
    else: return None
    cur=(s,cell,bits); best=cur
    for it in range(steps):
        s0,c0,b0=cur
        v_,mem,hard,vin,vout,cp=c0
        hard=set(hard); bits=list(b0); mem=list(mem); vin_,vout_=vin,vout
        r=rng.random()
        pool=[p for p in itertools.combinations(range(v_),2)]
        soft=set()
        for t in mem:
            soft|={key(t[0],t[1]),key(t[1],t[2]),key(t[0],t[2])}
        if r<0.6:
            p=rng.choice(pool)
            if p in soft: continue
            hard^={p}
        elif r<0.75:
            bits[rng.randrange(4)]^=1
        elif r<0.85:
            cp=list(cp); cp[rng.randrange(4)]=(rng.randrange(v_),rng.randrange(v_))
        else:
            i=rng.randrange(len(mem)); t=tuple(rng.sample(range(v_),3)); mem2=mem[:i]+[t]+mem[i+1:]
            ps=set(); ok=True
            for tt in mem2:
                for p in (key(tt[0],tt[1]),key(tt[1],tt[2]),key(tt[0],tt[2])):
                    if p in ps: ok=False
                    ps.add(p)
            if not ok: continue
            mem=mem2; hard={p for p in hard if p not in ps}
        c1=(v_,mem,hard,vin_,vout_,cp)
        s,info=score(c1,bits,CMAX)
        if s<0: continue
        T=0.6*(1-it/steps)+0.02
        if s>=s0 or rng.random()<math.exp((s-s0)/T):
            cur=(s,c1,bits)
            if s>best[0]: best=cur
    return best
if __name__=='__main__':
    seed=int(sys.argv[1]); v=int(sys.argv[2]); m=int(sys.argv[3]); steps=int(sys.argv[4]); CMAX=int(sys.argv[5])
    t0=time.time(); b=run(seed,v,m,steps,CMAX)
    if b is None: print('no init'); sys.exit()
    s,cell,bits=b
    print('seed',seed,'v',v,'m',m,'c<=',CMAX,'best (k1(N=2),k1(N=3)) score',round(s,3),'time',round(time.time()-t0),flush=True)
    if s>=30: print(cell,bits)
