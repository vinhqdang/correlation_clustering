"""annealing over orientations and hard statuses for a FIXED member structure (list of vertex triples)."""
import sys, random, itertools, math, time
import numpy as np
from w9b_lib import *
from w9b_sigma import max_core
from w9b_octtree import enum_sigma

def ladder(N,glue):
    """diamond ladder; glue[j] = permutation of (b,c,d) giving (p',q',a') of next period"""
    vid={}
    def V(name):
        if name not in vid: vid[name]=len(vid)
        return vid[name]
    members=[]; pqa=None
    for j in range(N):
        if j==0: p,q,a=V((0,'p')),V((0,'q')),V((0,'a'))
        else: p,q,a=pqa
        r,b,c,d=V((j,'r')),V((j,'b')),V((j,'c')),V((j,'d'))
        for t in [(p,q,a),(q,r,b),(p,r,c)]: members.append(tuple(t))
        D=(b,c,d)
        if j<N-1:
            pqa=tuple(D[g] for g in glue[j])
            # D is the next A: do not add it now
        else: members.append(D)
    return len(vid),members

def make_inst(n,members,orient,hardbits,hardlist):
    H=[]
    for t,o in zip(members,orient):
        H.append(tuple(itertools.permutations(t))[o])
    soft=set()
    for (a,b,c) in H: soft|={key(a,b),key(b,c),key(a,c)}
    E=set()
    for (a,b,c) in H: E.add(key(a,b)); E.add(key(b,c))
    for p,bit in zip(hardlist,hardbits):
        if bit and p not in soft: E.add(p)
    return Inst(n,E,H)

def anneal(n,members,steps,rng,cmax,objective='k1'):
    m=len(members)
    hardlist=[p for p in itertools.combinations(range(n),2)]
    def sc(orient,hb):
        I=make_inst(n,members,orient,hb,hardlist)
        if not I.valid() or I.cmax()>cmax: return -1,None,None
        S=enum_sigma(I,limit=4000)
        if len(S)<2 or len(S)>4000: return -1,None,None
        Sg=np.array(S,dtype=np.int8)
        k1=max_core(Sg,singleton_only=True,maxk=8)
        return k1*10+min(len(S),30)*0.01,k1,I
    # init: planted partition
    perms6=list(itertools.permutations(range(3)))
    for _ in range(20000):
        lab=[rng.randrange(max(2,n//2)) for _ in range(n)]
        orient=[]; ok=True
        for t in members:
            together=[(i,j) for i,j in itertools.combinations(range(3),2) if lab[t[i]]==lab[t[j]]]
            if not together: ok=False;break
            if len(together)==3:
                orient.append(rng.randrange(6)); continue
            (i,j)=together[0]; k=3-i-j   # vertices i,j together, k apart; centre must be i or j ; ends: other of (i,j) and k
            cands=[]
            for oi,pm in enumerate(itertools.permutations(t)):
                # pm=(a,b,c): a~b together or b~c together and the other pair apart
                a,b,c=pm
                ab=lab[a]==lab[b]; bc=lab[b]==lab[c]; ac=lab[a]==lab[c]
                if (ab and not bc and not ac) or (bc and not ab and not ac): cands.append(oi)
            if not cands: ok=False;break
            orient.append(rng.choice(cands))
        if not ok: continue
        hb=[1 if lab[u]==lab[v] else 0 for (u,v) in hardlist]
        s,k,I=sc(orient,hb)
        if s>=0: break
    else: return None
    cur=(s,orient,hb); best=(s,k,orient,hb,I)
    for it in range(steps):
        s0,o0,h0=cur
        o1=list(o0); h1=list(h0)
        if rng.random()<0.25:
            o1[rng.randrange(m)]=rng.randrange(6)
        else:
            for _ in range(rng.choice([1,1,2])):
                h1[rng.randrange(len(h1))]^=1
        s,k,I=sc(o1,h1)
        if s<0: continue
        T=0.5*(1-it/steps)+0.01
        if s>=s0 or rng.random()<math.exp((s-s0)/T):
            cur=(s,o1,h1)
            if s>best[0]: best=(s,k,o1,h1,I)
    return best

if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); steps=int(sys.argv[3]); cmax=int(sys.argv[4])
    rng=random.Random(seed)
    perms=list(itertools.permutations(range(3)))
    glue=[rng.choice(perms) for _ in range(N-1)]
    n,members=ladder(N,glue)
    t0=time.time()
    r=anneal(n,members,steps,rng,cmax)
    if r is None: print('no feasible init'); sys.exit()
    s,k,o,h,I=r
    print('ladder N',N,'members',len(members),'n',n,'glue',glue,'best k1',k,'score',round(s,2),'time',round(time.time()-t0))
    if k>=3: print(I.H,sorted(I.E))

def strip_struct(L,R):
    vid={}
    def V(i,j):
        i%=L
        if (i,j) not in vid: vid[(i,j)]=len(vid)
        return vid[(i,j)]
    members=[]
    for j in range(R):
        for i in range(L):
            members.append((V(i,j),V(i+1,j),V(i,j+1)))
    return len(vid),members
