"""max size of a minimal inconsistent LIST family for relations closed under a given uniform majority m (all-distinct output rule), via Baker-Pixley:
<X> = {t : all 2-projections lie in the m-closure of the 2-projections of X}.   Exhaustive for k<=4, randomized for larger k."""
import itertools, random, sys
D=(0,1,2)
LISTS=[(0,),(1,),(2,),(0,1),(0,2),(1,2)]
perms=list(itertools.permutations(D))
def make_m(pos):
    def m(x,y,z):
        if x==y or x==z: return x
        if y==z: return y
        return (x,y,z)[pos[perms.index((x,y,z))]]
    return m
def closure2(pairs,m):
    S=set(pairs)
    while True:
        new=set()
        L=list(S)
        for a in L:
            for b in L:
                for c in L:
                    r=(m(a[0],b[0],c[0]),m(a[1],b[1],c[1]))
                    if r not in S: new.add(r)
        if not new: return S
        S|=new
def consistent(lists,tuples,m):
    k=len(lists)
    P={}
    for i in range(k):
        for j in range(i+1,k):
            P[(i,j)]=closure2({(t[i],t[j]) for t in tuples},m)
    # solve binary CSP with domains lists
    asg=[None]*k
    def rec(i):
        if i==k: return True
        for v in lists[i]:
            ok=True
            for j in range(i):
                if (asg[j],v) not in P[(j,i)]: ok=False;break
            if ok:
                asg[i]=v
                if rec(i+1): return True
        return False
    return rec(0)
def search_k(k,m,trials=None,rng=None):
    best=None
    def lists_iter():
        if trials is None: yield from itertools.product(LISTS,repeat=k)
        else:
            while True: yield tuple(rng.choice(LISTS) for _ in range(k))
    cnt=0
    for ls in lists_iter():
        cnt+=1
        if trials is not None and cnt>trials: break
        opts=[]
        for g in range(k):
            o=[]
            for i in range(k):
                o.append([v for v in D if (v not in ls[i])] if i==g else list(ls[i]))
            opts.append(o)
        if trials is None:
            gens=[list(itertools.product(*opts[g])) for g in range(k)]
            for tup in itertools.product(*gens):
                if not consistent(ls,tup,m): return ls,tup
        else:
            for _ in range(20):
                tup=tuple(tuple(rng.choice(opts[g][i]) for i in range(k)) for g in range(k))
                if not consistent(ls,tup,m): return ls,tup
    return None
if __name__=='__main__':
    pos=tuple(int(c) for c in sys.argv[1])
    m=make_m(pos)
    for k in (3,4):
        r=search_k(k,m)
        print('k',k,'core exists' if r else 'no core', r)
    rng=random.Random(1)
    for k in (5,6,7,8):
        r=search_k(k,m,trials=20000,rng=rng)
        print('k',k,'random search:', 'core found' if r else 'none', r if r else '')
