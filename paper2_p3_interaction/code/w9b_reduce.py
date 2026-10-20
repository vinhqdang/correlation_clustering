import pickle, itertools, sys
sys.path.insert(0,'.')
P,T=pickle.load(open('/tmp/w9b_scopes.pkl','rb'))
D=(0,1,2)
tup=list(itertools.product(D,repeat=3)); tid={t:i for i,t in enumerate(tup)}
def mask(R): 
    m=0
    for t in R: m|=1<<tid[t]
    return m
def lift_pair(R,i,j):
    m=0
    for t in tup:
        if (t[i],t[j]) in R: m|=1<<tid[t]
    return m
def lift_unary(L,i):
    m=0
    for t in tup:
        if t[i] in L: m|=1<<tid[t]
    return m
all_masks=set()
tri_masks={mask(R) for R in T}
all_masks|=tri_masks
for R in P:
    for i,j in itertools.permutations(range(3),2): all_masks.add(lift_pair(R,i,j))
for k in (1,2):
    for L in itertools.combinations(D,k):
        for i in range(3): all_masks.add(lift_unary(set(L),i))
full=(1<<27)-1
print('total lifted',len(all_masks))
irred=[]
for m in tri_masks:
    a=full
    for o in all_masks:
        if o!=m and (o&m)==m: a&=o
    if a!=m: continue
    irred.append(m)
print('triangle meet-irreducible',len(irred))
pickle.dump(irred,open('/tmp/w9b_irred.pkl','wb'))
