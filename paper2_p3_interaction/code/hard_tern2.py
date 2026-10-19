import itertools,pickle
from hard_closure import cat,T,comp
c=pickle.load(open("/tmp/cat3.pkl","rb"))
B=set(cat)
full=frozenset(itertools.product(range(3),repeat=2))
for R in c:
    for a,b in itertools.permutations(range(3),2):
        # projection onto coords (a,b)
        B.add(frozenset((t[a],t[b]) for t in R))
        # slices: fix third coordinate to value v
        d=3-a-b
        for v in range(3):
            B.add(frozenset((t[a],t[b]) for t in R if t[d]==v))
B|={T(r) for r in B}; B.add(full)
S=set(B); ch=True
while ch:
    ch=False
    L=list(S)
    for R in L:
        for Q in L:
            for X in (comp(R,Q),R&Q):
                if X not in S: S.add(X); ch=True
EQ=frozenset((i,i) for i in range(3)); NE=frozenset((i,j) for i in range(3) for j in range(3) if i!=j)
print(len(S),EQ in S,NE in S)
bij=[r for r in S if len(r)==3 and len({x for x,_ in r})==3 and len({y for _,y in r})==3]
print("bijections",[sorted(b) for b in bij])
print("size>=6 nonfull:",[sorted(r) for r in S if 6<=len(r)<9][:6], sum(1 for r in S if 6<=len(r)<9))
for pi in itertools.permutations(range(3)):
    N=frozenset((x,y) for x in range(3) for y in range(3) if y!=pi[x])
    print("N_pi",pi,N in S)
