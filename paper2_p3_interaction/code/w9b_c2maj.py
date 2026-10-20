import pickle, itertools, sys
from w9b_poly import *
from w9b_scope import *
exec(open('w9b_c2lang.py').read().split("print('language size'")[0])
D=(0,1,2); perms=list(itertools.permutations(D))
good=[]
relsl=[(k,R) for k,(a,R) in rels.items() if len(R)]
for pos in itertools.product(range(3),repeat=6):     # output = argument at position pos[perm index]
    def mj(x,y,z):
        if x==y or x==z: return x
        if y==z: return y
        return (x,y,z)[pos[perms.index((x,y,z))]]
    ok=True
    for name,R in relsl:
        k=len(next(iter(R)))
        Rl=list(R)
        for t1 in Rl:
            for t2 in Rl:
                for t3 in Rl:
                    if tuple(mj(t1[i],t2[i],t3[i]) for i in range(k)) not in R: ok=False;break
                if not ok:break
            if not ok:break
        if not ok:break
    if ok: good.append(pos)
print('majorities (given as positions for the 6 all-distinct patterns) preserving the c=2 language:',len(good),'of 729')
pickle.dump(good,open('/tmp/w9b_c2_majs.pkl','wb'))
print(good[:10])
