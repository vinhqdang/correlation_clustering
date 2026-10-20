import sys, itertools; sys.path.insert(0,'.')
from r1_lib import *
E={frozenset((int(s[0]),int(s[1]))) for s in "03 05 12 14 23 24 25 34 35 45".split()}
H=[(0,5,4),(1,2,5),(0,3,2),(1,4,3)]; n=6
def iscl(S):
    A=lambda p:(frozenset(p) in E)!=(frozenset(p) in S)
    return not any(A((x,y)) and A((y,z)) and not A((x,z)) for x,y,z in itertools.permutations(range(n),3))
Sig=[]
for st in itertools.product('LRM',repeat=4):
    S=set()
    for t,s in zip(H,st):
        x,y,z=t; S.add(frozenset({'L':(x,y),'R':(y,z),'M':(x,z)}[s]))
    if iscl(S): Sig.append(''.join(st))
print(Sig)
def maj(a,b,c):
    if a==b or a==c: return a
    if b==c: return b
    return 'M'
print(''.join(maj(*cs) for cs in zip('LMLM','MLML','RRRR')))
