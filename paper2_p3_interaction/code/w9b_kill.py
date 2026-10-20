import sys,pickle
from w9b_nu2 import *
r=int(sys.argv[1]); which=sys.argv[2]
rels=load(which)
def has(rs):
    s=Search(r,rs); return s.solve(verbose=False) is not None
assert not has(rels)
keep=list(range(len(rels)))
# first coarse: try dropping in chunks
chunk=len(keep)//2
while chunk>=1:
    i=0
    while i<len(keep):
        trial=keep[:i]+keep[i+chunk:]
        if not has([rels[j] for j in trial]): keep=trial
        else: i+=chunk
    chunk//=2
print('minimal killing set size',len(keep))
for j in keep:
    k,R=rels[j]; print(k,sorted(map(tuple,R.tolist())))
