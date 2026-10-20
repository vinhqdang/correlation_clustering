import sys,time
from w9b_poly import *
r=int(sys.argv[1]); name=sys.argv[2]
rels=language(name)
print('language',name,len(rels),'relations; arity counts',{k:sum(1 for v in rels.values() if v[0]==k) for k in (1,2,3)},flush=True)
t0=time.time()
ps=PolySearch(r,rels,nu=True,conservative=True)
res=ps.solve(verbose=False)
print('NU arity',r,name,'->','NONE' if res is None else ('timeout' if res=='timeout' else 'FOUND'),'cegar clauses',ps.added,'t',round(time.time()-t0,1),flush=True)
