# Re-check explicit examples of the paper by exhaustive partitions
import sys; sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from r1_lib import *
def E_of(s): return set(frozenset((int(a),int(b))) for a,b in (x.strip() for x in s.split(',')))
def show(n,E,H=None,label=''):
    ap=allparts_cost(n,E); m=min(c for c,_ in ap)
    opts=[p for c,p in ap if c==m]
    print(label,'n',n,'OPT',m,'#opt',len(opts))
    return m,opts
def valid(n,E,H):
    used=set();ok=True
    for t in H:
        x,y,z=t
        if not(frozenset((x,y)) in E and frozenset((y,z)) in E and frozenset((x,z)) not in E): print('  not induced P3',t); ok=False
        P=pairs3(t)
        if P&used: print('  pair overlap',t); ok=False
        used|=P
    pm=[t for t in p3s(n,E) if not (pairs3(t)&used)]
    return ok, len(pm)==0
def mkE(lst): return set(frozenset((a//10,a%10)) if isinstance(a,int) and False else frozenset(a) for a in lst)
def edges(s):
    return set(frozenset((int(t[0]),int(t[1]))) for t in s.split())
# Remark domain 7
E=edges('02 04 12 23 24'); m,o=show(5,E,label='dom2'); print(o)
E=set(frozenset(e) for e in [(0,1),(0,2),(0,3),(0,4),(1,2),(1,4),(1,5),(2,4),(2,5),(3,5),(4,5)]); m,o=show(6,E,label='dom3'); print(o, valid(6,E,[(0,3,5)]))
# Prop vertexmax
E=edges('03 13 14 24 34'); m,o=show(5,E,label='vmax'); print(o, valid(5,E,[(0,3,1)]))
# Prop c4
E=edges('01 02 13 23'); m,o=show(4,E,label='c4'); print(o, valid(4,E,[(0,1,3)]))
# Prop cex a
E=edges('03 04 12 13 14 15 23 35 36'); H=[(4,1,5),(0,3,2),(1,3,6)]; print(valid(7,E,H)); m,o=show(7,E,label='cex a'); print(o)
E=edges('02 06 14 16 17 23 25 37 45 46 47 57'); H=[(0,6,4),(3,7,4),(2,5,4),(0,2,3),(6,1,7)]; print(valid(8,E,H)); m,o=show(8,E,label='cex b'); print(o)
E=edges('02 03 04 05 12 13 15 25 35 45'); H=[(3,5,4),(0,2,1)]; print(valid(6,E,H)); m,o=show(6,E,label='cex c'); print(o)
# NFbad
E=edges('01 06 12 23 29 34 45 56 67 78 79'); H=[(0,1,2),(2,3,4),(4,5,6),(6,7,8),(2,9,7)]; print(valid(10,E,H)); m,o=show(10,E,label='nfbad2'); print(o)
E=edges('01 08 12 14 23 34 45 49 56 67 78'); H=[(0,1,2),(2,3,4),(4,5,6),(6,7,8),(1,4,9)]; print(valid(10,E,H)); m,o=show(10,E,label='nfbad3'); print(o)
# G* flip gadget? 7 vertices 0,1,3,4,5,6,7 -> relabel to 0..6 skip; majority counterexample n=6
E=edges('03 05 12 14 23 24 25 34 35 45'); H=[(0,5,4),(1,2,5),(0,3,2),(1,4,3)]; print(valid(6,E,H)); m,o=show(6,E,label='maj'); print(m,len(H))
# FAR
for mm in (3,4,5):
    n=2*mm+1; E=set(frozenset((i,i+1)) for i in range(2*mm)); E.add(frozenset((1,2*mm-1)))
    H=[(2*i,2*i+1,2*i+2) for i in range(mm)]
    print('FAR',mm,valid(n,E,H))
    ap=allparts_cost(n,E); m=min(c for c,_ in ap); opts=[p for c,p in ap if c==m]
    def sep(p):
        for b in p:
            if 1 in b and 2*mm-1 in b: return False
        return True
    print(' OPT',m,'expected',mm+1,'#opt',len(opts),'all separate chord',all(sep(p) for p in opts))
