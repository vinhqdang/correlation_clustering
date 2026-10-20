import sys; sys.path.insert(0,'.')
from r2_lib import *
def E_(s): return set(frozenset(map(int,e)) for e in s.split())
def E2(s): return set(frozenset((int(e[0]),int(e[1]))) for e in s.split())
def info(name,n,E,H):
    cs=costs(n,E); L,prs,T=parts(n); opt=int(cs.min()); ks=np.nonzero(cs==opt)[0]
    used=set(); ok=True
    for t in H:
        x,y,z=t
        ok&= frozenset((x,y)) in E and frozenset((y,z)) in E and frozenset((x,z)) not in E
        P=pairs3(t); ok&= not (P&used); used|=P
    pm=all(pairs3(t)&used for t in p3s(n,E))
    print(name,"n",n,"packing ok",ok,"pairmax",pm,"|H|",len(H),"OPT",opt,"#opt",len(ks),[blocks(L[k]) for k in ks][:8])
info("domain2",5,E2("02 04 12 23 24"),[(1,2,3)])
info("domain3",6,E2("01 02 03 04 12 14 15 24 25 35 45"),[(0,3,5)])
info("vertexmax",5,E2("03 13 14 24 34"),[(0,3,1)])
info("c4",4,E2("01 02 13 23"),[(0,1,3)])
info("cexa",7,E2("03 04 12 13 14 15 23 35 36"),[(4,1,5),(0,3,2),(1,3,6)])
info("cexb",8,E2("02 06 14 16 17 23 25 37 45 46 47 57"),[(0,6,4),(3,7,4),(2,5,4),(0,2,3),(6,1,7)])
info("cexc",6,E2("02 03 04 05 12 13 15 25 35 45"),[(3,5,4),(0,2,1)])
info("gstar",8,E2("01 04 05 14 15 45 03 13 36 37 06 56"),[(6,3,7),(4,1,3),(5,0,3),(4,5,6),(6,0,1)])
for m in (3,4,5):
    n=2*m+1; E=set(frozenset((i,i+1)) for i in range(2*m)); E.add(frozenset((1,2*m-1)))
    H=[(2*i,2*i+1,2*i+2) for i in range(m)]
    info("FAR%d"%m,n,E,H)
info("NFbad2",10,E2("01 06 12 23 29 34 45 56 67 78 79"),[(0,1,2),(2,3,4),(4,5,6),(6,7,8),(2,9,7)])
info("NFbad3",10,E2("01 08 12 14 23 34 45 49 56 67 78"),[(0,1,2),(2,3,4),(4,5,6),(6,7,8),(1,4,9)])
# shrink-rule example
info("shrink_l1_X2",5,E2("01 13 03 04 14 34 0 4".split()[0:0] or "01 13 0 5".split()[:0] or "01 13"),[(0,1,3)])
E=E2("01 13 02 23"); info("c4again",4,E,[(0,1,3)])
# tightness l=1: vertices 0,1,3,x=2,x'=4 edges 01 13 02 04 32 34 24
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(0,4),(3,2),(3,4),(2,4)]); info("shrink_l1",5,E,[(0,1,3)])
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2)]); info("shrink_Gminus",4,E,[(0,1,3)])
E=set(frozenset(e) for e in [(0,1),(1,3)]); print("G-2 OPT", int(costs(3,E).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(0,4),(3,2),(3,4)]); print("G-x OPT",int(costs(5,E).min()), "n=5 with x'")
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2)]); print("G-x(l=1) vertices 0,1,3,2: OPT",int(costs(4,E).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4)]); print("X={2,4}", int(costs(5,E).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4),(0,5),(3,5),(2,5),(4,5)]); print("X size3 OPT",int(costs(6,E).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4)]); print("X size2 OPT",int(costs(5,E).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4),(0,5),(3,5),(2,5),(4,5)]);
E.discard(frozenset((0,5)))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2)]); print("X size1 OPT",int(costs(4,E).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4)]); print("X size2 again",int(costs(5,E).min()), "minus x:", int(costs(4,set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2)])).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4),(0,5),(3,5),(2,5),(4,5)]); 
print("X size3 OPT",int(costs(6,E).min()),"G-x (size2)", int(costs(5,set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4)])).min()))
E=set(frozenset(e) for e in [(0,1),(1,3),(0,2),(3,2),(0,4),(3,4),(2,4),(0,5),(3,5),(2,5),(4,5),(0,6),(3,6),(2,6),(4,6),(5,6)]); print("X size4 OPT",int(costs(7,E).min()))
