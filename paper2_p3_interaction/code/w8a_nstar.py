import sys, time
from w8a_lib import *
from w6a_trig2 import make
def conv(n, adj):
    return G8(n, [(a,b) for a in range(n) for b in range(a+1,n) if adj[a]>>b&1])
tri = [(0, 1), (1, 2), (2, 0)]
K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
for name,F in (('triangle',tri),('K4',K4)):
    n, adj, H, cl, y, estar, ed2 = make(F)
    g = conv(n, adj)
    H=[tuple(t) for t in H]
    print(name,'n=',n,'|H|=',len(H),'valid',valid_packing(g,H),'pairmax',pair_maximal(g,H),'tw(I)',tw_exact(Igraph(H)))
    t=time.time(); nu,Hn=nu_mod(g); print('  nu_mod =',nu,'(%.1fs)'%(time.time()-t))
    t=time.time(); opt,part=opt_ilp(g); print('  OPT(ILP)=',opt,'(%.1fs)'%(time.time()-t))
    Hn=list(Hn); print('  max packing pairmax',pair_maximal(g,Hn),'tw(I(max packing)) <=',tw_exact(Igraph(Hn)),'|I comps|',len(list(nx.connected_components(Igraph(Hn)))))
