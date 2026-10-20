"""Independent check of Prop FAR (paper) and of the cost-sharing semantics of a tunnel."""
from w7_opt import *
for m in (3,4,5,6):
    n=2*m+1
    E=[(i,i+1) for i in range(2*m)]+[(1,2*m-1)]
    opt,sols=solve(n,E)
    sep=all(not any((C>>1&1) and (C>>(2*m-1)&1) for C in s) for s in sols)
    print('FAR m=%d n=%d OPT=%d (=m+1: %s) #optima=%d, every optimum deletes chord (1,%d): %s'%(m,n,opt,opt==m+1,len(sols),2*m-1,sep))
# cost sharing: members (0,1,2),(3,4,5); pendant hard neighbours 6 (of 2), 7 (of 5); tunnel 0-3
def inst(tunnel,stubs,pend=True):
    E=[(0,1),(1,2),(3,4),(4,5)]
    if pend: E+=[(2,6),(5,7)]
    n=8
    if tunnel: E.append((0,3))
    if stubs: E+= [(0,8),(3,9)]; n=10
    return n,E
for name,(t,s,p) in {'no tunnel':(0,0,1),'private stubs':(0,1,1),'one tunnel 0-3':(1,0,1),'tunnel, pendants removed':(1,0,0)}.items():
    n,E=inst(t,s,p); print('cost-sharing',name,'OPT =',solve(n,E,enum=False)[0],' (|H|=2)')
