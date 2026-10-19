from a5_lib import *
E=[(0,1),(1,2),(1,3),(3,4),(4,5),(4,7),(4,8),(5,7),(6,7),(7,8)]
adj=mk(9,E); H=[(0,1,2),(3,4,5),(6,7,8)]
print('packing',is_packing(adj,H),'pair-max',pair_maximal(adj,H))
pk=packed(H); opt,opts=optima(adj); comp=lambda v:v//3
print('OPT',opt,'|H|',len(H),'#optima',len(opts))
for lab in opts:
    cl={}
    for v in range(9): cl.setdefault(lab[v],[]).append(v)
    U=U_set(adj,pk,lab)
    print(list(cl.values()),'U',U,'inserted between different members:',[p for p in U if not adj[p[0]][p[1]] and comp(p[0])!=comp(p[1])])
