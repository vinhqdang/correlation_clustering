"""Counterexample to the 'break = node hit' abstraction (RM2): two locally infeasible members M1, M2 (each repaired by a break) are adjacent to the same member M0;
the two breaks lock M0 and force a third break.  RM2 predicts excess 2; the true excess is 3 (exact subset DP over all partitions)."""
import sys
sys.path.insert(0,'/home/user/correlation_clustering/paper2_p3_interaction/code')
from w9_lib import *
from w7_dpy import opt
n=11
E=[(0,1),(0,2),(0,4),(0,5),(0,6),(1,3),(2,8),(2,10),(3,4),(3,7),(3,8),(7,8)]
H=[(2,0,1),(4,3,1),(5,0,6),(2,8,7)]
inst=Inst(n,E,H); cons=inst.compile()
print('packing valid',inst.valid_packing(),'pair-maximal',inst.pair_maximal(),'I(H) forest',nx.is_forest(inst.Ig()),'I(H) edges',sorted(inst.Ig().edges()),'c =',inst.cmax(),'soft triangles',len(inst.soft_triangles()))
mm=len(H)
for i in (1,2):
    print('member',i,H[i],'ZE-infeasible alone:',not ze_feasible_subset(inst,[i],cons=cons),'| feasible with O:',feasible_subset_with(inst,[i],{i:[3]},cons),'with X:',feasible_subset_with(inst,[i],{i:[7]},cons))
print('RM2 value (hit both regions by one O each): 2')
allowed={0:[1,2,4],1:[3],2:[3],3:[1,2,4]}
print('M1=O and M2=O and M0,M3 zero-excess feasible:',Solver(inst,allowed=allowed).solve() is not None)
r=Solver(inst).solve(want_all=True)
print('exact CSP optimum (U=empty):',r[0]-mm,'extra; optimal assignments:',[{k:NAMES[v] for k,v in s.items()} for s in r[1]])
O,sols=opt(n,sorted(E),maxsol=1)
print('exact OPT over all partitions =',O,'= |H| + ',O-mm)
