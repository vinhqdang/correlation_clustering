"""Claim (5) hub counterexample + claim (7) table."""
import itertools
from a12_lib import *
from a12_t3 import Inst, isP3, ZE
n=11
E=[(0,1),(0,2),(0,4),(0,5),(0,6),(1,3),(2,8),(2,10),(3,4),(3,7),(3,8),(7,8)]
H=[(2,0,1),(4,3,1),(5,0,6),(2,8,7)]
Es=set(K(*e) for e in E)
print('edges',len(Es),'packing valid',all(K(a,b) in Es and K(b,c) in Es and K(a,c) not in Es for a,b,c in H),
      'pairs disjoint',len(soft_set(H))==12)
print('P3 count',len(p3s(n,Es)),'pair-maximal',all(any(p in soft_set(H) for p in prs(t)) for t in p3s(n,Es)))
print('OPT(subset DP)',opt_dp(n,Es))
# maximum packing size of modification-disjoint P3s (lower bound check): greedy/exhaustive
P=p3s(n,Es)
best=0
def rec(i,used,k):
    global best
    best=max(best,k)
    for j in range(i,len(P)):
        ps=prs(P[j])
        if not any(p in used for p in ps): rec(j+1,used|set(ps),k+1)
rec(0,frozenset(),0); print('max packing',best,'(OPT>=this)')
ins=Inst(n,Es,H)
print('I(H) edges',[(i,j) for i in range(4) for j in range(i+1,4) if set(H[i])&set(H[j])], 'soft triangle',has_soft_triangle(H))
# excess model: min cost over states with U = empty, F empty/any cross edits
import itertools
def best_model(U_max=0,F_allowed=True):
    best=10**9
    cands=[p for p in ins.hard]
    for k in range(0,U_max+1):
        for W in itertools.combinations(cands,k):
            Ew=Es^set(W)
            for xs in itertools.product(['L','R','M','O','X'],repeat=4):
                E2=set(Ew); c=k
                for t,s in zip(H,xs): E2=apply_state(E2,t,ST[s]); c+=1+(s=='O')+2*(s=='X')
                if c<best and is_cluster(n,E2): best=c
    return best
for k in range(0,3): print('min total edits using <=',k,'hard edits and states:',best_model(k))
# RM2 prediction: members individually ZE-infeasible given everything? compute per-member lists
def ze_feasible(E,H,members,ignore_far_groups=False):
    pass
# alone-feasibility: Sigma over ZE states of the whole instance, is member i forced out?
for i in range(4):
    # allow all other members any of 5 states, member i restricted to ZE
    ok=False
    for xs in itertools.product(['L','R','M','O','X'],repeat=4):
        if xs[i]=='O' or xs[i]=='X': continue
        E2=set(Es)
        for t,s in zip(H,xs): E2=apply_state(E2,t,ST[s])
        if is_cluster(n,E2): ok=True;break
    print('member',i,H[i],'can be ZE (others arbitrary states):',ok)
# all optimal state vectors with U empty
sols=[]
for xs in itertools.product(['L','R','M','O','X'],repeat=4):
    E2=set(Es)
    for t,s in zip(H,xs): E2=apply_state(E2,t,ST[s])
    if is_cluster(n,E2): sols.append((sum(1+(s=='O')+2*(s=='X') for s in xs),xs))
sols.sort(); print('cheapest state vectors w/ no hard edits:',sols[:6])
# hard-vertex view of M1, M2: unary lists
def lists_for(i):
    out=[]
    a,b,c=H[i]
    for u in range(n):
        if u in H[i]: continue
        S=frozenset(x for x in (a,b,c) if K(x,u) in Es); out.append((u,sorted(S)))
    return out
print('M1 lists',lists_for(1)); print('M2 lists',lists_for(2))
