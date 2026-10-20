import sys; 
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6c_lib import *
import itertools
def necklace(N):
    # c_0..c_N = 0..N ; z_j = N+j (j=1..N)
    n = 2*N+1; E=[]
    for j in range(1,N+1):
        E += [(j-1,j),(j-1,N+j),(j,N+j)]
    return adj_from_edges(n,E)
def max_packing(A):
    P3 = p3_list(A); import scipy.sparse as sp
    from scipy.optimize import milp, LinearConstraint, Bounds
    pairs = sorted({q for t in P3 for q in pairs_of(t)}); idx={p:i for i,p in enumerate(pairs)}
    rows=[];cols=[]
    for k,t in enumerate(P3):
        for q in pairs_of(t): rows.append(idx[q]); cols.append(k)
    M = sp.csr_matrix(([1]*len(rows),(rows,cols)),shape=(len(pairs),len(P3)))
    r = milp(-np.ones(len(P3)),constraints=LinearConstraint(M,-np.inf,1),integrality=np.ones(len(P3)),bounds=Bounds(0,1))
    sel=[P3[k] for k in range(len(P3)) if r.x[k]>0.5]
    return sel
for N in range(2,9):
    A = necklace(N); n=len(A)
    c,lab = ilp_opt(A)
    H = max_packing(A)
    print(N,n,'OPT',c,'maxpack',len(H),'gap',c-len(H))
