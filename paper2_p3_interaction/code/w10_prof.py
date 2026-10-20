"""Guest (profile) lemma check for all even m=8..60 and every group u: host h=s_{4u+2}, explicit cluster partition P_a (w10_chain.explicit_partition).
Classes A0={a<=2u-1}, A1={2u}, A2={2u+1}, A3={a>=2u+2} have pairwise distinct host clusters, and the four guest allowed-sets
 W_pad = h + P(A0)+P(A3), W_X = h + P(A0), W_Y = h + P(A3), W_all = everything   have effective forbidden sets exactly {2u,2u+1}, {a>=2u}, {a<=2u+1}, {} (restricted to existing positions)."""
import sys
from w10_chain import explicit_partition
bad=0; checked=0
for m in range(8,61,2):
    U=m//2
    parts=[explicit_partition(m,a) for a in range(m)]
    for u in range(U):
        h=('s',4*u+2)
        P={}
        for a in range(m):
            for c in parts[a]:
                if h in c: P[a]=frozenset(c-{h})
        cls={'A0':[a for a in range(m) if a<=2*u-1],'A1':[a for a in range(m) if a==2*u],'A2':[a for a in range(m) if a==2*u+1],'A3':[a for a in range(m) if a>=2*u+2]}
        # within a class clusters equal
        for nm,L in cls.items():
            assert len({P[a] for a in L})<=1,(m,u,nm)
        def eff(allowed_classes):
            Wset={h}
            for nm in allowed_classes:
                for a in cls[nm][:1]: Wset|=set(P[a])
            return {a for a in range(m) if not set(P[a])<=Wset}
        exp={'pad':({'A0','A3'},{2*u,2*u+1}), 'X':({'A0'},{a for a in range(m) if a>=2*u}), 'Y':({'A3'},{a for a in range(m) if a<=2*u+1}), 'all':({'A0','A1','A2','A3'},set())}
        for nm,(cl,ex) in exp.items():
            checked+=1
            if eff(cl)!=ex: bad+=1; print('BAD',m,u,nm,eff(cl),ex)
print('checked',checked,'bad',bad)
