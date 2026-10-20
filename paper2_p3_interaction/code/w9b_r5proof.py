"""The single triangle-scope relation Rk that kills every conservative NU, and its Boolean restriction R5 = {x+y+z != 1}."""
import pickle, itertools
from w9b_bool import bool_nu
Rk={(0,0,0),(0,0,2),(0,1,1),(0,2,0),(1,0,1),(1,1,0),(1,1,1),(1,2,1),(2,0,0),(2,2,0)}
sub={t for t in Rk if all(x in (0,1) for x in t)}
print('Rk restricted to {L,R}^3 =',sorted(sub))
R5={t for t in itertools.product((0,1),repeat=3) if sum(t)!=1}
print('equals R5 = {x+y+z != 1}:',sub==R5)
res=pickle.load(open('/tmp/w9b_killers2.pkl','rb'))
for d,l,R,cfg in res:
    if set(R)==Rk:
        print('realised as triangle scope relation with',d,'dropped pair(s), config',cfg); break
for r in range(3,9): print('Boolean NU arity',r,'for R5:',bool_nu([R5],r) is not None)
