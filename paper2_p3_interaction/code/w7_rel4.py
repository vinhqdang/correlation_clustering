import pickle, itertools
compat=pickle.load(open('/tmp/w7_compat.pkl','rb'))
allA=None; allB=None; prod=True
for ty,res in compat.items():
    pa={a for a,b in res}; pb={b for a,b in res}
    isprod=len(res)==len(pa)*len(pb)
    prod&=isprod
    print(ty,len(pa),len(pb),isprod)
    allA=pa if allA is None else allA&pa
    allB=pb if allB is None else allB&pb
print('all product:',prod,'rules good as A-side for every type:',len(allA),' as B-side:',len(allB))
both=None
for ty,res in compat.items():
    pa={a for a,b in res}|set(); pb={b for a,b in res}
    s=pa&pb if False else None
# a rule that is good on A-side of all types AND B-side of all types
print('rules ok on both sides of all types:',len(allA&allB))
print()
side={}  # (pair e, side) -> intersection of allowed rules over eps and partner types
for ty,res in compat.items():
    eA,eB,eps=ty
    pa={a for a,b in res}; pb={b for a,b in res}
    for e,p in ((eA,pa),(eB,pb)):
        side[e]=side[e]&p if e in side else set(p)
for e in side: print(e,len(side[e]))
import itertools
for r in range(1,4):
    for sub in itertools.combinations(['ab','bc','ac'],r):
        inter=set.intersection(*[side[e] for e in sub])
        print(sub,len(inter))
