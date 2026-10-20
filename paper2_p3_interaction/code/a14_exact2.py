"""Audit 14: EXACT full-problem test (ALL pairs editable) of the comparator on a minimal layout.
Two chains P,Q (U=5 groups, m=10): g0 pad, g1 slot (block 1), g2 separator, g3 slot (block 2), g4 pad.  Each chain's tile groups carry lock guest + anchor with a private stub tunnel.
One OR-switch comparator at the separator g2 between P (first) and Q (second): clause (block(P) <= 1) or (block(Q) > 1), i.e. block(P) <= block(Q).
OPT <= |H|+2 expected iff exists tiles with block(P) <= block(Q).  Instance built with the notes' guest calculus (w10_lib), then validated independently here.
usage: a14_exact2.py p "T1" "T2" l   with T1,T2 subsets of {1,3}"""
import sys, os, time, ast, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w10_lib import *
from w10_exact import exact_decide

def build(p, T1, T2, U=5, sepg=2):
    m=2*U; W=World(p); tiles={'P':set(T1),'Q':set(T2)}
    ch={t:W.add_chain(t,m) for t in tiles}
    zn=lambda t,u:('z',t,u)
    for t,c in ch.items():
        for u in range(U):
            host=('s',4*u+2); F={2*u,2*u+1}
            if u in tiles[t]:
                c.perm_guest(zn(t,u),host,set()); lock_guest(c,('y',t,u),('x',t,u),zn(t,u),host,F)
                W.B.hard(('x',t,u),('stub',t,u))
            else:
                eff=c.perm_guest(zn(t,u),host,F); assert eff==F
    u=sepg; host=('s',4*u+2)
    F1={a for a in range(m) if a>=2*u}; F2={a for a in range(m) if a<=2*u+1}
    cmp_guest(ch['P'],host,F1,ch['Q'],host,F2,('c',),zn('P',u),zn('Q',u))
    return W

if __name__=='__main__':
    p=int(sys.argv[1]); T1=ast.literal_eval(sys.argv[2]); T2=ast.literal_eval(sys.argv[3]); ls=[int(x) for x in sys.argv[4].split(',')]
    W=build(p,T1,T2); I=Inst(W)
    blk=lambda g:(1 if g==1 else 2)
    exp=any(blk(a)<=blk(b) for a in T1 for b in T2)
    print('p',p,'T1',T1,'T2',T2,'expected YES(l=2):',exp,'| n',I.n,'|H|',len(I.members),'valid',I.valid,'c',I.c,'maxmatch',max(W.used.values()),flush=True)
    for l in ls:
        t0=time.time(); r,sol=exact_decide(I.n,I.E,I.members,l,timeout=40000)
        print('RESULT p',p,'T1',T1,'T2',T2,'l',l,'OPT<=|H|+l ?',r,'expected(l=2)',exp,round(time.time()-t0,1),flush=True)
