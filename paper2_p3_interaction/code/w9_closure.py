"""Finite closure of rooted-subtree profiles for the soft-triangle-free forest CSP (states = partitions of the 3 vertices of a member)."""
import itertools, sys
INF=99
CAP=int(sys.argv[1]) if len(sys.argv)>1 else 3
# states: L R M O X with y=(ab,bc,ac)
ST={'L':(0,1,0),'R':(1,0,0),'M':(1,1,1),'O':(0,0,0),'X':(0,0,1)}
XC={'L':0,'R':0,'M':0,'O':1,'X':2}
def view(s,role):
    ab,bc,ac=ST[s]
    return {'a':(ab,ac),'b':(ab,bc),'c':(bc,ac)}[role]
LISTS={'full':set(ST),'LO':{'L','O'},'OX':{'O','X'},'RO':{'R','O'},'O':{'O'},'L':{'L'},'R':{'R'},'X':{'X'},'M':{'M'},'empty':set()}
def compat(vu,vv,eps):
    for i in range(2):
        for j in range(2):
            if eps[i][j]:
                if vu[i]!=vv[j]: return False
            else:
                if vu[i] and vv[j]: return False
    return True
VIEWS=[(0,0),(0,1),(1,0),(1,1)]
EPS=[((a,b),(c,d)) for a in (0,1) for b in (0,1) for c in (0,1) for d in (0,1)]
def sat(x): return INF if x>=INF else min(x,CAP+1) if x>CAP else x
def addc(a,b):
    if a>=INF or b>=INF: return INF
    s=a+b
    return s if s<=CAP else CAP+1
def profile(role,lst,children):
    """children: list of (child_role_at_u, profile, eps) ; child attached at vertex of u with role child_role_at_u (role of u's vertex), eps 2x2 between u's other two vertices (ordered as in view) and child's"""
    prof={v:INF for v in VIEWS}
    for s in ST:
        if s not in LISTS[lst]: continue
        c=XC[s]
        for (qr,pr,eps) in children:
            vu=view(s,qr)
            best=INF
            for vv in VIEWS:
                pv=pr[VIEWS.index(vv)]
                if pv<INF and compat(vu,vv,eps):
                    best=min(best,pv)
            c=addc(c,best)
            if c>=INF: break
        if c>=INF: continue
        w=view(s,role)
        if c<prof[w]: prof[w]=c
    return tuple(prof[v] for v in VIEWS)
def norm(p):
    m=min(p)
    if m>=INF: return (INF,)*4,INF
    return tuple(x-m if x<INF else INF for x in p),m
# closure
seen={}   # normalized profile -> sample derivation
frontier=[]
roles=['a','b','c']
def others(role): return [r for r in roles if r!=role]
count=0
def add(p,der):
    global count
    n,m=norm(p)
    if n not in seen:
        seen[n]=der; frontier.append(n); return True
    return False
# normalized representation reused as an actual profile (base 0): sound since only differences matter for min-based recursion? NOT exactly because saturation; use base 0 profile.
# generate leaves
for role in roles:
    for lst in LISTS:
        add(profile(role,lst,[]),('leaf',role,lst))
changed=True
rounds=0
while frontier:
    rounds+=1
    cur=list(seen.keys())
    newf=[]
    frontier=[]
    for role in roles:
        o=others(role)
        for lst in LISTS:
            # attach 0,1,2 children at the two other vertices
            ch_opts=[None]+[(pr,eps) for pr in cur for eps in EPS]
            for c1 in ch_opts:
                for c2 in ch_opts:
                    if c1 is None and c2 is None: continue
                    # only consider pairs where each present child plausible; to limit growth, restrict c2 index >= c1? no: vertices differ
                    ch=[]
                    if c1 is not None: ch.append((o[0],c1[0],c1[1]))
                    if c2 is not None: ch.append((o[1],c2[0],c2[1]))
                    p=profile(role,lst,ch)
                    add(p,(role,lst,ch))
    print('round',rounds,'profiles',len(seen),flush=True)
    if rounds>=6: break
print(len(seen))
for k in sorted(seen): print(k)
