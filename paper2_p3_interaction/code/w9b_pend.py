"""pendant guards for a given chain instance and list family; exact minimum certificate (guard-set search) and intended-certificate checks."""
import itertools, sys, time
import numpy as np
from w9b_lib import *
from w9b_octtree import enum_sigma
from w9b_chain_min import guards_of_pendant
LIST2ROLE={(0,):(0,0),(1,):(1,0),(2,):(2,1),(1,2):(0,1),(0,2):(1,1),(0,1):(2,0)}   # list -> (role index ab/bc/ac, required togetherness)
def add_pendants(I,fam,endchoice):
    """fam: list of (member, list tuple).  endchoice: tuple of 0/1 per family element choosing which end gets the edge (for value 0).
    returns new Inst and intended T (far pair list)"""
    E=set(I.E); nv=I.n; T=[]
    for (i,L),ec in zip(fam,endchoice):
        role,val=LIST2ROLE[tuple(sorted(L))]
        a,b,c=I.H[i]
        pr={0:(a,b),1:(b,c),2:(a,c)}[role]
        w=nv; nv+=1
        if val==0:
            x,y=pr if ec==0 else pr[::-1]
            E.add(key(w,x)); T+= [key(w,x),key(w,y)]
        else:
            E.add(key(w,pr[0])); E.add(key(w,pr[1])); T+=[key(w,pr[0]),key(w,pr[1])]
    return Inst(nv,E,I.H),T
def min_certificate_size(I,limit_budget=40):
    pendants=[v for v in range(I.n) if not I.mem[v]]
    S=enum_sigma(I,limit=10**6,drop_far=True)
    if not S: return None
    Sg=np.array(S,dtype=np.int8); K=len(S); full=(1<<K)-1
    def mask(i,allowed):
        m=0
        for j in np.nonzero(np.isin(Sg[:,i],list(allowed)))[0]: m|=1<<int(j)
        return m
    d={}
    for w in pendants:
        for pairs,i,allowed in guards_of_pendant(I,w):
            d[pairs]=d.get(pairs,full)&mask(i,allowed)
    items=list(d.items())
    for budget in range(1,limit_budget):
        found=[None]
        def rec(start,cur,pairs):
            if found[0]: return
            if cur==0: found[0]=set(pairs); return
            for j in range(start,len(items)):
                pj,mj=items[j]; np_=pairs|pj
                if len(np_)>budget: continue
                nxt=cur&mj
                if nxt==cur: continue
                rec(j+1,nxt,np_)
                if found[0]: return
        rec(0,full,frozenset())
        if found[0]: return budget,found[0]
    return None
