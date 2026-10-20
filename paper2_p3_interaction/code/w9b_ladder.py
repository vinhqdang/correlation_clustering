"""Diamond ladder: A_j,B_j,C_j,D_j with D_j = A_{j+1}; annealing over orientations / glue permutation / hard statuses (non-periodic)."""
import sys, random, itertools, math, time
import numpy as np
from w9b_lib import *
from w9b_sigma import max_core
from w9b_octtree import enum_sigma

def ladder_vertices(N):
    # returns list of member vertex-sets (unordered triples) and n
    vid={}
    def V(name):
        if name not in vid: vid[name]=len(vid)
        return vid[name]
    mem=[]
    # period j has vertices p,q,r,a,b,c,d ; A_j uses (p,q,a) except for j>0 where A_j = D_{j-1} with identification
    cur=None
    sets=[]
    for j in range(N):
        if j==0:
            p,q,a=V((j,'p')),V((j,'q')),V((j,'a'))
        else:
            p,q,a=glue[j]
        r,b,c,d=V((j,'r')),V((j,'b')),V((j,'c')),V((j,'d'))
        A=(p,q,a);B=(q,r,b);C=(p,r,c);D=(b,c,d)
        sets+= [A,B,C,D]
        glue[j+1]=None
    return None
def build(N,sigmas,orients,hard):
    """sigmas[j] permutation of (b,c,d) giving (p',q',a') of next period (j=0..N-2). orients: list of permutation index per member. hard: dict pair->bit (pairs by vertex ids)"""
    vid={}
    def V(name):
        if name not in vid: vid[name]=len(vid)
        return vid[name]
    sets=[]
    pqa=None
    for j in range(N):
        if j==0: p,q,a=V((0,'p')),V((0,'q')),V((0,'a'))
        else: p,q,a=pqa
        r,b,c,d=V((j,'r')),V((j,'b')),V((j,'c')),V((j,'d'))
        sets+=[(p,q,a),(q,r,b),(p,r,c),(b,c,d)]
        if j<N-1:
            bcd=(b,c,d); s=sigmas[j]; pqa=(bcd[s[0]],bcd[s[1]],bcd[s[2]])
    # the glued member A_{j+1} is the same set as D_j: drop duplicates (A_{j+1} == D_j as vertex set)
    members=[]
    seen=set()
    for t in sets:
        k=frozenset(t)
        if k in seen: continue
        seen.add(k); members.append(t)
    n=len(vid)
    H=[]
    for i,t in enumerate(members):
        perm=list(itertools.permutations(t))[orients[i%6 if False else 0]] if False else None
    return n,members,vid
