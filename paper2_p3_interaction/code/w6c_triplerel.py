"""Triple constraints of the zero-excess CSP and closure under the 'M-majority' (coordinatewise majority of member states, M if all three differ).
Sigma = intersection of the triple relations R_tau, each over <= 3 distinct members; so closure of every R_tau implies closure of Sigma.
Role of a pair in its member (a,b,c): ab (edge), bc (edge), ac (non-edge); status maps: ab: L0 R1 M1; bc: L1 R0 M1; ac: L0 R0 M1.  Hard pairs are constants."""
import itertools
ROLES = ['h0', 'h1', 'ab', 'bc', 'ac']
F = {'ab': {0: 0, 1: 1, 2: 1}, 'bc': {0: 1, 1: 0, 2: 1}, 'ac': {0: 0, 1: 0, 2: 1}}   # states L=0,R=1,M=2
V = {(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)}
def maj(a, b, c):
    if a == b or a == c: return a
    if b == c: return b
    return 2
def relation(roles):
    soft = [i for i, r in enumerate(roles) if r in F]
    out = set()
    for s in itertools.product(range(3), repeat=len(soft)):
        st = []; k = 0
        for r in roles:
            if r == 'h0': st.append(0)
            elif r == 'h1': st.append(1)
            else: st.append(F[r][s[k]]); k += 1
        if tuple(st) in V: out.add(s)
    return out, len(soft)
bad = []
for roles in itertools.product(ROLES, repeat=3):
    X, k = relation(roles)
    if k == 0: continue
    Xs = set(X)
    for x1, x2, x3 in itertools.product(X, repeat=3):
        y = tuple(maj(x1[i], x2[i], x3[i]) for i in range(k))
        if y not in Xs:
            bad.append((roles, x1, x2, x3, y)); break
print('role patterns checked:', sum(1 for r in itertools.product(ROLES, repeat=3) if relation(r)[1] > 0))
print('patterns not closed under M-majority:', len(bad))
for b in bad[:10]: print(b)
