from sage.all import *
import itertools

"""
Cs = [
        codes.HammingCode(GF(2), 1),
        codes.HammingCode(GF(2), 2),
        codes.HammingCode(GF(2), 3),
        codes.HammingCode(GF(2), 4),
        codes.HammingCode(GF(2), 5),
        codes.GolayCode(GF(2),extended=False),
    ]
"""

Cs = [
        (3,1,1), # H2
        (7,4,1), # H3
        (15,11,1), # H4
        #(31,26,1), # H5
        (23,12,7),  # Golay
]

# Troncation sur k bits : (n'-k')/2
# H2 + H3 : -> n' = 

# Limites : 
# * n-k -> taille de SM, 
# * k -> nb de buckets

dimok = range(15,25)
#lengthok = range(35,45) # on vise n-k = 20?
elems = []

for i in range(1,10):
    for comb in itertools.combinations_with_replacement(Cs, i):
        k = sum(C[1] for C in comb)
        n = sum(C[0] for C in comb)
        r = sum(C[2] for C in comb)
        trunc = (n-k)/2 + k
        nns = n - r
        gap = nns - trunc
        if k in dimok:
            elems.append((comb, n, k, r, gap))
            print("[{}] : [{},{}]{}, gap={}".format(comb, n, k, r, gap))

for (comb, n, k, r, gap) in sorted(elems, key=lambda e : e[4]):
    print("[{}] : [{},{}]{}, gap={}".format(comb, n, k, r, gap))

