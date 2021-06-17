from sage.all import *

F = GF(2)

with open(sys.argv[1], "r") as fd:
    lines = fd.read().strip().split("\n")

n = int(lines[1])
k = n//2
seed = int(lines[3])

G = []
for i in range(k):
    v = lines[5+i]
    v += "0" * i + "1" + "0" * (k - i - 1)
    G.append(vector(F,v))

G = Matrix(F, G)
C = LinearCode(G)

G = C.generator_matrix()
weights = [0 for i in range(n)]

supp = list(range(n))
ctr = 0
niter = 1000
for i in range(niter):
    while True:
        Iset = sample(supp,k)
        Gis = G.matrix_from_columns(Iset)
        try:
            Gis_inv = Gis.inverse()
            break
        except ZeroDivisionError:
            continue

    Glw = Gis_inv * G
    for r in range(k):
        l = Glw.row(r)
        weights[vector(l).hamming_weight()] += 1

print("w, proba theorique, proba experimentale, écart")
for w in range(200, 400, 5):
    occ = weights[w]
    pt = log((k * binomial(k, w - 1)) / 2**k, 2).n()
    pe = log(occ/(niter), 2).n()
    d = abs(pt - pe)
    print("{}, {}, {}".format(w, pt,  pe, d))
