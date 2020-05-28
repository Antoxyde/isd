from sage.all import *
import itertools

H = [
[1,1,0,1,1,1,1,1,1,0],
[0,1,0,0,1,0,1,0,0,1],
[1,0,1,1,1,0,1,1,1,0],
[0,0,1,0,1,1,0,1,0,0],
[0,0,0,1,0,0,1,1,0,1],
[1,0,1,0,1,1,0,1,1,0],
[1,0,0,0,0,1,1,0,0,0],
[0,0,0,1,1,0,0,1,1,1],
[1,1,0,1,0,1,1,0,0,0],
[1,0,0,1,0,1,1,0,0,0]
]

# Micmac to get G from H

Hm = matrix(GF(2),H)
Hm = Hm.transpose()
I = identity_matrix(GF(2),10)
H = block_matrix([I,Hm],nrows=1)
C = LinearCode(H)
G = C.generator_matrix()
print("G")
print(G)
print()


n, k = C.length(), C.dimension()

# Stern parameters
sigma = 5
p = 1

# One iteration of the algorithm
I = list(range(k))
I1 = list(range(k//2))
I2 = list(range(k//2, k))
L = sample(list(range(k)),sigma)

Z = G.matrix_from_columns(list(range(k, n))) # Redundant part
# Z1 = Z.matrix_from_rows(I1)
# Z2 = Z.matrix_from_rows(I2)

for comb1 in itertools.combinations(I1, p):
    delta1 = sum(Z[i] for i in comb1)
    for comb2 in itertools.combinations(I2, p):
        delta2 = sum(Z[i] for i in comb2)
	
	# if delta1|L = delta2|L , delta1+delta2|L will vanish
        if matrix(delta1).matrix_from_columns(L) == matrix(delta2).matrix_from_columns(L):
                v2 = delta1 + delta2
                wt = 2*p + v2.hamming_weight()
                v = [0]*k
                for i in comb1+comb2:
                         v[i] = 1

                v.extend(v2)
                v = matrix(GF(2), v)
                # Sanity checks
                print(v,wt,wt == v[0].hamming_weight(), C.parity_check_matrix() * v.transpose() == 0)
