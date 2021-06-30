from math import comb

n = 1280
k = n//2
w = 220
nsec = 172800

niterps = 58000000 
proba = 2**(k-n) * sum(comb(n - k, i - 1) for i in range(144, w+1))

print(nsec * niterps * proba)
