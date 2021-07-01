from math import log, comb

n,k = 1280,640
GV = 144
p = 2
l = 20
m = 5

"""
prange = lambda w : sum(k * comb(n-k, i-1) for i in range(GV, w+1)) / 2**(n-k)
stern = lambda w : sum(comb(k//2, p)**2 * comb(n-k-l, i-2*p) for i in range(GV, w+1)) / 2**(n-k)
semibc = lambda w : ((l+1) * comb(k//2, p)**2 * sum(comb(n-k-l, i-2*p) for i in range(GV, w+1))) / 2**(n-k)
"""

prange = lambda w : k * comb(n-k, w-1) / 2**(n-k)
#stern = lambda w : (comb(k//2, p)**2 * comb(n-k-l, w-2*p)) / 2**(n-k)
stern = lambda w : (comb(k//2, p)**2) * 2**(-l) * (comb(n-k-l, w-2*p) / 2**(n-k-l))
stern_tweaked = lambda w : (m * comb(k//2, p)**2*2**(-l) * comb(n-k-l, w-2*p)) / 2**(n-k-l)
semibc = lambda w : ( comb(k//2, p)**2*2**(-l) * (l * comb(n-k-l, w-2*p-1)  + comb(n-k-l, w-2*p) )) / 2**(n-k-l)

# Stern : l=18, p=2 basique, c-c = 
# Stern tweaked : l=18, p=2, nwin = 5, c-c = 32
for w in range(215, 230, 1):
    print("w = {}".format(w))
    print("- Prange : proba=2^{}, 48h expectation=2^{}".format(log(prange(w), 2), log(prange(w) * 172_800 * 182000 , 2)))
    print("- Stern : proba=2^{},  48h expectation=2^{}".format(log(stern(w), 2), log(stern(w) * 172_800 * 990, 2)))
    print("- Stern tweaked : proba=2^{}, 48h expectation=2^{}".format(log(stern_tweaked(w), 2), log(stern_tweaked(w) * 172_800 * 500, 2)))
    print("- SemiBC : proba=2^{}, 48h expectation=2^{}".format(log(semibc(w), 2), log(semibc(w) * 172_800 * 200, 2)))
