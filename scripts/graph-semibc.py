from math import log, comb

n,k = 1280,640
GV = 144
p = 2
l = 20
m = 1

"""
prange = lambda w : sum(k * comb(n-k, i-1) for i in range(GV, w+1)) / 2**(n-k)
stern = lambda w : sum(comb(k//2, p)**2 * comb(n-k-l, i-2*p) for i in range(GV, w+1)) / 2**(n-k)
semibc = lambda w : ((l+1) * comb(k//2, p)**2 * sum(comb(n-k-l, i-2*p) for i in range(GV, w+1))) / 2**(n-k)
"""
# Esperances de nb de mots de poids w sur une itération
prange = lambda w : k * comb(n-k, w-1) / 2**(n-k)
#stern = lambda w : (comb(k//2, p)**2 * comb(n-k-l, w-2*p)) / 2**(n-k)

stern = lambda w : (comb(k//2, p)**2) * (comb(n-k-l, w-2*p) / 2**(n-k))

stern_tweaked = lambda w : (m * comb(k//2, p)**2) * (comb(n-k-l, w-2*p) / 2**(n-k))
semibc = lambda w : ( comb(k//2, p)**2 * (l * comb(n-k-l, w-2*p-1)  + comb(n-k-l, w-2*p) )) / 2**(n-k)

# Stern : l=18, p=2 basique, c-c = 
# Stern tweaked : l=18, p=2, nwin = 5, c-c = 32

pts_stern_x, pts_sbc_x = [],[]
pts_stern_y, pts_sbc_y = [],[]
for w in range(200, 250, 1):
    
    t = 1/(stern_tweaked(w) * 1750) 
    pts_stern_x.append(w)
    pts_stern_y.append(log(t, 2))
    t = 1/(semibc(w) * 158) 

    pts_sbc_x.append(w)
    pts_sbc_y.append(log(t, 2))


import matplotlib.pyplot as plt

plt.plot(pts_stern_x, pts_stern_y, color="red", label="Stern's algorithm")
plt.plot(pts_sbc_x, pts_sbc_y, color="blue", label="Semi-BC algorithm")

plt.xlabel('Weight', fontsize=25)
plt.ylabel('log2 expected running time in seconds', fontsize=25)
plt.legend(fontsize=20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.show()
