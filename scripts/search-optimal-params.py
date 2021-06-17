#!/usr/bin/python3

from math import log, comb

# Fonctionne (mêmes valeurs que papier)
def proba_success(p1, p2, l1, l2, q1, q2, w, n, k1, k2):
    k = k2 + k2
    denom = min( comb(n, w), pow(2,k) )
    #denom = comb(n, w)
    print(n-k-l1-l2, w-p1-p2-q1-q2)
    num = comb(n - k1 - k2 - l1 - l2, w - p1 - p2 - q1 - q2)

    if num < 0:
        print("{}, {}".format(n - k1 - k2 - l1 - l2, w - p1 - p2 -q1 - q2))
        raise ValueError 

    num *= comb(k1, p1)
    num *= comb(k2, p2)
    num *= comb(l1, q1)
    num *= comb(l2, q2)

    print("num = {}, denom = {}".format(log(num,2), log(denom, 2)))
    
    val = num / denom

    if val == 0.0: return 0
    return log(val, 2) / n


def proba_success2(n, k, w):

    num = k * comb(n - k, w - 1)
    denom = pow(2, k) #comb(n, w)

    val = num / denom
    return log(val, 2)
    #if val == 0.0: return 0
    #return log(val, 2) / n

# Valeurs différentes que papier
# Version asymptotique
def iteration_cost(p1, p2, l1, l2, q1, q2, w, n, k1, k2):
    first = comb(k1, p1) * comb(l1, q1) 
    second = comb(k2, p2) * comb(l2, q2)
    num = first + second + first * second * pow(2, -l1 - l2)

    if num == 0.0: return 0
    return log(num, 2) / n

#"""
## Version basique
#def L(k,p):
#    return sum( comb(i,k) for i in range(1,p + 1))
#
#def iteration_cost(p1, p2, l1, l2, q1, q2, w, n, k1, k2):
#    k = k1 + k2
#    val = 0.5 * pow(n - k, 2) * (n + k) + (l1 + l2) * (L(k1,p1) + L(k2,p2) - k1)
#    val += min(1, q1) * comb(k1, p1) * L(l1,q1) + min(1,q2) * comb(k2, p2) * L(l2, q2)
#    val += (2 * (w - p1 - p2 - q1 - q2 + 1) * (p1 + p2) * comb(k1, p1) * comb(k2, p2) * comb(l1, q1) * comb(l2, q2) ) * pow(2, -l1 - l2)
#    return log(val, 2) / n
#"""

def run(wrange, qrange, prange, lrange, n, k1, k2, ps_bound=1, mem_bound=1000000):
    bv = []
    bsp = bic = bmc = 0
    for w in wrange: 
        lowest_bitops = 10000
        for q1 in qrange: 
            for p1 in prange:
                for l1 in lrange:
                    if q1 <= l1 and (n - k1 - k1 - l1 - l1) >= 0 and (w - p1 - p1 - q1 - q1) >= 0:

                        ic = iteration_cost(p1, p1, l1, l1, q1, q1, w, n, k1, k2)
                        sp = proba_success(p1, p1, l1, l1, q1, q1, w, n, k1, k2)

                        structsize = p1 * 16 + 32  + q1*16
                        tabz = 5
                        #mc2 = structsize * tabz * (comb(k1, p1) * comb(l1, q1)) # #S + #T = 2*#S
                        mc2 = (comb(k1, p1) * comb(l1, q1))
                        mc = log(mc2, 2)
                        #mc = 2*l1 + 1
                        #print("w={},av. bitops=2^{},success prob=2^{}, iteration cost=2^{},params [q,p,l]={} memory cost : 2^{}".format(w, n*(bic - bsp), n*sp, n*ic, [q1,p1,l1], mc))

                        if sp != 0.0 and n*sp < ps_bound and mc < mem_bound :#and n*sp <= 1.5:
                            bitops = n * (ic - sp) 
                            if bitops < lowest_bitops:
                                lowest_bitops = bitops
                                bic = ic # best iteration cost
                                bsp = sp # best success prob.
                                bmc = mc
                                bv = [q1,p1,l1] # best values
        print("w={},av. bitops=2^{},success prob=2^{}, iteration cost=2^{},params [q,p,l]={} memory cost : 2^{}".format(w, n*(bic - bsp), n*bsp, n*bic, bv, bmc))

if __name__ == '__main__':
    """ 
    # 256 bits paper verification
    n,k1,k2 = 6624,2565, 2564
    ps = proba_success(8,8,47,47,1,1,117,6624,2565,2564) * n
    ic = iteration_cost(8,8,47,47,1,1,117,6624,2565,2564) * n
    print("256 bits : prob. success=2^{},iteration cost=2^{},bitops=2^{}".format(ps, ic, ic - ps))
    
    # 1000 bits paper verification
    n = 30332
    k = 22968
    w = 494
    k1 = k2 = 11484
    l1 = l2 = 140
    p1 = p2 = 27
    q1 = q2 = 0
    ps = proba_success(p1,p2,l1,l2,q1,q2,w,n,k1,k2) * n
    ic = iteration_cost(p1,p2,l1,l2,q1,q2,w,n,k1,k2) * n
    print("1000 bits : prob. success=2^{},iteration cost=2^{},bitops=2^{}".format(ps, ic, ic - ps))
    
    print("Best values 256 bits : ")
    n,k1,k2 = 6624,2565, 2564
    wrange = [117]
    qrange = range(5)
    prange = range(10)
    lrange = range(10,60)
    ps = proba_success(8,8,47,47,1,1,117,6624,2565,2564) * n
    ic = iteration_cost(8,8,47,47,1,1,117,6624,2565,2564) * n

    run(wrange,qrange,prange,lrange,n, k1, k2)
    """

    #print("Best values decodingchallenge : ")
    #n,k1,k2 = 1280,320,320
    #k = k1 + k2
    #wrange = range(140,250, 5)
    #qrange = [1] #range(10)
    #prange = [2] #range(20)
    #lrange =  [18] #range(50)
    #run(wrange,qrange,prange,lrange,n, k1, k2)
    
    n,k1,k2 = 1280,320,320
    k = k1 + k2
    ps = proba_success(2,2,9,9,0,0,215,1280,320,320) * n 
    ic = iteration_cost(2,2,9,9,0,0,215,1280,320,320) * n 
    mc = log(comb(k, 2), 2)
    print("Stern 2020 : w = 215, sp=2^{}, ic=2^{}, bitops=2^{}, mc=2^{}".format(ps, ic, ic - ps, mc))

    n,k1,k2 = 1280,320,320
    k = k1 + k2
    ps = proba_success(2, 2, 9, 9, 0, 0, 215, 1280, 320, 320) * n
    ic = iteration_cost(2,2,9,9,0,0,215,1280,320,320) * n 
    mc = log(comb(k, 2), 2)
    print("Stern 2020 : w = 215, sp=2^{}, ic=2^{}, bitops=2^{}, mc=2^{}".format(ps, ic, ic - ps, mc))
    
    """
    ic = log(k, 2)
    ps = log( k * comb(n - k, 215 - 1), 2) - (n -k)
    print("Prange bitops={}".format(ic - ps))

    print("Best values decodingchallenge : ")
    wrange = [215] #range(, 5)
    qrange = range(10)#range(2)
    prange = range(20)#range(6,8)
    lrange = range(50)#range(20,23)
    # 2^7 bits par elem dans S
    # Max 6GB / iter : 6 * 2^33
    # -> Membound a 6*2^33 / 2^7 = 6 * 2^26 = 2^28.58
    # 8
    run(wrange,qrange,prange,lrange,n, k1, k2, ps_bound=1, mem_bound=28)
    """
