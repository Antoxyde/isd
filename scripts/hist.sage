import tqdm
import itertools
import matplotlib.pyplot as plt
import random

F2 = GF(2)

def to_int(vec):
    v = 0
    for i,j in enumerate(vec[::-1]):
        v += int(j)*2**i
    return v

def get_hist_double(C1, C2):
    k1,n1 = C1.dimension(), C1.length()
    k2,n2 = C2.dimension(), C2.length()
    
    weights = [ 0 for _ in range(n1+n2) ]
    buckets = [ [] for _ in range(2**(k1+k2)) ]

    for _ in range(1<<15):
        ee = random.randint(0, 1 << (n1+n2))
        ee = vector(F2, [(ee >> i) & 1 for i in range(n1+n2)])
        e1,e2 = vector(F2, ee[:n1]), vector(F2, ee[n1:])
        key1 = to_int(C1.decode_to_message(e1))
        key2 = to_int(C2.decode_to_message(e2))
        key = (key1 << k2) | key2
        buckets[key].append(ee)

    for bucket in buckets:
        for e1,e2 in itertools.combinations(bucket, 2):
            weights[(e1+e2).hamming_weight()] += 1
    return weights
    
def get_hist(C):
    k,n = C.dimension(), C.length()
    weights = [ 0 for _ in range(n) ]
    buckets = [ [] for _ in range(2**k) ]

    for i in range(1, 1<<15):
        x = random.randint(0, 1 << n)
        ee = vector(GF(2), [(x >> i) & 1 for i in range(n)])
        key = to_int(C.decode_to_message(ee))
        buckets[key].append(ee)

    for bucket in buckets:
        for e1,e2 in itertools.combinations(bucket, 2):
            #print("{}+{} -> {}".format(to_int(e1), to_int(e2), (e1+e2).hamming_weight()))
            weights[(e1+e2).hamming_weight()] += 1


    return weights

C = LinearCode(MatrixSpace(GF(2), 6, 20).random_element())
h1 = get_hist(C)
h12 = get_hist_double(C, C)
m1 = sum(i*j for i,j in enumerate(h1))/sum(h1)
m12 = sum(i*j for i,j in enumerate(h12))/sum(h12)
print("\t Av. ibd C1", m1.n())
print("\t Av. ibd C1+C1", m12.n())
print()
"""
for C,weights, av in [(C1,h1, m1),("double " + str(C1),h12, m12)]:
    fig, ax = plt.subplots()
    ax.set_ylabel('Number of combinations', fontsize = 25)
    ax.set_xlabel('Weight', fontsize = 25)
  
    labels = map(str, range(len(weights)))
    plt.xticks(range(len(weights)), labels, fontsize=15)
    plt.yticks(fontsize=15)
    plt.bar(list(range(len(weights))), weights)
    plt.title("intra-bucket distance for {}".format(str(C)), fontsize=25)
    plt.show()
"""
