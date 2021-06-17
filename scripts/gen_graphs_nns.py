import glob
import matplotlib.pyplot as plt
import re
import sys
from math import log

from dataclasses import dataclass
from typing import List

"""
Generate graphs about code candidates for a good LSH
Code candidates are in runs/nns_stats/ , and have been created with nns_test.c
"""

@dataclass
class Entry:
    n: int
    k: int
    nbvec: int
    gv: int
    minw: int
    dinter_av: float
    weight_distrib: List[int]

def load_run(runfile):
    """
    Format of a run file : 
    0 n                                  40
    1 k                                  18
    2 Number of vector in F^2_n hashed   2^25
    3 Theoretical min weight (GV bound)  7
    4 Practical min weight               5
    5 Average internal distance          10.347873
    6 Number of empty bucket             0
    7 Min elements / bucket              80
    8 Max elements / bucket              132
    9 G                                  2748380001 339A280002 F5BB8C0004 7F60200008 B6C940010 2466440020 388C80040 23AF000080 3B98680100 5D1A8C0200 A6B4800400 657B080800 BCE95C1000 838542000 6FC0A04000 AE52F08000 CF6A690000 9F04B60000
    10 Weight distrib                      518 6416 74579 498966 2793868 10885378 36594202 90744016 196123338 317623440 439023658 450499928 355444122 199853965 46700133 614078 1493 0 0
    """
    with open(runfile, "r") as fd:
        try:
            lines = fd.read().strip().split("\n")
            n = int(re.findall("\d+", lines[0])[0])
            k = int(re.findall("\d+", lines[1])[0])
            nbvec = int(lines[2][-2:])
            gv = int(re.findall("\d+", lines[3])[0])
            minw = int(re.findall("\d+", lines[4])[0])
            dinter_av = float(re.findall("\d+\.\d+", lines[5])[0])
            weight_distrib = list(map(int, re.findall("\d+", lines[10])))
        except ValueError as e:
            return None
        except IndexError as e:
            return None
    return Entry(n,k,nbvec,gv,minw,dinter_av, weight_distrib)

def minw_dinter():

    runs = [data for runfile in glob.glob("runs/nns_stats/*_k18_n40") if (data := load_run(runfile)) is not None]
    xs, ys = [],[]
    for run in runs:
        ys.append(run.minw)
        xs.append(run.dinter_av)
    plt.scatter(xs,ys)
    plt.show()

def histo_weight():
    
    runs = [data for runfile in glob.glob("runs/nns_stats/*_k18_n40") if (data := load_run(runfile)) is not None]
    weight1 = runs[0].weight_distrib
    weight2 = runs[1].weight_distrib

    labels = map(str, range(19))

    plt.xticks(range(len(weight1)), labels)
    plt.bar(range(len(weight1)), weight1, width=0.3)
    plt.bar(range(len(weight2)), weight2, width=0.3)
    plt.show()


def hist_golay():
    golay = [0, 23928832, 263217152, 464257024, 2321285120, 1378263040, 4134789120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # sage
    
    labels = map(str, range(len(golay)))
    plt.xticks(range(len(golay)), labels)
    plt.bar(range(len(golay)), golay)
    plt.title("Inter bucket distance repartition for [23,12] Golay Code over GF(2)", fontsize=25)
    plt.show()


def min_max():
    runs = [data for runfile in glob.glob("runs/nns_stats/*_k18_n40") if (data := load_run(runfile)) is not None]

    plt.xlabel('Code minimal distance', fontsize=25)
    plt.ylabel('Average inter-bucket distance', fontsize=25)

    start = 3 
    stop = 8
    colors = ["blue", "red", "green", "orange", "black", "cyan"]
    for w in range(start,stop+1):
        xs,ys = [],[]
        for run in runs:
            if run.minw == w:
                xs.append(run.minw)
                ys.append(run.dinter_av)
        plt.scatter(xs, ys, color=colors[w-start], s=50, label='d={} : {} codes'.format(w, len(xs)))
    leg = plt.legend(fontsize=20)
    plt.title("Average distance per bucket by minimal distance for random linear [40,18] codes", fontsize=25)
    plt.show()


def comp_length():
    colors = ["blue", "red", "green", "orange", "black", "cyan"]
    for i,n in enumerate([40,42,44,46]): 
        runs = [data for runfile in glob.glob("runs/nns_stats/*_k18_n{}".format(n)) if (data := load_run(runfile)) is not None]
        xs,ys = [],[]
        for run in runs:
            xs.append(run.n)
            ys.append(run.dinter_av)
        plt.scatter(xs, ys, color=colors[i], s=50, label="n={} : {} codes".format(n, len(xs)))
    plt.title("Average distance per bucket by code length for random code of dimension 18")
    plt.show()

def stats():
    runs = [data.dinter_av for runfile in glob.glob("runs/nns_stats/*_k18_n40") if (data := load_run(runfile)) is not None]
    print("nb of runs : ", len(runs))
    print("av : ", sum(runs)/len(runs))
    print("min : ", min(runs))
    print("max : ", max(runs))
          
min_max()
histo_weight()
hist_golay()
minw_dinter()
comp_length()
stats()
