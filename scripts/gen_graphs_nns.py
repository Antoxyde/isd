import glob
import matplotlib.pyplot as plt
import re
import sys

def load_run(runfile):
    """
    ...
    207 n                                  40
    208 k                                  18
    209 Number of vector in F^2_n hashed   2^25
    210 Theoretical min weight (GV bound)  7
    211 Practical min weight               6
    212 Average internal distance          10.355862
    213 Number of empty bucket             0
    214 Min elements / bucket              77
    215 Max elements / bucket              128
    216 G                                  37D9AC0001 A6514C0002 D67AAC0004 64C6480008 72ECC40010 D4EF140020 97F67C0040 3BF74C0080 81C27C0100 B13480200 874CF40400 EB22580800 13B4201000 8A3B942000 BDA8784000 3B108000 C64B710000 D5BE3E0000
    """
    with open(runfile, "r") as fd:
        try:
            lines = fd.read().strip().split("\n")
            n = int(re.findall("\d+", lines[207])[0])
            k = int(re.findall("\d+", lines[208])[0])
            nbvec = int(lines[209][-2:])
            gv = int(re.findall("\d+", lines[210])[0])
            minw = int(re.findall("\d+", lines[211])[0])
            dinter_av = float(re.findall("\d+\.\d+", lines[212])[0])
        except ValueError as e:
            return None
            

    return (n,k,nbvec,gv,minw,dinter_av)

datas = [data for runfile in glob.glob("runs/nns_stats/*") if (data := load_run(runfile)) is not None]
print(datas)
