import glob
import matplotlib.pyplot as plt
import re

REG = re.compile('niter=(\d+), time=(\d+).(\d+), wt=(\d{3})')

def load_run(runfile):
    xi,yi = [],[]
    with open(runfile, "r") as fd:
        for line in fd.read().strip().split("\n"):
            _, time, _, wt = REG.search(line).groups()
            xi.append(int(time))
            yi.append(int(wt))
    return xi,yi

for runfile in glob.glob("runs/prange*"):
    plt.plot(*load_run(runfile))

plt.xlabel("Time in second")
plt.ylabel("Hamming weigth")
plt.title("Prange runs (24h each)")
plt.ylim(220,240)
plt.show()



