import glob
import matplotlib.pyplot as plt
import re

REG = re.compile('niter=(\d+), time=(\d+).(\d+), wt=(\d{3})')

def load_run(runfile):
    xi,yi = [],[]
    with open(runfile, "r") as fd:
        for line in fd.read().strip().split("\n"):
            if not any([line.startswith('#'), line.startswith('0'), line.startswith('1')]):
                _, time, _, wt = REG.search(line).groups()
                xi.append(int(time))
                yi.append(int(wt))
    return xi,yi

for runfile in glob.glob("runs/prange_24_*"):
    plt.plot(*load_run(runfile))

plt.xlabel("Time in second")
plt.ylabel("Hamming weigth")
plt.title("Prange runs (24h each)")
plt.ylim(220,260)
plt.show()

for runfile in glob.glob("runs/prange_48_*"):
    plt.plot(*load_run(runfile))

plt.xlabel("Time in second")
plt.ylabel("Hamming weigth")
plt.title("Prange runs (48h each)")
plt.ylim(220,260)
plt.show()


