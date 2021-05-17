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


def get_average_at(runs, t):

    som = 0
    for run in runs:
        for ti,di in zip(run[0],run[1]):
            if ti > t :
                som += di_last
                break
            elif ti == run[0][-1]:
                som += di
                break
            di_last = di
    return som / len(runs)

pranges = [load_run(runfile) for runfile in glob.glob("runs/prange_48h/*")]
stern = [load_run(runfile) for runfile in glob.glob("runs/stern_48h/*")]
stern_3pass = [load_run(runfile) for runfile in glob.glob("runs/stern_multiwin_48h/*")]

print("Number of runs with Prange48h : {}".format(len(pranges)))
print("Number of runs with Stern48h : {}".format(len(stern)))
print("Number of runs with Stern multiwin 48h : {}".format(len(stern_3pass)))

print("Min found with Prange48h : {}".format(min(run[1][-1] for run in pranges)))
print("Min found with Stern48h : {}".format(min(run[1][-1] for run in stern)))
print("Min found with Stern multiwin 48h : {}".format(min(run[1][-1] for run in stern_3pass)))

print("Moy low prange : {}".format( sum(run[1][-1] for run in pranges)/len(pranges)))
print("Moy low stern : {}".format( sum(run[1][-1] for run in stern)/len(stern)))
print("Moy low stern multiwin : {}".format( sum(run[1][-1] for run in stern_3pass)/len(stern_3pass)))

def plot_average_weight_per_hour():

    prange_averages, stern_averages, stern_3pass_averages = [],[],[]
    for t in range(12*48):
        prange_averages.append(get_average_at(pranges, t*3600/12))
        stern_averages.append(get_average_at(stern, t*3600/12))
        stern_3pass_averages.append(get_average_at(stern_3pass, t*3600/12))

    plt.plot([x/12 for x in range(0, 12*48)], prange_averages, label='Prange')
    plt.plot([x/12 for x in range(0, 12*48)], stern_averages, label='Stern')
    plt.plot([x/12 for x in range(0, 12*48)], stern_3pass_averages, label='Stern multi-fenêtre')
    plt.rc('font', family='serif', size=40)
    plt.xlabel('Time (h)', fontsize=40)
    plt.ylabel('Hamming weight', fontsize=40)
    plt.rcParams.update({'font.size': 40})
    plt.legend()
    plt.show()

plot_average_weight_per_hour()
