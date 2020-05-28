#!/usr/bin/env python3.5

import sys
import random
import os
import math

prefix = os.getcwd() + "/Challenges/LW/"

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def usage():
    eprint("ERROR the script expects 2 integer arguments: 'n' and 'seed'.")
    eprint("\b - 'n' is an integer corresponding to the size of the matrix: H will be of size n//2 * n")
    eprint("\b - 'seed' is an integer corresponding to the initial value of the random seed")
    eprint("This script generates an instance of the low weight codeword problem.")
    eprint("This matrix H is given in systematic form. The identity part is omitted.")
    eprint("The instance is stored in 'Challenges/LW/LW_n_seed'.")


def main(n, seed):
    random.seed(seed)
    text = ""
    text += "# n\n"
    text += str(n)
    text += "\n"
    text += "# seed\n"
    text += str(seed)
    text += "\n"
    text += "# H^transpose (each line corresponds to column of H, the identity part is omitted)\n"
    for i in range(n-n//2):
        line = ""
        for j in range(n//2):
            line += str(random.randint(0,1))
        line += "\n"
        text += line
    filename = prefix + "LW_" + str(n) + "_" + str(seed)
    file  = open(filename, "w")
    file.write(text)
    file.close()

    def binomial(n, r):
        p = 1
        for i in range(1, min(r, n - r) + 1):
            p *= n
            p //= i
            n -= 1
        return p

    s = 0
    j = 0
    k = n//2
    while True:
        s += binomial(n, j)
        if s >= pow(2, n - k):
            break
        j += 1


    print("GV bound : {}".format(j))

if __name__ == "__main__":
    # execute only if run as a script
    if len(sys.argv)!=3:
        usage()
        exit(1)
    try:
        n = int(sys.argv[1])
        seed = int(sys.argv[2])
    except:
        usage()
        exit(1)
    main(n, seed)
