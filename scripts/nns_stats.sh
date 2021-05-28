#!/bin/bash

nbvec=25
k=18
for n in 42 44 46 48; do
	for run in $(seq 0 20); do
	    fname="scripts/runs/nns_stats/run${run}_k${k}_n${n}"
	    echo "$fname"
	    ./nns_test "$n" "$k" "$nbvec" > "$fname"
	done
done
done





