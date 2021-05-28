#!/bin/bash

nbvec=25
for k in 18; do
    #for n in $(seq $((2*$k)) 2 $((2*$k + 6))); do
    for n in 40; do
        for run in $(seq 20 1000); do
            fname="runs/nns_stats/run${run}_k${k}_n${n}"
            echo "$fname"
            ../nns_test "$n" "$k" "$nbvec" > "$fname"
        done
    done
done





