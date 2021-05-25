#!/bin/bash

nbvec=25
for run in $(seq 1 20); do
    for k in $(seq 10 2 20); do
        for n in $(seq $((2*$k)) 2 $((3*$k))); do
            fname="runs/nns_stats/run${run}_k${k}_n${n}"
            echo "$fname"
            ../nns_test "$n" "$k" "$nbvec" > "$fname"
        done
    done
done





