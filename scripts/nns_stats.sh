#!/bin/bash

for k in $(seq 6 10); do
    for n in $(seq $(($k + 10)) $((k + 16)) ); do 
        nbvec=$(($k*2+4))
        for run in $(seq 0 20); do
            fname="scripts/runs/nns_stats/run${run}_k${k}_n${n}_double"
            echo "$fname"
            ./nns_test rand "$n" "$k" "$nbvec" > "$fname"
        done
    done
done
