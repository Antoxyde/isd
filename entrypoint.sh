#!/bin/bash

nrun="5"
ninstances=16
fpath="scripts/runs"
algname="stern"
prog="main_stern"
uniqueness="$RANDOM" # Avoid colliding results if two instances of the script are launched at the same time

for i in $(seq 1 $(( ninstances - 1 ))); do
    (./${prog} > "${fpath}/${algname}_48h_r${nrun}_${i}_${uniqueness}") &
done

# Don't launch the last one in backgroud because the job is killed when this script exits.
./${prog} > "${fpath}/${algname}_48h_r${nrun}_${ninstances}_${uniqueness}"
