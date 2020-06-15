#!/bin/bash

ninstances=16
fname="scripts/runs/stern_48"

case "$1" in
    "stern")
        prog="main_stern";;
    "prange")
        prog="main_prange";;
    *)
        echo "Usage: ./entrypoint <stern or prange>" && exit 1;;
esac


for i in $(seq 1 $ninstances); do
    (./${prog} > "${fname}_${i}") &
done
cd -
