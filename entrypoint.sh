#!/bin/bash

ninstances=16
fname="scripts/runs/stern_48"

#case "$1" in
#    "stern")
#        prog="main_stern";;
#    "prange")
#        prog="main_prange";;
#    *)
#        echo "Usage: ./entrypoint <stern or prange>" && exit 1;;
#esac
prog="main_stern"


for i in $(seq 1 $(( ninstances - 1 ))); do
    (./${prog} > "${fname}_${i}") &
done

sleep 10
./${prog} > "${fname}_${ninstances}"

