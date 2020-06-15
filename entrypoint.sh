#!/bin/bash

ninstances=16
prog="main_stern"
fname="scripts/runs/stern_48"

for i in $(seq 1 $ninstances); do
    (./${prog} > "${fname}_${i}") &
done
cd -
