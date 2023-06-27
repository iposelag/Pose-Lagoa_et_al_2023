#!/bin/bash
#for n_features in {10..500}
#do
./mrmr -i ../../../data/mrmr/copd_ctrl_train.csv -n 500 -t 1 -v 16235 -s 225 > ../../../data/mrmr/mRMR/output.txt
#done
