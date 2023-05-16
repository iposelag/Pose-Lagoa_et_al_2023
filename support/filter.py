#!/bin/python
import os
os.chdir("../../data/ML_out/")
for mrmr in range(10,501):
    name_old = "three_class/mRMR/mrmr/output_"+str(mrmr)+".txt"
    name_new = "three_class/mRMR/mrmr/"+str(mrmr)+"_output"+".txt"
    file = open(name_old)
    content = file.readlines()
    i = 8 + mrmr
    k = i + mrmr + 1
    lines = content[i:k]
    with open(name_new, 'w') as f:
        for line in lines:
            f.write(line)
