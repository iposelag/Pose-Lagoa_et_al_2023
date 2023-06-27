#!/bin/python
import os
os.chdir(".")
for mrmr in range(1,10):
    name_old = "output/output_"+str(mrmr)+".txt"
    name_new = "output_filtered/"+str(mrmr)+"_output"+".txt"
    file = open(name_old)
    content = file.readlines()
    i = 8 + 500
    k = i + 501
    lines = content[i:k]
    with open(name_new, 'w') as f:
        for line in lines:
            f.write(line)
