#!/usr/bin/env python
# tandemDuplications.csv is generated from sim.R script
out = open('tandemDuplications.txt', 'w')
with open('tandemDuplications.csv') as f:
    f.readline()
    for line in f:
        l = line.rstrip().split('\t')
        start = l[2]
        end = int(l[3])+1
        size = l[4]
        out.write('duplication\t{}\t{}\t{}\n'.format(start, size, end))
out.close()
