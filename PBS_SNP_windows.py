#!/usr/bin/env python3
# Author : Diana Aguilar
'''
SNP windows of PBS
'''
import sys
import pandas as pd
import numpy as np

filename=sys.argv[1]
outfilename=sys.argv[2]
window_min=int(sys.argv[3])

outfile= open(outfilename,"w+")
infile = open(filename, "r")
outfile.write("\t".join(["chr","w_start","w_end","Nsites","Fst01","Fst02","Fst12","PBS0","PBS1","PBS2"])+"\n")

def calculate_fst_pbs(window):
        '''
        Order in vector:  "A01", "B01", "A02", "B02", "A12", "B12"
        Returns a list with [Fst01, Fst02, Fst12, pbs0, pbs1, pbs2]
        '''
        sum_a01 = np.sum([float(x[0]) for x in window])
        sum_b01 = np.sum([float(x[1]) for x in window])
        sum_a02 = np.sum([float(x[2]) for x in window])
        sum_b02 = np.sum([float(x[3]) for x in window])
        sum_a12 = np.sum([float(x[4]) for x in window])
        sum_b12 = np.sum([float(x[5]) for x in window])
        Fst01=sum_a01/sum_b01
        Fst02=sum_a02/sum_b02
        Fst12=sum_a12/sum_b12
        p0 = -np.log(1-Fst01)
        p1 = -np.log(1-Fst02)
        p2 = -np.log(1-Fst12)
        pbs0 = (p0 + p1 - p2) / 2
        pbs1 = (p0 + p2 - p1) / 2
        pbs2 = (p1 + p2 - p0) / 2

        return [Fst01, Fst02, Fst12, pbs0, pbs1, pbs2]

start_chr="None"
window=[]
for line in infile:
        column=line.split("\n")[0].split("\t")
        chromosome=column[0]
        SNP_pos=column[1]
        values=column[2:]
        if start_chr == chromosome:
                window.append(values)
                last_chr = chromosome
                last_pos = SNP_pos
                if len(window) == 1:
                        w_start = SNP_pos
        else:
                start_chr = chromosome
                window = [values]
                w_start = SNP_pos
        if len(window) == window_min:
                w_end = SNP_pos
                outfile.write('\t'.join([chromosome,w_start,w_end,str(len(window))]+[str(i) for i in calculate_fst_pbs(window)])+"\n")
                window = []
        else:
                last_pos = SNP_pos

infile.close()
outfile.close()
