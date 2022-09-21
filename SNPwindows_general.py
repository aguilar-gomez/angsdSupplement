#!/usr/bin/env python3
# Author : Diana Aguilar
'''
SNP windows of feature
'''
import sys
import pandas as pd
import numpy as np

filename=sys.argv[1]
outfilename=sys.argv[2]
window_min=int(sys.argv[3])
feature=int(sys.argv[4])

outfile= open(outfilename,"w+")
infile = open(filename, "r")
outfile.write("\t".join(["chr","w_start","w_end","Nsites",feature])+"\n")


def calculate_mean(window):
        '''
        Order in vector:  feature
        '''
        window = np.mean(window)
        return [window]


start_chr = "None"
window =[]
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
                outfile.write('\t'.join([chromosome,w_start,w_end,str(len(window))]+[str(i) for i in calculate_mean(window)])+"\n")
                window = []
        else:
                last_pos = SNP_pos

infile.close()
outfile.close()
