#!/usr/bin/env python3
# Author : Diana Aguilar
'''
SNP windows of feature
'''
import math 
import sys
import pandas as pd
import numpy as np
from itertools import tee


filename=sys.argv[1]
outfilename=sys.argv[2]
window_min=int(sys.argv[3])
feature=sys.argv[4]

outfile= open(outfilename,"w+")
infile = open(filename, "r")
outfile.write("\t".join(["chr","w_start","w_end","Nsites",feature])+"\n")


def calculate_mean(window):
        '''
        Order in vector:  feature
        '''
        window = np.mean([float(x[0]) for x in window])
        return [window]


start_chr = "None"
window =[]
nwindow=0
first_it, second_it = tee(infile)
next(second_it)
for line in first_it:
        nextline=next(second_it)
        column=line.split("\n")[0].split("\t")
        chromosome=column[0]
        SNP_pos=column[1]
        values=column[2:]
        ncolumn=nextline.split("\n")[0].split("\t")
        nchromosome=ncolumn[0]
        nSNP_pos=ncolumn[1]
        nvalues=ncolumn[2:]
        #print(line,nextline)
        if start_chr == chromosome:
                window.append(values)
                last_chr = chromosome
                last_pos = SNP_pos
                if len(window) == 1 & nwindow == 0:
                        w_start = SNP_pos
                        nwindow+=1
                elif len(window) == 1:
                        w_start = w_end+1

        else:
                start_chr = chromosome
                window = [values]
                w_start = SNP_pos
                nwindow=0
        if len(window) == window_min:
                if nchromosome == chromosome:
                        w_end = math.floor((int(SNP_pos)+int(nSNP_pos))/2)
                else:
                        wend=SNP_pos
                outfile.write('\t'.join([chromosome,str(w_start),str(w_end),str(len(window))]+[str(i) for i in calculate_mean(window)])+"\n")
                window = []
        else:
                last_pos = SNP_pos

infile.close()
outfile.close()
