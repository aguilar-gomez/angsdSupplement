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
nindividuals=int(sys.argv[4])

outfile= open(outfilename,"w+")
infile = open(filename, "r")
outfile.write("\t".join(["chr","w_start","w_end","Watterson","Pairwise","TajimasD"])+"\n")

n=2*nindividuals


def calculate_Tajimas(window):
        '''
        Order in vector:  Watterson     Pairwise        thetaSingleton  thetaH  thetaL
        Returns a list with [sumW, sum_P, TajimasD]
        '''
        tW = np.sum([np.e**float(x[0]) for x in window])
        tP = np.sum([np.e**float(x[1]) for x in window])
        a1 = sum([1/x for x in range(1,n)])
        a2 = sum([1/x**2 for x in range(1,n)])
        b1 = (n+1)/(3*(n-1))
        b2 = 2*(n**2+n+3)/(9*n*(n-1))
        c1 = b1-1/a1
        c2 = b2-((n+2)/(a1*n))+(a2/a1**2)
        e1 = c1/a1
        e2 = c2/(a1**2+a2)
        S = tW*a1
        TajimasD = (tP-tW)/(np.sqrt(e1*S+e2*S*(S-1))) 

        return [tW, tP, TajimasD]


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
                outfile.write('\t'.join([chromosome,w_start,w_end,str(len(window))]+[str(i) for i in calculate_Tajimas(window)])+"\n")
                window = []
        else:
                last_pos = SNP_pos

infile.close()
outfile.close()
