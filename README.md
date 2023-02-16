# angsdSupplement


## Tajima' s D 

- Takes the output of angsd thetas calculations 
```
angsd  -doSaf 1 
realSFS saf2theta $pop.saf.idx -sfs $pop.folded.sfs -outname $pop
thetaStat print $pop.thetas.idx > $pop.thetas.persite.txt
```

Arguments:
- inFile: name of infile, thetaStat print output( $pop.thetas.persite.txt )
- outFile: name of output file
- nSites: number of sites per window
- nIndiv: number of individuals 

Usage: 
```
Tajimas_SNPmidpoint.py inFile outFile nSites nIndiv
```
