# angsdSupplement


#Tajimas D 

- Takes the output of angsd thetas calculations 
```
angsd  -doSaf 1 
realSFS saf2theta $pop.saf.idx -sfs $pop.folded.sfs -outname $pop
thetaStat print $pop.thetas.idx > $pop.thetas.persite.txt
```

Usage: Tajimas_SNPmidpoint.py $pop.thetas.persite.txt output nSites PopSize
