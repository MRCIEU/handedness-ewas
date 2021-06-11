#!/bin/bash

## If run on bluecrystalp3 
##   (note requires a high memory node -- 184Gb RAM) 
##   qsub -q himem -l nodes=1:ppn=16,walltime=30:00:00 -I
##   module add languages/python-3.7.7
##   module add apps/plink-2.00

## Requires ldpred
##   pip install ldpred==1.0.10
##   (most recent version crashed with "KeyError: 'raw_snps'")
## (optional) check if ldpred works (3.5 minutes running time)
##   ~/.local/bin/ldpred-unittest

## ldpred script
LDPRED=~/.local/bin/ldpred
## summary stats file in STANDARD format (see LDpred docs)
SSF=handedness-gwas-summary-stats-STANDARD.txt
## size of random selection to estimate LD
R=2500
## GWAS size
N=733711 


## With 16 cores, required ~26 hours total running time and 184Gb RAM
mkdir output-children
bash calculate-pgs-scores.sh $LDPRED children $R $SSF $N output-children/


## With 16 cores, required ~26 hours total running time and 184Gb RAM
mkdir output-mothers
bash calculate-pgs-scores.sh $LDPRED mothers $R $SSF $N output-mothers/

