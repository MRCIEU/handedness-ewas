#!/bin/bash

## extracts genotypes for ALSPAC children and mothers along with SNP info
## 
## creates children.{bed,fam,bim} 
##         mother.{bed,fam,bim} files
##         snps.bim.gz 
## 
## to be run on bluecrystalp3 where the data is located
## (I used an interactive session qsub -l nodes=1:ppn=8,walltime=6:00:00 -I)

ALSPACDIR=/panfs/panasas01/shared/alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-06-25/data

module load apps/plink-1.90

## 1. Obtain SNP frequencies, rsid and coordinates: 'all.bim.gz'
##   These files provide, for each SNP, the rsid, chromosomal position, and the first and second allele genotypes. For example:
##    1 rs58108140 0 10583 A G
##    1 rs189107123 0 10611 G C
##    1 rs180734498 0 13302 T C
##    1 rs144762171 0 13327 C G
##   They are concatenated into a single file:
cat $ALSPACDIR/genotypes/bestguess/*.bim | gzip -c > snps.bim.gz

## 2. Make a list of genotype files (data is split by chromosome)
touch mergefile.txt
for i in $ALSPACDIR/genotypes/bestguess/data_chr*.bed; do
  BASE=`dirname $i`/`basename $i .bed`
  echo "$BASE.bed $BASE.bim $BASE.fam" >> mergefile.txt
done

## remove the files for chr01
tail -n +2 "mergefile.txt" > "mergefile.tmp" && mv "mergefile.tmp" "mergefile.txt"

## 3. Extract child genotypes
plink --bfile $ALSPACDIR/genotypes/bestguess/data_chr01 \
      --merge-list mergefile.txt \
      --keep $ALSPACDIR/derived/unrelated_ids/children_unrelated.txt \
      --remove $ALSPACDIR/derived/unrelated_ids/children_exclusion_list.txt \
      --make-bed \
      --out children

## 4. Extract mother genotypes
plink --bfile $ALSPACDIR/genotypes/bestguess/data_chr01 \
      --merge-list mergefile.txt \
      --keep $ALSPACDIR/derived/unrelated_ids/mothers_unrelated.txt \
      --remove $ALSPACDIR/derived/unrelated_ids/mothers_exclusion_list.txt \
      --make-bed \
      --out mothers
