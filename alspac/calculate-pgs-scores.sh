#!/bin/bash 

## With 16 cores, required ~26 hours total running time and 184Gb RAM
## for ~8K sample genotypes and R=2500.

LDPRED=$1 ## ldpred script
GF=$2     ## genotype file prefix (.bed,.bim,.fam)
R=$3      ## number of random samples to estimate LD
SSF=$4    ## GWAS summary statistics file (STANDARD format)
N=$5      ## GWAS sample size
OUTPUT=$6 ## output directory

## PGS scores are saved to ${OUTPUT}*.sscore file(s).

if [ ! -f "${GF}.bed" ]; then 
    echo "Genotype data file ${GF}.bed does not exist."
fi

if [ ! -f "${LDPRED}" ]; then 
    echo "LDpred script ${LDPRED} does not exist."
fi

if [ ! -f "${SSF}" ]; then 
    echo "Summary statistics file ${SSF} does not exist."
fi

## collect memory use every 60s
## stdbuf -oL top -d 60 -b | grep -e ldpred --line-buffered > top.txt

## 1. coordinate data
## make a random selection of $R samples
if [ ! -f "${OUTPUT}random.ids" ]; then
  shuf -n ${R} ${GF}.fam | awk '{print $1}' > ${OUTPUT}random.ids
  plink2 --bfile ${GF} \
         --keep ${OUTPUT}random.ids \
         --make-bed \
         --out ${OUTPUT}data
  ## ldpred is supposed to take a list of IDs as --ilist argument
  ## to specify a subset, but it wasn't clear that this was used.
  ## To be sure, a subset dataset was created using plink.
fi

${LDPRED} coord \
  --gf=${OUTPUT}data          `## genotype file prefix (.bed,.bim,.fam)` \
  --ssf-format=STANDARD       `## summary stats file format` \
  --ssf=${SSF}                `## summary stats file` \ 
  --N=${N}                    `## samples size of GWAS` \
  --out=${OUTPUT}coord        `## output file prefix` \
  > ${OUTPUT}ldpred-coord.o \
  2> ${OUTPUT}ldpred-coord.e
## 7 hours

## 2. adjust SNP effects for LD

## calculate 'LD radius' = number SNPs/12000
RADIUS=`grep -e common ${OUTPUT}ldpred-coord.o | sed 's/[^0-9]\+//g' | awk '{print $1/12000}' | cut -d . -f 1`

${LDPRED} gibbs \
  --cf=${OUTPUT}coord           `## 'ldpred coord' output` \
  --ldr=${RADIUS}               `## number SNPs/12000` \
  --N=${N}                      `## sample size of GWAS` \
  --ldf=${OUTPUT}ldfile         `## prefix of LD file (to be generated)` \
  --f=0.5                       `## assumed fraction of causal SNPs` \
  --out=${OUTPUT}ldpred-weights `## output file prefix` \
  > ${OUTPUT}ldpred-gibbs.o \
  2> ${OUTPUT}ldpred-gibbs.e
## 19 hours and 113Gb RAM

## 3. Calculate PGS from adjusted effects
for WEIGHTS in ${OUTPUT}ldpred-weights*.txt; do
    BASE=${OUTPUT}`basename $WEIGHTS .txt`
    date
    awk '{print $3,$4,$7}' ${WEIGHTS} > ${OUTPUT}weights.txt
    plink2 \
	--bfile ${GF} `## genotype file prefix (.bed,.bim,.fam)` \
	--score ${OUTPUT}weights.txt \
	--out ${BASE} \
	> ${BASE}-plink2.o \
	2> ${BASE}-plink2.e
done
## 15 minutes

