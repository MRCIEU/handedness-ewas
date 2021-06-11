stop("do the following first")
## 1 copy the following files to bluecrystalp3
##   pgs.dir/handedness-gwas-summary-stats-STANDARD.txt
##   alspac.dir/genotype/children.{bed,bim,fam}
##   alspac.dir/genotype/mothers.{bed,bim,fam}
##   calculate*-pgs-scores.sh
## 2 run the following on bluecrystalp3
##   bash calculate-handedness-pgs-scores.sh
##    (note that each score requires 26 hours of computation,
##     so might want to run once for children, once for mothers)
## 3 copy output-mothers and output-children to pgs.dir/

mothers.pgs.file <- file.path(pgs.dir, "output-mothers/ldpred-weights_LDpred_p5.0000e-01.sscore")
children.pgs.file <- file.path(pgs.dir, "output-children/ldpred-weights_LDpred_p5.0000e-01.sscore")
stopifnot(file.exists(mothers.pgs.file))
stopifnot(file.exists(children.pgs.file))

## load PGS for mothers and children
mothers.pgs <- read.table(mothers.pgs.file, header=T, sep="\t", comment.char="", stringsAsFactors=F)
children.pgs <- read.table(children.pgs.file, header=T, sep="\t", comment.char="", stringsAsFactors=F)

## add PGS to alspac.table
idx <- match(paste0(alspac.table$aln,"M"), mothers.pgs$IID)
alspac.table$pgs.mom <- mothers.pgs$SCORE1_AVG[idx]
idx <- match(paste0(alspac.table$aln,alspac.table$qlet), children.pgs$IID)
alspac.table$pgs.child <- children.pgs$SCORE1_AVG[idx]

