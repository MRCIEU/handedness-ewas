## This script calculates polygenic scores for handedness
## based on a recent handedness GWAS (with LD adjustments using LDpred)
## and calculates DNA methylation scores for handedness
## in ALSPAC using EWAS summary statistics calculated in the NTR study
## in children and parents.
## It then tests associations between these scores and handedness.
## The results are saved to output.dir/scoring-stats.csv.

source("settings.r")

## directory for all data as well as output files
data.dir <- file.path(project.dir, "working", "data")

## directory for ALSPAC data
alspac.dir <- "~/work/alspac/data"

## directory for ALSPAC genotype data
geno.dir <- file.path(data.dir, "alspac")

## intermediate files directory 
dir.create(object.dir <- file.path(data.dir, "objects"))

## directory for output files
dir.create(output.dir <- file.path(data.dir, "output"))

## directory for files needed to calculate polygenic scores
dir.create(pgs.dir <- file.path(data.dir, "pgs-files"))

## contains coefficients for GWAS and EWAS models of
## handedness generated in NTR 
ntr.dir <- file.path(data.dir, "ntr")

## https://github.com/explodecomputer/alspac
library(alspac)
alspac::setDataDir(alspac.dir)

## https://github.com/perishky/eval.save
library(eval.save)
eval.save.dir(object.dir)

source("load-list-function.r")
## out: load.list()

source("load-aries-function.r")
## out: load.aries()

source("my-load-aries-function.r")
## out: my.load.aries()

## load/install R packages
source("load-packages-function.r")
## out: load.packages()

source("scatterplot-function.r")
## out: scatterplot()

## load microarray annotation
source("load-annotation.r")
## out: annotation

## load alspac data
source("load-alspac.r")
## out: alspac.table
 
## https://github.com/perishky/ewaff
library(ewaff)
## https://github.com/perishky/dmrff
library(dmrff)

## extract genotype data for ALSPAC
## (includes script to be run on bluecrystalp3, takes ~6 hours)
source("extract-genotypes.r")
## out: geno.dir/genotype/children.{bed,fam,bim} 
##      geno.dir/genotype/mother.{bed,fam,bim} 
##      geno.dir/genotype/snps.bim.gz


## prep gwas summary statistics file for LDpred
source("prep-gwas-summary-statistics.r") ## 10 minutes
## in: ntr.dir/data/Meta_RvL_no23andMe_noALSPAC_noNTR1.txt.sorted.tsv
##     geno.dir/genotype/snps.bim.gz
## out: pgs.dir/handedness-gwas-summary-stats-STANDARD.txt

## calculate handedness PGS in mothers and children
## (includes script to be run on bluecrystalp3, takes ~26 hours per score)
source("load-pgs.r")
## in: pgs.dir
## out: alspac.table$pgs.mom/pgs.child

alspac.table$pgs.mom <- -scale(alspac.table$pgs.mom)
alspac.table$pgs.child <- -scale(alspac.table$pgs.child)

## see how handedness relates to pgs
source("check-pgs.r")
## in: alspac.table$pgs.mom/child
##     and alspac.table$handed/handed.mom/handed.partner

library(foreign)
library(readstata13)
## Add ALSPAC ancestry PCs to alspac.table
source("load-alspac-ancestry.r")
## in: data.dir/alspac/ancestry-pca/
## out: alspac.table$ancestry.child/mom.pc1-10

## Get lists of unrelated children and mothers
source("load-unrelated.r")
## in: data.dir/alspac/unrelated/
## out: unrelated.children, unrelated.mothers
unrelated.mothers <- sub("M", "", unrelated.mothers)

## thresholds for building methylation scores
thresholds <- c(1e-1, 1e-3, 1e-5)

## load R package for calculating methylation scores
## https://github.com/perishky/meffonym
library(meffonym)

## add NTR meth score models to meffonym
## so that scores can be calculated using
## meffonym.score(methylation,modelname)
source("load-pms-models.r")
## in: ntr.dir/data/NTR_ewas_hand01_.*
## out: models "child.handedness0.1/0.001/1e-5"
##         and "adult.handedness0.1/0.001/1e-5"
##      added to meffonym

adult.time.points <- c("antenatal","FOM")#,"FOF","parents")
child.time.points <- c("cord","F7","15up","F24")
time.points <- c(child.time.points,adult.time.points)

## calculate handedness methylation scores in mothers and children
source("load-pms.r")
## in: time.points, meffonym.models(), meffonym.score(),
##     aries.dir, my.load.aries()
## out: alspac.table$[time.points]x[ms.models]


## prep datasets for PGS models
source("prep-pgs-model-datasets.r")
## in: alspac.table, time.points, child.time.points, adult.time.points
## out: pgs.dat

## extract covariate data used for EWAS models
## (will be used to construct datasets for PMS models)
source("extract-ewas-covariate-datasets.r")
## in: aries.dir (to load aries samplesheets),
##     eval.ret() (to load saved EWAS R data files)
## out: ewas.dat

## construct datasets for PGS models
## (but restricted to samples with methylation)
sub.dat <- sapply(time.points, function(time) {
    pgs.dat[[time]][match(ewas.dat[[time]]$id, pgs.dat[[time]]$id),]
}, simplify=F)

## prep datasets for PMS models
source("prep-pms-model-datasets.r")
## in: sub.dat, ewas.dat, time.points, child.time.points, adult.time.points
##     alspac.table (with PMS variables)
## out: pms.dat

## remove related individuals from the ALSPAC dataset with genetic data (pgs.dat)
for (time in c("cord","F7","15up","F24")) 
    pgs.dat[[time]] <- pgs.dat[[time]][which(pgs.dat[[time]]$id %in% unrelated.children),]

for (time in c("antenatal","FOM")) 
    pgs.dat[[time]] <- pgs.dat[[time]][which(pgs.dat[[time]]$id %in% unrelated.mothers),]

## check datasets
sapply(pgs.dat, dim)
sapply(sub.dat, dim)
sapply(pms.dat, function(dats) sapply(dats, nrow))
sapply(pms.dat, function(dats) sapply(dats, ncol))
sapply(pgs.dat, function(dat) "id" %in% colnames(dat))
sapply(sub.dat, function(dat) "id" %in% colnames(dat))
sapply(pms.dat, function(dats) sapply(dats, function(dat) "id" %in% colnames(dat)))
## > sapply(pgs.dat, dim)
##      cord   F7 15up  F24 antenatal  FOM
## [1,] 7873 7873 7873 7873      7848 7848
## [2,]   14   14   14   14        13   13
## > sapply(sub.dat, dim)
##      cord  F7 15up F24 antenatal FOM
## [1,]  688 742  641 430       792 707
## [2,]   14  14   14  14        13  13
## > sapply(pms.dat, function(dats) sapply(dats, nrow))
##                       cord  F7 15up F24 antenatal FOM
## child.handedness0.001  688 742  641 430       792 707
## child.handedness0.1    688 742  641 430       792 707
## child.handedness1e-05  688 742  641 430       792 707
## > sapply(pms.dat, function(dats) sapply(dats, ncol))
##                       cord F7 15up F24 antenatal FOM
## child.handedness0.001   45 45   44  44        43  43
## child.handedness0.1     45 45   44  44        43  43
## child.handedness1e-05   45 45   44  44        43  43
## > sapply(pgs.dat, function(dat) "id" %in% colnames(dat))
##      cord        F7      15up       F24 antenatal       FOM 
##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE 
## > sapply(sub.dat, function(dat) "id" %in% colnames(dat))
##      cord        F7      15up       F24 antenatal       FOM 
##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE 
## > sapply(pms.dat, function(dats) sapply(dats, function(dat) "id" %in% colnames(dat)))
##                       cord   F7 15up  F24 antenatal  FOM
## child.handedness0.001 TRUE TRUE TRUE TRUE      TRUE TRUE
## child.handedness0.1   TRUE TRUE TRUE TRUE      TRUE TRUE
## child.handedness1e-05 TRUE TRUE TRUE TRUE      TRUE TRUE


sapply(pgs.dat, dim)
##      cord   F7 15up  F24 antenatal  FOM
## [1,] 7873 7873 7873 7873      7848 7848
## [2,]   14   14   14   14        13   13

## function to fit logistic regression models for PGS/PMS
## dat: dataset
## var: variable(s) of interest, either "pgs" or "pms" or both
fit.handed.model <- function(dat, var, show=var) {
    dat$id <- NULL
    dat <- na.omit(dat)
    null <- dat
    null[[var]] <- NULL
    full <- glm(handed ~ ., dat, family=binomial(link="logit"))
    base <- glm(handed ~ ., null, family=binomial(link="logit"))
    n <- length(residuals(full))
    stats <- coef(summary(full))
    colnames(stats) <- c("estimate","se","z","p")
    stats <- stats[which(rownames(stats)==show),]
    c(stats,
      n=n,
      var.full=var(predict(full)),
      var.base=var(predict(base)))
}
        
## fit PGS models
stats.pgs <- t(sapply(pgs.dat, fit.handed.model, var="pgs"))
stats.pgs <- data.frame(time=rownames(stats.pgs),
                       model="handedness ~ pgs + pgs covs",
                       pms="",
                       score="pgs",
                       stats.pgs)

## fit PGS models restricted to samples with methylation
stats.sub <- t(sapply(sub.dat, fit.handed.model, var="pgs"))
stats.sub <- data.frame(time=rownames(stats.sub),
                       model="handedness ~ pgs + pgs covs",
                       pms="",
                       score="pgs",
                       stats.sub)

## fit PMS+PGS models, report PMS
stats.pms <- do.call(rbind, lapply(names(pms.dat), function(time) {
    stats <- t(sapply(pms.dat[[time]], fit.handed.model, var="pms"))
    data.frame(time=time,
               model="handedness ~ pms + pgs + pms covs + ewas covs",
               pms=names(pms.dat[[time]]),
               score="pms",
               stats)
})) 
stats.pms$pms <- sub(".handedness", " p<", stats.pms$pms)

## fit PGS+PMS models, report PGS
stats.pgss <- do.call(rbind, lapply(names(pms.dat), function(time) {
    stats <- t(sapply(pms.dat[[time]], fit.handed.model, var="pgs"))
    data.frame(time=time,
               model="handedness ~ pms + pgs + pms covs + ewas covs",
               pms=names(pms.dat[[time]]),
               score="pgs",
               stats)
})) 
stats.pgss$pms <- sub(".handedness", " p<", stats.pgss$pms)


## Calculate variance explained according to:
##  > Sang Hong Lee, Michael E Goddard, Naomi R Wray, and Peter M Visscher 
##  > A Better Coefficient of Determination for Genetic Profile Analysis
##  > Genetic Epidemiology 36 : 214–224 (2012)
stats <- rbind(stats.pgs, stats.sub, stats.pms, stats.pgss)
rownames(stats) <- NULL
stats$r2.full <- stats$var.full/(stats$var.full + pi^2/3)
stats$r2.base <- stats$var.base/(stats$var.base + pi^2/3)
stats$r2.diff <- stats$r2.full-stats$r2.base

## save results
write.csv(stats, file=file.path(output.dir, "scoring-stats.csv"),
          row.names=F)

