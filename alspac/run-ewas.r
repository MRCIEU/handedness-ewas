## This script performs EWAS of handedness in ALSPAC
## study children (at birth, age 7 and age 15)
## and the parents (mothers at the FOM clinic,
## fathers at the FOF clinic).

source("settings.r")

dir.create(object.dir <- file.path(project.dir, "working", "data", "objects"))
dir.create(output.dir <- file.path(project.dir, "working", "data", "output"))

options(mc.cores=10)

## https://github.com/explodecomputer/alspac
library(alspac)
alspac::setDataDir(alspac.dir)

## https://github.com/perishky/eval.save
library(eval.save)
eval.save.dir(object.dir)

source("load-list-function.r",echo=T)
## out: load.list()

source("load-aries-function.r",echo=T)
## out: load.aries()

source("my-load-aries-function.r",echo=T)
## out: my.load.aries()

## load/install R packages
source("load-packages-function.r")
## out: load.packages()

source("scatterplot-function.r")
## out: scatterplot()

## load microarray annotation
source("load-annotation.r")
## out: annotation

source("load-alspac.r")
## out: alspac.table

#library(devtools)
#install_github("perishky/ewaff")
library(ewaff)
#install_github("perishky/dmrff")
library(dmrff)

## run handedness EWAS 
source("ewas.r",echo=T)
## in: alspac.table, methylation, counts.blood
## out: ewas.ret

## run handedness DMR analyses 
source("dmr.r")
## in: methylation, ewas.ret
## out: dmr.ret, output.dir/alspac*-pre.rda, output.dir/dmr-*.csv

library(markdown)
library(knitr)
library(tableone)
source("knit-report-function.r")
source("kableone-function.r")
## out: knit.report()

options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))
options(markdown.HTML.stylesheet=file.path(getwd(), "style.css"))
options(markdown.HTML.header=file.path(getwd(), "collapsible.html"))


knit.report("report.rmd", file.path(output.dir, "report.html"))






