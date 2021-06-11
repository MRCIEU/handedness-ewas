alspac.table <- eval.save({
    if (!dictionaryGood("current")) 
        createDictionary("Current", "current")
    if (!dictionaryGood("useful"))
        createDictionary("Useful_data", "useful")
    variables <- list(## misc covariates
                      "bestgest"=c(obj="^bestgest"),# gestational age
                      "kz021"=c(obj="^cp_"), # sex
                      "c645a"=c(obj="^c_"), # maternal education
                      "b032"=c(obj="^b_"), # parity
                      "b670"=c(obj="^b_"), # prenatal smoking
                      "months"=c(obj="^mums_age"), # mother age at birth
                      "e041"=c(obj="^e_"), # type of delivery (2m age)
                      "kz030"=c(obj="^kz_"), # birthweight
                      ## ethnicity
                      c800=list(cat1="Current"), ## mom
                      c801=list(cat1="Current"), ## partner
                      c804=list(cat1="Current"), ## child
                      ## bmi
                      "FKMS1040"=c(obj="^F24_"), ## bmi 24 years
                      "FJMR022a"=c(obj="^tf4"), ## bmi 17.5 years
                      "dw042"=c(obj="d_"), ## bmi mothers at 12 weeks gestation
                      "fm1ms111"=c(obj="FOM1_"), ## bmi mothers at 18 years postnatal
                      "ff1ms111"=c(obj="FOF1_"), ## bmi fathers at 20 years postnatal
                      ## maternal prenatal smoking
                      "a200"=c(obj="^a_"),#Number of cigarettes per at 8 weeks gestation
                      "b650"=c(obj="^b_"),#Have ever smoked
                      "b663"=c(obj="^b_"),#Smoked regularly before pregnancy,antenatal,mum
                      "b665"=c(obj="^b_"),#Smoked during first trimester,antenatal,mum
                      "b667"=c(obj="^b_"),#Smoked during second trimester,antenatal,mum
                      "b669"=c(obj="^b_"),#Smoking before pregnancy,antenatal,mum
                      "b670"=c(obj="^b_"),#Smoking during first timester,antenatal,mum
                      "b671"=c(obj="^b_"),#Smoking during second trimester,antenatal,mum
                      "c482"=c(obj="^c_"),#Smoking during third trimester,antenatal,mum
                      ## child handedness
                      "kj666"=c(obj="^kj_"), # child handedness scale 1=left,3=right
                      "kj660"=c(obj="^kj_"),
                      "kj661"=c(obj="^kj_"),
                      "kj662"=c(obj="^kj_"),
                      "kj663"=c(obj="^kj_"),
                      "kj664"=c(obj="^kj_"),
                      "kj665"=c(obj="^kj_"),
                      "f7ws013"=c(obj="^f07_"), # right=1,left=2 (writing)
                      "f8at013"=c(obj="^f08_"), # left=1,right=2 (attention session)
                      "fdcm022"=c(obj="^f10_"), # right=1,left=2 (computer task)
                      "feat013"=c(obj="^F11_"), # left=1,right=2 (attention session)
                      ## partner handedness
                      "pf2000"=c(obj="^pf_"), # partner write left=1,right=2,either=3
                      "pf2001"=c(obj="^pf_"), 
                      "pf2002"=c(obj="^pf_"), 
                      "pf2003"=c(obj="^pf_"), 
                      "pf2004"=c(obj="^pf_"), 
                      "pf2005"=c(obj="^pf_"), 
                      "pf2006"=c(obj="^pf_"), 
                      "pf2007"=c(obj="^pf_"), 
                      "pf2008"=c(obj="^pf_"), 
                      "pf2009"=c(obj="^pf_"), 
                      "pf2010"=c(obj="^pf_"),
                      ## mother handedness
                      "h121"=c(obj="^h_"), # mother handedness scale 1=left 23=right
                      "h121a"=c(obj="^h_"), # mother handedness (1=left,2=mixed,3=right)
                      "h120a"=c(obj="^h_"),
                      "h120b"=c(obj="^h_"),
                      "h120c"=c(obj="^h_"),
                      "h120d"=c(obj="^h_"),
                      "h120e"=c(obj="^h_"),
                      "h120f"=c(obj="^h_"),
                      "h120g"=c(obj="^h_"),
                      "h120h"=c(obj="^h_"),
                      "h120i"=c(obj="^h_"),
                      "h120j"=c(obj="^h_"),
                      "h120k"=c(obj="^h_"))    
    alspac.vars <- findVars(names(variables))
    alspac.vars <- alspac.vars[which(alspac.vars$name %in% names(variables)),]
    stopifnot(all(names(variables) %in% alspac.vars$name))
    alspac.vars <- do.call("filterVars", c(list(alspac.vars), variables))
    extractVars(alspac.vars)
}, "alspac-table")


set.cond <- function(x,condition,value) {
    original <- x
    x[which(condition)] <- value
    if (is.numeric(x) & length(unique(x)) > 10)
        print(quantile(x, na.rm=T))
    else
        print(table(original=original, new=x, useNA="always"))
    x
}
map <- function(x, ...) {
    key <- list(...)
    x[which(!x %in% names(key))] <- NA
    xp <- unlist(key)[match(x, names(key))]
    names(xp) <- names(x)
    print(table(original=x, mapped=xp, useNA="ifany"))
    xp
}
or.na <- function(...) {
    x <- sapply(list(...), function(x) x)
    is.na <- rowSums(!is.na(x)) == 0
    ret <- rowSums(x, na.rm=T) > 0
    ret[is.na] <- NA
    ret
}
and.na <- function(...) {
    x <- do.call(cbind, list(...))
    rowSums(x, na.rm=T) == ncol(x)
}


alspac.table$sex <- map(alspac.table$kz021, "1"="M", "2"="F")
alspac.table$maternal.education <- set.cond(alspac.table$c645a,alspac.table$c645a<0,NA) 
alspac.table$parity <- set.cond(alspac.table$b032,alspac.table$b032<0,NA)
alspac.table$parity <- set.cond(alspac.table$parity,alspac.table$b032>2,3)
alspac.table$maternal.smoking <- set.cond(alspac.table$b670, alspac.table$b670<0,NA)
alspac.table$maternal.smoking <- alspac.table$maternal.smoking > 0
alspac.table$white <- map(alspac.table$c804, "1"=T, "2"=F)
alspac.table$maternal.age <- alspac.table$months/12
alspac.table$caesarean <- map(alspac.table$e041, "1"=T, "2"=F)
alspac.table$birthweight <- set.cond(alspac.table$kz030, alspac.table$kz030<0,NA)
alspac.table$gestational.age <- set.cond(alspac.table$bestgest, alspac.table$bestgest<0,NA)

alspac.table$bmi.mother.ant <- set.cond(alspac.table$dw042,alspac.table$dw042<0,NA)
alspac.table$bmi.mother.fom <- set.cond(alspac.table$fm1ms111,alspac.table$fm1ms111<0,NA)
alspac.table$bmi.24 <- set.cond(alspac.table$FKMS1040,alspac.table$FKMS1040<0,NA)
alspac.table$bmi.17 <- set.cond(alspac.table$FJMR022a,alspac.table$FJMR022a<0,NA)
alspac.table$bmi.partner.fof <- set.cond(alspac.table$ff1ms111,alspac.table$ff1ms111<0,NA)



handed <- as.matrix(alspac.table[,grepl("^kj66[0-5]+", colnames(alspac.table))])
handed[which(handed < 1)] <- NA
handed[which(handed == 1)] <- -1
handed[which(handed == 3)] <- 0
handed[which(handed == 2)] <- 1
handed <- rowSums(handed, na.rm=T)

alspac.table$handed <- NA
alspac.table$handed[which(handed < -3)] <- "left"
alspac.table$handed[which(handed > 3)] <- "right"
alspac.table$handed[which(handed %in% -3:3)] <- "mixed"
table(alspac.table$handed)
## left mixed right 
##  756  7382  7507 


partner <- as.matrix(alspac.table[,grepl("^pf20", colnames(alspac.table))])
partner[which(partner < 1)] <- NA
partner[which(partner == 1)] <- -1
partner[which(partner == 3)] <- 0
partner[which(partner == 2)] <- 1
partner <- rowSums(partner)

alspac.table$handed.partner <- NA
alspac.table$handed.partner[which(partner < -6)] <- "left"
alspac.table$handed.partner[which(partner > 6)] <- "right"
alspac.table$handed.partner[which(partner %in% -6:6)] <- "mixed"
table(alspac.table$handed.partner)
## left right mixed
##  424  4686  325

                 
alspac.table$handed.mom <- map(alspac.table$h121a, "1"="left","2"="mixed","3"="right")

table(alspac.table$handed.mom)
# left right mixed 
#  724  8528  393

## check mom calculation is correct
mom <- as.matrix(alspac.table[,grepl("^h120", colnames(alspac.table))])
mom[which(mom < 1)] <- NA
mom[which(mom == 1)] <- -1
mom[which(mom == 3)] <- 0
mom[which(mom == 2)] <- 1
mom <- rowSums(mom)
handed.mom <- rep(NA, nrow(alspac.table))
handed.mom[which(mom < -6)] <- "left"
handed.mom[which(mom > 6)] <- "right"
handed.mom[which(mom %in% -6:6)] <- "mixed"

table(alspac.table$handed.mom, handed.mom, useNA="a")
  ##      handed.mom
  ##       left right <NA>
  ## left   701     0   23
  ## right    0  8528    0
  ## <NA>     0     0 6393



fisher.test(with(alspac.table, table(handed.mom=="left", handed.partner=="left")),
            alternative="greater")$p.value
fisher.test(with(alspac.table, table(handed.mom=="left", handed=="left")),
            alternative="greater")$p.value
fisher.test(with(alspac.table, table(handed.partner=="left", handed=="left")),
            alternative="greater")$p.value
## > fisher.test(with(alspac.table, table(handed.mom=="left", handed.partner=="left")),
## +             alternative="greater")$p.value
## [1] 0.206
## > fisher.test(with(alspac.table, table(handed.mom=="left", handed=="left")),
## +             alternative="greater")$p.value
## [1] 5.299594e-9
## > fisher.test(with(alspac.table, table(handed.partner=="left", handed=="left")),
## +             alternative="greater")$p.value
## [1] 0.00318378













