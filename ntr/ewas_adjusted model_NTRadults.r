# EWAS of left-handedness: adjusted model, NTR adults
# V.Odintsova, J.van Dongen

library(gee)

# Load DNA methylation data
# Load phenotype and covariate data

# Variabels for adjusted model
# familynumber
# sex
# age (age at DNA bloodsampling)
# Neut_Perc (neutrophils percentage)
# Mono_Perc (monocytes percentage)
# Eos_Perc (eosinophils percentage)
# Array_rownum (array rownumber)
# Sample_Plate (sample plate)
# hand01 (left-handedness variable: right-handed = 0, left-handed = 1)
# smoking
# BMI
                                                
mydata <- data.frame(familynumber, sex,age,Mono_Perc,Eos_Perc,Neut_Perc,Array_rownum,Sample_Plate,hand01,
    smoking, BMI)

mydata <- mydata[order(mydata$familynumber),]
DNAm <- DNAm[rownames(mydata),]

# GEE model for a continuous outcome (DNA methylation level) correcting for clustering of samples
#------------------------------------------------------------------------------------------------------------------------------------------------------

geemodel <- function (data)
{
gee(CpGi~hand01+sex+age+Mono_Perc+Eos_Perc+Neut_Perc+Array_rownum+Sample_Plate+smoking+BMI, data=data, id=familynumber, family=gaussian, corstr="exchangeable", maxiter=100, na.action=na.omit,silent=TRUE)
}
#------------------------------------------------------------------------------------------------------------------------------------------------------

Nestimates <-nrow(coeff)
Ncpgs      <-ncol(DNAm)
pval       <-matrix(NA,Ncpgs,Nestimates)
estimate   <-matrix(NA,Ncpgs,Nestimates)
RobustSE   <-matrix(NA,Ncpgs,Nestimates)
RobustZ    <-matrix(NA,Ncpgs,Nestimates)
N          <- matrix(NA,Ncpgs,1)
Ncpgs
Nestimates

for (i in 1:Ncpgs)
{
    mydata$CpGi <- DNAm[,i]
N[i,] <- length(which(!is.na(mydata$CpGi)))
coeff <- summary(geemodel(mydata))$coefficients
for (j in 1: Nestimates)
{
estimate[i,j] <- coeff[j,1]
RobustSE[i,j] <- coeff[j,4]
RobustZ[i,j]  <- coeff[j,5]
pval[i,j]     <- 2*pnorm(-abs(coeff[j,5]))
}
}

colnames(estimate)   <- rownames(coeff)
colnames(RobustSE)   <- rownames(coeff)
colnames(RobustZ)    <- rownames(coeff)
colnames(pval)       <- rownames(coeff)
rownames(estimate)   <- colnames(DNAm)
rownames(RobustSE)   <- colnames(DNAm)
rownames(RobustZ)    <- colnames(DNAm)
rownames(pval)       <- colnames(DNAm)
rownames(N)          <- colnames(DNAm)

save(estimate,pval,RobustSE,RobustZ,N, file= "EWASresults.RData")



