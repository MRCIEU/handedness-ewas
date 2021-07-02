# EWAS of left-handedness: adjusted model, NTR children (buccal cells)
# V.Odintsova, J.van Dongen

library(gee)

# Load DNA methylation data
# Load phenotype and covariate data

# Variabels for adjusted model
# familynumber
# sex
# age (age at DNA bloodsampling)
# gestational age
# birth weight
# maternal smoking during pregnancy (0=no, 1=yes)
# Epi (epithelium cells proportion)
# NK (natural killer cells proportion)
# Array_rownum (array rownumber)
# Sample_Plate (sample plate)
# hand01 (left-handedness variable: right-handed = 0, left-handed = 1)


mydata <- data.frame(FamilyNumber, hand01,GA, bweight, matsmok01, age,sex,Sample_Plate,Array_rownum,Epi,NK)  # add your phenotype here

rownames(mydata) <- rownames(covariates)
head(mydata)


mydata <- mydata[which(!is.na(mydata$hand01)),]
mydata <- mydata[which(!is.na(mydata$GA)&!is.na(mydata$bweight)&!is.na(mydata$matsmok01)),]
summary(mydata)
nrow(mydata)

mydata <- mydata[order(mydata$FamilyNumber),]
DNAm  <- DNAm [rownames(mydata),]

# gee model for a continuous outcome (DNA methylation level) correcting for clustering of samples by family


geemodel <- function (data)
{
gee(CpGi~hand01+sex+age+GA+bweight+matsmok01+Epi+NK+Array_rownum+Sample_Plate, data=data, id=FamilyNumber, family=gaussian, corstr="exchangeable", maxiter=100, na.action=na.omit,silent=TRUE)
}




# Testrun for 1 cpg to determine the number of estimates in output
#------------------------------------------------------------------------------------------------------------------------------------------------------
mydata$CpGi <- DNAm [,1]
coeff       <- summary(geemodel(mydata))$coefficients             
#------------------------------------------------------------------------------------------------------------------------------------------------------

# Create objects to save the output from gee
#------------------------------------------------------------------------------------------------------------------------------------------------------
Nestimates <-nrow(coeff)
Ncpgs      <-ncol(DNAm )
pval       <-matrix(NA,Ncpgs,Nestimates)
estimate   <-matrix(NA,Ncpgs,Nestimates)
RobustSE   <-matrix(NA,Ncpgs,Nestimates)
RobustZ    <-matrix(NA,Ncpgs,Nestimates)
N          <- matrix(NA,Ncpgs,1)
Ncpgs
Nestimates
#------------------------------------------------------------------------------------------------------------------------------------------------------
# create function to remove outliers from CpGi
iqr.trim<-function(methylation, iqr.limit=3) {
  quantiles <- quantile(methylation, probs = c(0.25, 0.75), na.rm = T)
  iqr <- quantiles[2] - quantiles[1]
  low <- quantiles[1]
  upper <- quantiles[2]
  methylation[which(methylation < low - iqr.limit * iqr)] <- NA
  methylation[which(methylation > upper + iqr.limit * iqr)] <- NA
  methylation
}
# Run geemodel for all CpGs and save output
#------------------------------------------------------------------------------------------------------------------------------------------------------
for (i in 1:Ncpgs)
{
    mydata$CpGi <- iqr.trim(DNAm [,i]) # trimming outliers
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
#------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------
colnames(estimate)   <- rownames(coeff)
colnames(RobustSE)   <- rownames(coeff)
colnames(RobustZ)    <- rownames(coeff)
colnames(pval)       <- rownames(coeff)
rownames(estimate)   <- colnames(DNAm )
rownames(RobustSE)   <- colnames(DNAm )
rownames(RobustZ)    <- colnames(DNAm )
rownames(pval)       <- colnames(DNAm )
rownames(N)          <- colnames(DNAm )
#------------------------------------------------------------------------------------------------------------------------------------------------------

save(estimate,pval,RobustSE,RobustZ,N, file="EWASresults.RData")


