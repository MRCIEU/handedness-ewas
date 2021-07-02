## Within-pair discordant MZ twins analysis (adjusted model, children, buccal cell DNAm)
# Authors: Suderman M., Odintsova V., van Dongen J.

library(limma)

# Load methylation data (betas)
# Load phenotype and covariate data (df)
# Variables for basic model
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
# handdiscord (handedness discordance: discordant MZ twins = "disc handedness", concordant MZ twins = "conc handedness")


# remove cases without methylation data
OL <- intersect(rownames(df),rownames(betas))
df <- df[OL,]
beta <- beta[OL,]
table(rownames(df)==rownames(beta))

# sortby Family Number
df <- df[order(df$familynumber),]
betas <- betas[rownames(df),]
n.na <- apply(is.na(betas), 2, sum)
idx <- which(n.na < nrow(betas)-1)
betas <- betas[,idx]

# calculate residuals on the total sample
model <- model.matrix(~Epi + NK + Array_rownum + Sample_Plate, data=mydata)
fit <- lmFit(t(betas), model)
betas.res <- t(residuals(fit,t(betas)))
save(betas.res , file= paste(phenotype,date,"betaresiduals.RData",sep="_"))

# select twins from discordant pairs (handedness discordance variable)
mydata <- mydata[which(!is.na(mydata$handdiscord) & mydata$handdiscord == 'disc handedness'),]
left.data <- mydata[which(mydata$hand01 == "0"),]
right.data <- mydata[which(mydata$hand01 == "1"),]
right.idx = match(rownames(right.data), rownames(betas))
left.idx = match(rownames(left.data), rownames(betas))


# Testrun for 1 cpg to determine the number of estimates in output
#------------------------------------------------------------------------------------------------------------------------------------------------------
#mydata$CpGi <- betas[,1]
        
#------------------------------------------------------------------------------------------------------------------------------------------------------
dim(betas)
dim(betas.res)
sum(is.na(betas))
sum(is.na(betas.res))
table(apply(is.na(betas.res), 2, sum))

#------------------------------------------------------------------------------------------------------------------------------------------------------
Ncpgs <- ncol(betas)
results <- matrix(NA,nrow=Ncpgs,ncol=6)
#------------------------------------------------------------------------------------------------------------------------------------------------------
#check that family number corresponds to left and right twin in the pair
print(left.data$FamilyNumber)
print(right.data$FamilyNumber)
dim(left.idx)
dim(right.idx)
table(right.data$FamilyNumber==left.data$FamilyNumber)# should give TRUE

# Testrun for 1 cpg to determine the number of estimates in output
df$CpGi <- betas[,1]
table(apply(is.na(betas.res), 2, sum))
Ncpgs <- ncol(betas)
results <- matrix(NA,nrow=Ncpgs,ncol=6)
table(right.data$familynumber==left.data$familynumber)# should give TRUE

# Run model for all CpGs and save output
for (i in 1:Ncpgs) {
    meth.left <- betas.res[left.idx,i]
    meth.right <- betas.res[right.idx,i]
    n <- sum(!is.na(meth.left)&!is.na(meth.right))
    if (n > 2) {
      stats <- t.test(meth.left, meth.right, paired=T)
      results[i,] <- c(stats$estimate, stats$statistic, stats$conf.int, stats$p.value, n)
  }
  else results[i,ncol(results)] <- NA
}

save(results, file= "ttestresults.RData",sep="_")

# END
