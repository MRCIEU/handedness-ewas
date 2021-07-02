## Within-pair discordant MZ twins analysis (adjusted model, adults, blood DNAm)
# Authors: Suderman M., Odintsova V., van Dongen J.

library(limma)

# Load methylation data (betas)
# Load phenotype and covariate data (df)
# Variables for basic model
# familynumber
# sex
# age (age at DNA bloodsampling)
# smoking
# BMI
# Neut_Perc (neutrophils percentage)
# Mono_Perc (monocytes percentage)
# Eos_Perc (eosinophils percentage)
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
model <- model.matrix(~ BMI + smoking + Mono_Perc + Eos_Perc + Neut_Perc + Array_rownum + Sample_Plate, data=mydata)
fit <- lmFit(t(betas), model)
betas.res <- t(residuals(fit,t(betas)))

# select twins from discordant pairs (handedness discordance variable)
df <- df[which(!is.na(df$handdiscord)
                       & mydata$handdiscord == 'disc handedness'),]   # discordant MZ twins

left.data <- df[which(df$hand01 == 0),]
right.data <- df[which(df$hand01 == 1),]
right.idx = match(rownames(right.data), rownames(betas.res))
left.idx = match(rownames(left.data), rownames(betas.res))

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

