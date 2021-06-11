# GWAS follow-up: left-handedness (LH) / T2D
# V.Odintsova, J.van Dongen

# load GWAS summary statistics (GWAShand)
# load meta-analysis EWAS summary statistics (EWAShand)
# CHR = chromosome
# POS = position
# PVAL = p-value

# include fraction of SNPs at p <5x10-08
GWAShand <- GWAShand[which(GWAShand$PVAL < 5e-08),]

# detect CpGs located in 1Mb window from SNPs
w1Mb   <- 1000000
samechrom <- matrix(0,nrow=nrow(EWAShand),ncol=nrow(GWAShand))
distance <-  matrix(0,nrow=nrow(EWAShand),ncol=nrow(GWAShand))

for (j in 1:nrow(EWAShand)) {
  samechrom[j,which(GWAShand$CHR==EWAShand$CHR[j])] <- 1
  distance[j,] <- abs(EWAShand$POS[j] - GWAShand$POS)
}

withinw1Mb <-  matrix(0,nrow=nrow(EWAShand),ncol=nrow(GWAShand))
for (j in 1:nrow(EWAShand))
{
withinw1Mb[j,which(samechrom[j,]==1 & distance[j,] < w1Mb)] <- 1
}

CpGs_w1Mb <- rownames(EWAShand[which(rowSums(withinw1Mb) > 0),])

# label CpGs located in 1Mb window and other CpGs in meta-analysis EWAS summary statistics
pathway <- rep("other",nrow(EWAShand))
pathway[which(rownames(EWAShand) %in% CpGs_w1Mb)] <- "handCpGs_w1Mb"
table(pathway)
pvalues <- data.frame(results$P.value,pathway)
handCpGs_w1Mb <- pvalues[,1]
handCpGs_w1Mb [which(!pvalues$pathway=="handCpGs_w1Mb")] <- NA
All_other_CpGs <- pvalues[,1]
All_other_CpGs[which(!pvalues$pathway=="other")] <- NA
pvals <- cbind(handCpGs_w1Mb,All_other_CpGs)
handCpGs_w1Mb <- rep(0,nrow(pvals))
handCpGs_w1Mb[which(pathway=="handCpGs_w1Mb")] <- 1
dataenrich <- data.frame(results$P.value,abs(results$t),handCpGs_w1Mb)
colnames(dataenrich) [1:2] <- c("P.value","t")
rownames(dataenrich) <- rownames(results)

# regression model
linReg0a <- lm(formula = t ~ handCpGs_w1Mb, data = dataenrich)
summary(linReg0a)
out <- summary(linReg0a$coefficients)

# bootstrapping to calculate p-value
library(simpleboot)
boot0a <- lm.boot(linReg0a,R = 2000)
summary(boot0a)




