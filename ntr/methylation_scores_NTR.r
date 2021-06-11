# Left-handedness Methylation scores calculation: based on weights from ALSPAC adults EWAS and DNA methylation data from NTR adults
# V.Odintsova, J.van Dongen

# Load data
# load methylation data NTR adults with samples in rows, CpGs in columns (beta)
# load covariates NTR adults (covariates)
# load summary statistics ALSPAC adults EWAS (all_weights)

# Create weight datasets with CpGs at 3 p-value thresholds (10^-1,10^-3,10^-5)
model1 <- all_weights[which(all_weights$P < 10^-1),]
model2 <- all_weights[which(all_weights$P < 10^-3),]
model3 <- all_weights[which(all_weights$P < 10^-5),]

# Calculate scores
calc_scores <- function (model) {
  OL <- intersect(colnames(beta),rownames(model))
  weights <- model[OL,] 
  met <- beta[,OL] 
  Nsubjects <- nrow(met)
  NCpGs <- ncol(met)
  products <- matrix(NA,nrow=Nsubjects, ncol=NCpGs)
  for (i in 1:Nsubjects) {
    for (j in 1:NCpGs)  {
      products[i,j] <- met[i,j] * weights[j,"Estimate"]  }}
  score <- rowSums(products,na.rm=T)
  return(score)
}

Nscores <- 3
scores <-data.frame(matrix(NA,nrow=nrow(beta), ncol=Nscores))
rownames(scores) <- rownames(beta)

scores[1] <- calc_scores(model1)
scores[2] <- calc_scores(model2)
scores[3] <- calc_scores(model3)

# Standardize the scores
scaled.scores <- scale(scores)
scaled.scores <- as.data.frame(scaled.scores)
colnames(scaled.scores) <- c("10^-1", "10^-3", "10^-5")

# Save



