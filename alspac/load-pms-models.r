## NTR meth score summary statistics (adults)
load(file.path(ntr.dir, "data", "NTR_ewas_hand01_adults_full.RData"))
adult.stats <- results

## NTR meth score summary statistics (children)
load(file.path(ntr.dir, "data", "NTR_ewas_hand01_F9_full.RData"), verbose=T)
child.stats <- results


for (threshold in thresholds) {
    stats <- child.stats[which(child.stats$pvalue < threshold),]
    meffonym.add.model(paste0("child.handedness",threshold),
                       c("intercept", stats$name),
                       c(0,stats$estimates),
                       "child handedness")
    stats <- adult.stats[which(adult.stats$pvalue < threshold),]
    meffonym.add.model(paste0("adult.handedness",threshold),
                       c("intercept", rownames(stats)),
                       c(0,stats$estimates),
                       "adult handedness")
}



## check
models <- meffonym.models()
models <- models[grep("handedness", models)]
sapply(models, function(model) length(meffonym.get.model(model)$coefficients))
## adult.handedness0.001   adult.handedness0.1 adult.handedness1e-05 
##                   415                 43791                     2 
## child.handedness0.001   child.handedness0.1 child.handedness1e-05 
##                  1322                103089                    26 
