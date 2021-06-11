
## calculate scores for each time point and each methylation model
scores <- sapply(time.points, function(time) {
    eval.save({
        cat(date(), "methylation score at time-point", time, "\n")
        aries <- my.load.aries(aries.dir, time)
        models <- meffonym.models()
        if (time %in% child.time.points) {
            pms.models <- models[grepl("child.handedness", models)]
            id <- aries$samples$alnqlet
        } else {        
            pms.models <- models[grepl("adult.handedness", models)]
            id <- aries$samples$aln
        }
        scores <- sapply(pms.models, function(model) {
            cat(time, model,
                length(intersect(rownames(aries$methylation),
                                 names(meffonym.get.model(model)$coefficients))), "\n")
            ms <- meffonym.score(aries$methylation, model)
            scale(ms$score)
        })
        rownames(scores) <- id
        scores
    }, paste("polymethylation-scores", time,sep="-"))
}, simplify=F)

## cord child.handedness0.001 546 
## cord child.handedness0.1 46698 
## cord child.handedness1e-05 7 
## F7 child.handedness0.001 546 
## F7 child.handedness0.1 46698 
## F7 child.handedness1e-05 7 
## 15up child.handedness0.001 546 
## 15up child.handedness0.1 46698 
## 15up child.handedness1e-05 7 
## F24 child.handedness0.001 1308 
## F24 child.handedness0.1 102269 
## F24 child.handedness1e-05 26 
## antenatal adult.handedness0.001 413 
## antenatal adult.handedness0.1 43628 
## antenatal adult.handedness1e-05 2 
## FOM adult.handedness0.001 413 
## FOM adult.handedness0.1 43628 
## FOM adult.handedness1e-05 2 

## add them to the alspac.table data frame
for (time in names(scores)) {
    colnames(scores[[time]]) <- paste(time,
                                      colnames(scores[[time]]),
                                      sep=".")
    if (time %in% child.time.points) 
        idx <- match(alspac.table$alnqlet, rownames(scores[[time]]))
    else
        idx <- match(alspac.table$aln, rownames(scores[[time]]))
    for (name in colnames(scores[[time]]))
        alspac.table[[name]] <- scores[[time]][idx,name]
}


## check
varnames <- colnames(alspac.table)[grep("handedness", colnames(alspac.table))]
sapply(varnames, function(varname) sum(!is.na(alspac.table[[varname]])))
##      cord.child.handedness0.001        cord.child.handedness0.1 
##                             905                             905 
##      cord.child.handedness1e-05        F7.child.handedness0.001 
##                             905                             969 
##          F7.child.handedness0.1        F7.child.handedness1e-05 
##                             969                             969 
##      15up.child.handedness0.001        15up.child.handedness0.1 
##                             970                             970 
##      15up.child.handedness1e-05       F24.child.handedness0.001 
##                             970                             569 
##         F24.child.handedness0.1       F24.child.handedness1e-05 
##                             569                             569 
## antenatal.adult.handedness0.001   antenatal.adult.handedness0.1 
##                             964                             964 
## antenatal.adult.handedness1e-05       FOM.adult.handedness0.001 
##                             964                             978 
##         FOM.adult.handedness0.1       FOM.adult.handedness1e-05 
##                             978                             978 

sapply(scores, nrow)
## cord        F7      15up       F24 antenatal       FOM 
##  914       978       981       575       987       992 

## should be about 1 because they are scaled
quantile(sapply(varnames, function(varname) var(alspac.table[[varname]], na.rm=T)))
##        0%       25%       50%       75%      100% 
## 0.9823845 0.9935760 1.0007303 1.0050169 1.0116309 
quantile(sapply(varnames, function(varname) mean(alspac.table[[varname]], na.rm=T)))
##           0%          25%          50%          75%         100% 
## -0.008973460 -0.003640274 -0.001909079  0.003245493  0.010631194 

