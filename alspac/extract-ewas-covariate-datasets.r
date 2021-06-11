aries.samples <- read.csv(file.path(aries.dir, "samplesheet", "samplesheet.csv"), stringsAsFactors=F)
ewas.dat <- sapply(time.points, function(time) {
    ewas <- eval.ret(paste("ewas-handed", time, "full", sep="-"))
    ewas.dat <- as.data.frame(ewas$design)
    ewas.dat <- ewas.dat[,which(colnames(ewas.dat) != "(Intercept)")]
    ewas.dat <- ewas.dat[,grep("handed",colnames(ewas.dat),invert=T)]
    idx <- match(rownames(ewas.dat), aries.samples$Sample_Name)
    if (time %in% child.time.points) {
        ewas.dat$sexM <- NULL  ## will be contributed by PGS model data
        ewas.dat$id <- aries.samples$alnqlet[idx]            
    }
    else {
        ewas.dat$id <- aries.samples$aln[idx]
    }
    ewas.dat
})
