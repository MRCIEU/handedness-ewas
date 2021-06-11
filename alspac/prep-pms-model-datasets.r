pms.dat <- sapply(time.points, function(time) {
    sub.dat <- sub.dat[[time]]
    ewas.dat <- ewas.dat[[time]]
    sub.dat$id <- NULL
    dat <- cbind(sub.dat, ewas.dat)
    if (time %in% child.time.points) {
        idx <- match(dat$id, alspac.table$alnqlet)
    } else {
        idx <- match(dat$id, alspac.table$aln)
    }
    pms.var.idx <- which(grepl(time, colnames(alspac.table))
                         & grepl("handedness", colnames(alspac.table)))
    pms.vars <- colnames(alspac.table)[pms.var.idx]
    names(pms.vars) <- sub(paste0(time,"."), "", pms.vars)
    lapply(pms.vars, function(var) {
        dat$pms <- alspac.table[[var]][idx]
        dat
    })
}, simplify=F)
