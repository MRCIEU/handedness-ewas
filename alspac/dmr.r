## dmr analysis
dmr.ret <- sapply(names(ewas.ret), function(time) {
    sapply(names(ewas.ret[[time]]), function(model) {
        eval.save({
            aries <- my.load.aries(aries.dir, time)
            ewas <- ewas.ret[[time]][[model]]
            sites <- rownames(ewas$table)
            with(ewas$table, dmrff(estimate, se, p.value,
                                   aries$methylation[sites,], chr, pos))
        }, paste("dmrff", time, model, sep="-"))
    }, simplify=F)
}, simplify=F) ## 85 minutes

## 'pre' objects for DMR meta-analysis
pre.ret <- sapply(names(ewas.ret), function(time) {
    sapply(names(ewas.ret[[time]]), function(model) {
        eval.save({
            aries <- my.load.aries(aries.dir, time)
            ewas <- ewas.ret[[time]][[model]]
            sites <- rownames(ewas$table)
            with(ewas$table, dmrff.pre(estimate, se, 
                                       aries$methylation[sites,], chr, pos))
        }, paste("pre", time, model, sep="-"))
    }, simplify=F)
}, simplify=F) ## 95 minutes

pre.fathers <- pre.ret[["FOF"]]
pre.mothers <- pre.ret[["FOM"]]

if (!file.exists(file.path(output.dir, "alspac-pre.rda")))
    save(pre.fathers, pre.mothers, file=file.path(output.dir, "alspac-pre.rda"))

if (!file.exists(file.path(output.dir, "alspac-7-pre.rda"))) {
    pre.7 <- pre.ret[["F7"]]
    save(pre.7, file=file.path(output.dir, "alspac-7-pre.rda"))
}

## save dmr statistics and cpg correlations
for (time in names(dmr.ret)) {
    for (model in names(dmr.ret[[time]])) {
        cat(date(), "dmrs at", time, "for model", model, "\n")
        dmrs <- dmr.ret[[time]][[model]]
        for (idx in which(dmrs$p.adjust < 0.05 & dmrs$n > 1)) {
            chr <- as.character(dmrs[idx,"chr"])
            start <- dmrs[idx,"start"] - 500
            end <- dmrs[idx,"end"] + 500

            ewas <- ewas.ret[[time]][[model]]
            if (!exists("aries") || aries$time != time)
                aries <- my.load.aries(aries.dir, ewas$time)

            stats <- ewas$table[which(ewas$table$chr == chr
                                      & ewas$table$pos >= start
                                      & ewas$table$pos <= end),]
            stats <- stats[order(stats$pos),]
            stats <- data.frame(TargetID=rownames(stats),
                                CHR=sub("chr", "", chr),
                                MAPINFO=stats$pos,
                                Pval=stats$p.value,
                                stringsAsFactors=F)
            filename <- paste("dmr-stats",
                              time,model,
                              chr,start,end, sep="-")
            filename <- paste(filename, "txt", sep=".")
            filename <- file.path(output.dir, filename)
            write.table(stats, col.names=T, row.names=F, quote=F,
                        sep="\t", file=filename)
            r <- cor(t(aries$methylation[stats$TargetID,]),use="p")
            filename <- gsub("-stats-","-corr-",filename)
            write.table(r, col.names=T, row.names=F, quote=F,
                        sep="\t", file=filename)
        }
    }
}

