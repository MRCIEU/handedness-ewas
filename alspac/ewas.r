ewas.ret <- sapply(c("cord","F7","15up","F24","antenatal","FOM","FOF"), function(time) {
    sapply(c("basic","full"), function(model) {
        eval.save({
            ## DNA methylation
            aries <- my.load.aries(aries.dir, time)
            ## select cell count estimate
            ref <- ifelse(time == "cord", "andrews-and-bakulski-cord-blood", "blood-gse35069")
            ## match DNA methylation to ALSPAC
            if (time %in% c("antenatal","FOM","FOF"))
                alspac.idx <- match(aries$samples$aln, alspac.table$aln)
            else
                alspac.idx <- match(aries$samples$alnqlet, alspac.table$alnqlet)
            ## basic covariates
            data <- data.frame(age=aries$samples$age,
                               sex=aries$samples$sex,
                               aries$cellcounts[[ref]])
            if (time == "cord")
                data$age <- alspac.table$gestational.age[alspac.idx]
            ## handed variable
            if (time %in% c("antenatal","FOM"))
                data$handed <- alspac.table$handed.mom[alspac.idx]
            else if (time == "FOF")
                data$handed <- alspac.table$handed.partner[alspac.idx]
            else
                data$handed <- alspac.table$handed[alspac.idx]
            if (length(unique(na.omit(data$sex))) < 2)
                data$sex <- NULL        
            data$handed[which(data$handed == "mixed")] <- NA
            ## full model covariates
            if (model == "full") {
                if (time %in% c("cord","F7")) {
                    data$birthweight <- alspac.table$birthweight[alspac.idx]
                    data$maternal.smoking <- alspac.table$maternal.smoking[alspac.idx]
                }
                if (time == "F7")
                    data$gestational.age <- alspac.table$gestational.age[alspac.idx]
                if (time == "15up")
                    data$bmi <- alspac.table$bmi.17[alspac.idx]
                else if (time == "F24")
                    data$bmi <- alspac.table$bmi.24[alspac.idx]
                else if (time == "antenatal")
                    data$bmi <- alspac.table$bmi.mother.ant[alspac.idx]
                else if (time == "FOM")
                    data$bmi <- alspac.table$bmi.mother.fom[alspac.idx]
                else if (time == "FOF")
                    data$bmi <- alspac.table$bmi.partner.fof[alspac.idx]
                if (time %in% c("15up","F24","antenatal","FOM","FOF"))
                    data$smoking <- aries$methylation["cg05575921",]
            }
            ## run EWAS
            ewaff.sites("methylation ~ .",
                        "handed",
                        methylation=aries$methylation,
                        data=data,
                        method="limma",
                        generate.confounders="smartsva",
                        n.confounders=20)
        }, paste("ewas", "handed", time, model, sep="-"))
    }, simplify=F) 
}, simplify=F) ## total: 4.5h (~20 minutes each)


## perform EWAS in parents
ewas.ret$parents <- sapply(c("basic","full"), function(model) {
    eval.save({
        aries <- my.load.aries(aries.dir, "parents")
        cols <- colnames(ewas.ret$FOM[[model]]$data)
        data <- rbind(ewas.ret$FOM[[model]]$data,
                      ewas.ret$FOF[[model]]$data[,cols])
        idx <- match(rownames(data), aries$samples$Sample_Name)        
        ewaff.sites("methylation ~ .",
                    "handed",
                    methylation=aries$methylation[,idx],
                    data=data,
                    method="limma",
                    generate.confounders="smartsva",
                    n.confounders=20)
    }, paste("ewas", "handed", "parents", model, sep="-"))
}, simplify=F) ## total: 30 minutes


## tell each EWAS it's ARIES time point and model 
for (i in 1:length(ewas.ret))
    for (j in 1:length(ewas.ret[[i]])) {
        ewas.ret[[i]][[j]]$time <- names(ewas.ret)[i]
        ewas.ret[[i]][[j]]$model <- names(ewas.ret[[i]])[j]
    }

## sort ewas stats by genomic position
for (time in names(ewas.ret)) {
    for (model in names(ewas.ret[[time]])) {
        ewas <- ewas.ret[[time]][[model]]$table
        anno = annotation[["450"]]
        if (time == "F24")
            anno = annotation[["epic"]]
        sites <- rownames(ewas)
        idx <- match(sites, rownames(anno))
        ewas$chr <- anno$chr[idx]
        ewas$pos <- anno$pos[idx]
        ewas <- ewas[order(ewas$chr, ewas$pos),]
        ewas.ret[[time]][[model]]$table <- ewas
    }
}

## generate ewas reports
for (time in names(ewas.ret)) {
    for (model in names(ewas.ret[[time]])) {
        summary <- eval.save({
            aries <- my.load.aries(aries.dir, time)
            ewas <- ewas.ret[[time]][[model]]
            sites <- rownames(ewas$table)
            summary <- ewaff.summary(ewas, chr=ewas$table$chr, pos=ewas$table$pos,
                                     methylation=aries$methylation[sites,])
            output.file <- file.path(output.dir, paste0("report-", time, "-", model, ".html"))
            ewaff.report(summary, output.file=output.file, author="M Suderman", study="Handedness")
            summary
        }, paste("ewas-summary", time, model, sep="-"))
    }    
} ## 45 minutes


## save EWAS summary statistics
for (ewas.by.time in ewas.ret) {
    for (ewas in ewas.by.time) {
        filename <- file.path(output.dir, paste0("ewas-", "handed-", ewas$time, "-", ewas$model, ".csv"))
        if (!file.exists(filename))
            write.csv(ewas$table, row.names=T,file=filename)
    }
}


## save summary statistics for top 100 associations
dir.create(top.dir <- file.path(output.dir, "top100"))
for (ewas.by.time in ewas.ret) {
    for (ewas in ewas.by.time) {
        filename <- file.path(top.dir, paste0("ewas-", "handed-", ewas$time, "-", ewas$model, ".csv"))
        stats <- ewas$table
        stats$fdr <- p.adjust(stats$p.value, "fdr")
        if (nrow(stats) < 500000)
            ann <- annotation[["450"]]
        else
            ann <- annotation[["epic"]]
        idx <- order(stats$p.value)[1:100]
        stats <- stats[idx,]
        idx <- match(rownames(stats), rownames(ann))
        stats$gene <- ann[["UCSC_RefGene_Name"]][idx]
        stats$region <- ann[["UCSC_RefGene_Group"]][idx]
        stats$direction <- ifelse(stats$estimate > 0, "--", "++")
        write.csv(stats[,c("chr","pos","gene","region","n",
                           "estimate","se","p.value","fdr","direction")],
                  row.names=T,file=filename)
    }
}


