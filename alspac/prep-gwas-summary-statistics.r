filename <- file.path(pgs.dir, "handedness-gwas-summary-stats-STANDARD.txt")
if (!file.exists(filename)) {
    stats <- read.table(file.path(ntr.dir, "data", "Meta_RvL_no23andMe_noALSPAC_noNTR1.txt.sorted.tsv"),
                        header=T, sep="\t", stringsAsFactors=F)   ## 3 minutes

    genotype.reference <- eval.save({
        read.delim(file.path(geno.dir, "genotype", "snps.bim.gz"),
                   sep=" ",
                   header=F,
                   stringsAsFactors=F)
    }, "genotype-reference") ## 5 minutes
    
    colnames(genotype.reference) <- c("chr", "rsid", "whoknows", "pos", "allele1", "allele2")
    
    idx <- match(stats$MarkerName, genotype.reference$rsid)
    ## What proportion are matched?
    sum(!is.na(idx))/nrow(stats)
    ## [1] 0.7979837
    stats$chr <- genotype.reference$chr[idx]
    stats$pos <- genotype.reference$pos[idx]
    stats$allele1.ref <- genotype.reference$allele1[idx]
    stats$allele2.ref <- genotype.reference$allele2[idx]
    
    stats <- stats[!is.na(stats$chr),]
    
    idx <- which(toupper(stats$Allele1) == stats$allele1.ref | toupper(stats$Allele1) == stats$allele2.ref)
    ## What proportion have matching alleles?
    length(idx)/nrow(stats)
    ## [1] 0.9370226
    stats <- stats[idx,]
    
    idx <- which(toupper(stats$Allele1) != stats$allele1.ref)
    stats$Freq1[idx] <- 1-stats$Freq1[idx]
    stats$Zscore[idx] <- -stats$Zscore[idx]
    
    stats$info <- 1
    
    stats <- stats[,c("chr","pos","allele1.ref", "allele2.ref", "Freq1", "info", "MarkerName", "P.value", "Zscore")] 
    
    colnames(stats) <- c("chr","pos","ref","alt","reffrq","info","rs","pval","effalt")
    
    ## Note that GWAS z-score is being passed on as the effect ('effalt').
    ## We only have z-score and p-value so we don't actually have an
    ## odds ratio which is what LDpred expects.
    ## Seems complicated to generate an estimate.
    ## https://seanharrisonblog.com/2020/04/11/estimating-an-odds-ratio-from-a-gwas-only-reporting-the-p-value/
    ## However, should be fine as it has a standard normal distribution.
    
    write.table(stats, file=filename,sep="\t", col.names=T, row.names=F, quote=F)
    ## STANDARD format (https://github.com/bvilhjal/ldpred/wiki/Q-and-A)
    ## chr pos ref alt reffrq info rs pval effalt
    ## chr1 1020428 C T 0.85083 0.98732 rs6687776 0.0587 -0.0100048507289348
}
