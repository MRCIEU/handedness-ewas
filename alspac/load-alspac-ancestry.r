
ancestry.children <- read.dta(file.path(alspac.dir,
                                        "ancestry-pca",
                                        "children",
                                        "alspac10k_pca.dta"))
ancestry.children$alnqlet <- with(ancestry.children, paste0(aln,qlet))
idx <- match(alspac.table$alnqlet,ancestry.children$alnqlet)
for (var in paste0("pc", 1:10))
    alspac.table[[paste("ancestry","child",var,sep=".")]] <- ancestry.children[[var]][idx] 

ancestry.mothers <- read.dta13(file.path(alspac.dir,
                                       "ancestry-pca",
                                       "Mothers",
                                       "ancestry_PC_Mum.dta"))
idx <- match(alspac.table$aln,ancestry.mothers$aln)
colnames(ancestry.mothers) <- sub("m", "", tolower(colnames(ancestry.mothers)))
for (var in paste0("pc", 1:10))
    alspac.table[[paste("ancestry","mom",var,sep=".")]] <- ancestry.mothers[[var]][idx] 
