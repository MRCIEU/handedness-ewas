pgs.dat <- list()
pgs.dat$child <- alspac.table[,c("handed",
                                 "pgs.child",
                                 "sex",
                                 paste0("ancestry.child.pc", 1:10))]
colnames(pgs.dat$child) <- c("handed","pgs","sex",paste0("pc",1:10))
pgs.dat$child$sex <- sign(pgs.dat$child$sex == "M")
pgs.dat$child$id <- alspac.table$alnqlet
pgs.dat$adult <- alspac.table[,c("handed.mom",
                                 "pgs.mom",
                                 paste0("ancestry.mom.pc", 1:10))]
colnames(pgs.dat$adult) <- c("handed","pgs",paste0("pc",1:10))
pgs.dat$adult$id <- alspac.table$aln
pgs.dat <- c(sapply(child.time.points, function(time) pgs.dat$child, simplify=F),
             sapply(adult.time.points, function(time) pgs.dat$adult, simplify=F))
pgs.dat <- lapply(pgs.dat, function(dat) {
    dat$handed <- sign(dat$handed == "left")
    dat
})
