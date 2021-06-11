stop("run the following command on bluecrystalp3")
## bash extract-genotypes.sh
##   out: children.{bed,fam,bim} 
##        mother.{bed,fam,bim} files
##        snps.bim.gz
## copy these files to geno.dir/genotype

stopifnot(file.exists(file.path(geno.dir, "genotype",
                                c("children.bed",
                                  "mothers.bed",
                                  "snps.bim.gz"))))
