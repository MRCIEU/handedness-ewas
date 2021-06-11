load.packages <- function(cran=NULL,bioc=NULL){
    installed <- installed.packages()[,"Package"]
    new.cran <- setdiff(cran, installed)
    new.bioc <- setdiff(bioc, installed)

    if (length(new.cran) > 0) {
        cat(date(), "Installing packages:", paste(new.cran, collapse=","), "\n")
        install.packages(new.cran, dependencies=TRUE)
    }

    if (length(new.bioc) > 0) {
        cat(date(), "Installing packages:", paste0(new.bioc, collapse=","), "\n")
        source("http://bioconductor.org/biocLite.R")
        biocLite(new.bioc, suppressUpdates=T, suppressAutoUpdate=T)
    }
    invisible(sapply(c(cran,bioc), library, character.only = TRUE))
}
