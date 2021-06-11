require(gdsfmt)

load.aries <- function(path, featureset, time=NULL, detectionp=F) {
    ## showfile.gds(closeall=T)
    ## ls.gdsn(node)
    samplesheet.filename <- file.path(path, "samplesheet", "samplesheet.csv")
    if (!file.exists(samplesheet.filename))
        stop("invalid path for ARIES data: ", path)
    
    samplesheet <- read.csv(samplesheet.filename, stringsAsFactors=F)

    if (!is.null(time)) {
        if (!all(time %in% samplesheet$time_point))
            stop("invalid time point for ARIES: ",
                 paste(setdiff(time,samplesheet$time_point), collapse="/"),
                 " (valid times = ",
                 paste(sort(unique(samplesheet$time_point)), collapse=", "), ")")
        
        idx <- which(samplesheet$time_point %in% time)
        samplesheet <- samplesheet[idx,,drop=F]
    }

    gds.filename <- file.path(path, "betas", paste(featureset, "gds", sep="."))
    if (!file.exists(gds.filename))
        stop(featureset, " is an invalid bead chip or feature set",
             " (possible values = ",
             paste(sub(".gds", "", list.files(file.path(path, "betas"), ".gds$")), collapse=", "),
             ")")
      
    gds.file <- openfn.gds(gds.filename)  
    samples <- read.gdsn(index.gdsn(gds.file, "col.names"))
    idx <- which(samples %in% samplesheet$Sample_Name)
    if (length(idx) == 0) {
        closefn.gds(gds.file)
        stop("there are no samples in '", featureset, "' format",
             ifelse(!is.null(time),
                    paste(" at time point", paste(time, collapse="/")),
                    ""))
    }
    start <- min(idx)
    count <- max(idx)-min(idx)+1
    
    matrix.node <- index.gdsn(gds.file, "matrix")
    meth <- read.gdsn(matrix.node,
                      start=c(1,start),
                      count=c(-1,count))
    
    rownames(meth) <- read.gdsn(index.gdsn(gds.file, "row.names"))
    colnames(meth) <- samples[start:(start+count-1)]
    closefn.gds(gds.file)
    
    samplesheet <- samplesheet[which(samplesheet$Sample_Name %in% colnames(meth)),]
    meth <- meth[,match(samplesheet$Sample_Name, colnames(meth))]

    control.matrix <- read.table(file.path(path, "control_matrix", paste(featureset, "txt", sep=".")),
                                 row.names=1, header=T, sep="\t")
    control.matrix <- control.matrix[colnames(meth),]

    cellcount.filenames <- list.files(file.path(path, "derived", "cellcounts"), ".txt$", full.names=T)
    counts <- sapply(cellcount.filenames, read.table, sep="\t", row.names=1, header=T, simplify=F)
    names(counts) <- sub(".txt", "", sapply(names(counts), basename))
    overlaps <- sapply(counts, function(counts) all(colnames(meth) %in% rownames(counts)))
    if (!all(overlaps))
        counts <- counts[-which(!overlaps)]
    counts <- lapply(counts, function(counts) counts[colnames(meth),])
    names(counts) <- gsub(" ", "-", names(counts))

    list(methylation=meth,
         samples=samplesheet,
         control.matrix=control.matrix,
         cellcounts=counts)
}
