my.load.aries <- function(path,time) {
    nicname <- time
    if (time == "parents")
        time <- c("FOM","FOF")
    chip <- ifelse("F24" %in% time, "epic", "450")
    annotation <- annotation[[chip]]
    aries <- load.aries(path,chip,time)
    xy.sites <- rownames(annotation)[which(annotation$chr %in% c("chrX","chrY"))]
    autosome.idx <- which(!rownames(aries$methylation) %in% xy.sites)
    aries$methylation <- aries$methylation[autosome.idx,]
    aries$methylation <- ewaff.handle.outliers(aries$methylation)$methylation
    aries$annotation <- annotation[match(rownames(aries$methylation), rownames(annotation)),]
    aries$time <- nicname
    aries
}
