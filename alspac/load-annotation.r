annotation <- eval.save({
    annotation.package <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    
    load.packages(annotation.package)
    data(list=annotation.package)
    data(Locations)
    data(Other)
    annotation <- as.data.frame(Locations)
    annotation <- cbind(annotation, as.data.frame(Other))
    
    
    
    annotation.package <- "IlluminaHumanMethylationEPICanno.ilm10b2.hg19"
    
    load.packages(annotation.package)
    data(list=annotation.package)
    data(Locations)
    data(Other)
    annotation.epic <- as.data.frame(Locations)
    annotation.epic <- cbind(annotation.epic, as.data.frame(Other))
    
   list("450"=annotation, "epic"=annotation.epic)
}, "annotation")
