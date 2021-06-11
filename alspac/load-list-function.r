load.list <- function(filename) {
    envir <- sys.frame(sys.nframe())
    object_names <- load(filename, envir=envir)
    L <- lapply(object_names, get, envir=envir)
    names(L) <- object_names
    L
}
