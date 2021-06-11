kableone <- function(x, ...) {
    capture.output(x <- print(x))
    knitr::kable(x, ...)
}
