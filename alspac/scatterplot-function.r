scatterplot <- function (x,y,title,xlab,ylab) {
    stopifnot(is.vector(x))
    stopifnot(is.vector(y))
    stopifnot(length(x) == length(y))

    fit <- lm(y~x)
    base <- lm(y~1)
    p.value.lm <- anova(fit, base)[2, "Pr(>F)"]
    stats.desc <- paste(xlab, "\np[lm]= ",
                        format(p.value.lm,  digits = 3), sep = "")
    if (is.factor(x) || length(unique(x)) <= 20) {
        x <- as.factor(x)
        p <- (ggplot(data.frame(x=x,y=y), aes(x = x, y = y)) + 
              geom_boxplot())
    }
    else {
        p <- (ggplot(data.frame(x=x,y=y), aes(x = x, y = y)) + 
              geom_point() + geom_smooth(method = lm))
    }
    (p + ggtitle(title) + xlab(stats.desc) + ylab(ylab))
}
