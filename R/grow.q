#
#  tree/grow.q by B. D. Ripley  Copyright (C) 1994-2003
#
tree <-
function(formula, data, weights, subset,
         na.action = na.pass, control = tree.control(nobs, ...),
         method = "recursive.partition",
         split = c("deviance", "gini"),
         model = FALSE, x = FALSE, y = TRUE, wts = TRUE, ...)
{
    if (is.data.frame(model)) {
	m <- model
	model <- FALSE
    } else {
        m <- match.call(expand = FALSE)
        m$method <- m$model <- m$control <- m$... <- m$x <- m$y <- m$wts <-
            m$split <- NULL
        m[[1]] <- as.name("model.frame.default")
        m <- eval.parent(m)
        if(method == "model.frame") return(m)
    }
    split <- match.arg(split)
    Terms <- attr(m, "terms")
    if(any(attr(Terms, "order") > 1))
        stop("Trees cannot handle interaction terms")
    Y <- model.extract(m, "response")
    if(is.matrix(Y) && ncol(Y) > 1)
        stop("Trees cannot handle multiple responses")
    ylevels <- levels(Y)
    w <- model.extract(m, "weights")
    if(!length(w)) w <- rep(1, nrow(m))
    if(any(yna <- is.na(Y))) {
        Y[yna] <- 1                     # an innocent value
        w[yna] <- 0
    }
    offset <- attr(Terms, "offset")
    if(!is.null(offset)) {
        if(length(ylevels))
            stop("Offset not implemented for classification trees")
        offset <- m[[offset]]
        Y <- Y - offset
    }
    X <- tree.matrix(m)
    xlevels <- attr(X, "column.levels")
    if(is.null(xlevels)) {
        xlevels <- rep(list(NULL), ncol(X))
        names(xlevels) <- dimnames(X)[[2]]
    }
    nobs <- length(Y)
    if(!is.null(control$nobs) && control$nobs < nobs) {
        stop("control$nobs < number of observations in data")
    }
    fit <- .C("BDRgrow1",
              as.double(X),
              as.double(unclass(Y)),
              as.double(w),
              as.integer(c(sapply(xlevels, length), length(ylevels))),
              as.integer(rep(1, nobs)),
              as.integer(nobs),
              as.integer(ncol(X)),
              node = integer(control$nmax),
              var = integer(control$nmax),
              cutleft = character(control$nmax),
              cutright = character(control$nmax),
              n = double(control$nmax),
              dev = double(control$nmax),
              yval = double(control$nmax),
              yprob = double(max(control$nmax * length(ylevels), 1)),
              as.integer(control$minsize),
              as.integer(control$mincut),
              as.double(max(0, control$mindev)),
              nnode = as.integer(0),
              where = integer(nobs),
              as.integer(control$nmax),
              as.integer(split=="gini"),
              as.integer(sapply(m, is.ordered)),
              NAOK = TRUE,
              PACKAGE = "tree")
    n <- fit$nnode
    frame <- data.frame(fit[c("var", "n", "dev", "yval")])[1:n,  ]
    frame$var <- factor(frame$var, 0:length(xlevels),
                        c("<leaf>", names(xlevels)))
    frame$splits <-
        array(unlist(fit[c("cutleft", "cutright")]),
              c(control$nmax, 2),
              list(character(0), c("cutleft", "cutright")))[1:n, , drop = FALSE]
    if(length(ylevels)) {
        frame$yval <- factor(frame$yval, 1:length(ylevels), ylevels)
        class(frame$yval) <- class(Y)
        frame$yprob <-
            t(array(fit$yprob, c(length(ylevels), control$nmax),
                    list(ylevels, character(0)))[, 1:n, drop = FALSE])
    }
    row.names(frame) <- fit$node[1:n]
    fit <- list(frame = frame, where = fit$where, terms = Terms,
                call = match.call())
    attr(fit$where, "names") <- row.names(m)
    if(n > 1) class(fit) <- "tree" else class(fit) <- c("singlenode", "tree")
    attr(fit, "xlevels") <- xlevels
    if(length(ylevels)) attr(fit, "ylevels") <- ylevels
    if(is.logical(model) && model) fit$model <- m
    if(x) fit$x <- X
    if(y) fit$y <- Y
    if(wts) fit$weights <- w
    fit
}
