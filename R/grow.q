#
#  tree/grow.q by B. D. Ripley  Copyright (C) 1994
#

na.pass <- function(x) x

tree <-
function(formula = formula(data), data = sys.parent(), weights, subset, 
	na.action = na.pass, control = tree.control(nobs, ...), method = 
	"recursive.partition", split=c("deviance", "gini"), 
	model = NULL, x = FALSE, y = TRUE, wts = TRUE, ...)
{
  if(is.null(model)) {
    model <- match.call(expand = FALSE)
    model$method <- model$model <- model$control <- 
      model$... <- model$x <- model$y <- model$wts <- model$split <- NULL
    model[[1]] <- as.name("model.frame.default")
    model <- eval(model, sys.parent())
    if(method == "model.frame") return(model)
  }
  split <- match.arg(split)
  Terms <- attr(model, "terms")
  if(any(attr(Terms, "order") > 1))
    stop("Trees cannot handle interaction terms")
  Y <- model.extract(model, "response")
  if(is.matrix(Y) && ncol(Y) > 1)
    stop("Trees cannot handle multiple responses")
  ylevels <- levels(Y)
  w <- model.extract(model, "weights")
  if(!length(w)) w <- rep(1, nrow(model))
  if(any(yna <- is.na(Y))) {
    Y[yna] <- 1                         # an innocent value
    w[yna] <- 0
  }
  offset <- attr(Terms, "offset")
  if(!is.null(offset)) {
    if(length(ylevels))
      stop("Offset not implemented for classification trees")
    offset <- model[[offset]]
    Y <- Y - offset
  }
  X <- tree.matrix(model)
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
            as.double(control$mindev),
            nnode = as.integer(0),
            where = integer(nobs),
            as.integer(control$nmax),
            as.integer(split=="gini"),
            as.integer(sapply(model, is.ordered)),
            NAOK = TRUE)
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
  attr(fit$where, "names") <- row.names(model)
  if(n > 1) class(fit) <- "tree"
  else class(fit) <- "singlenode"
  attr(fit, "xlevels") <- xlevels
  if(length(ylevels)) attr(fit, "ylevels") <- ylevels
  if(x) fit$x <- X
  if(y) fit$y <- Y
  if(wts) fit$weights <- w
  fit
}
