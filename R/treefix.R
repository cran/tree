# file treefix.R copyright (C) 1994-2000 B. D. Ripley
#
prune.tree <-
    function(tree, k = NULL, best = NULL, newdata, nwts,
             method = c("deviance", "misclass"),
             loss = 1-diag(nc), eps = 1e-3)
{
    if(inherits(tree, "singlenode")) stop("Can't prune singlenode tree")
    if(!inherits(tree, "tree")) stop("Not legitimate tree")
    method <- match.arg(method)
    nc <- length(attr(tree, "ylevels"))
    if(method == "misclass" & !nc)
        stop("misclass only for classification trees")
    frame <- tree$frame
    node <- row.names(frame)
    nodes <- as.numeric(node)
    nnode <- length(node)
    ndim <- ceiling(nnode/2)

    if(is.null(y <- tree$y))
        y <- model.extract(model.frame.tree(tree), "response")
    if(is.null(w <- tree$weights))
        w <- model.extract(model.frame.tree(tree), "weights")
    if(is.null(w)) w <- rep(1, length(y))
    if(method == "misclass") {
        Z <- .C("VR_dev1",
                as.integer(nnode),
                as.integer(nodes),
                integer(nnode),
                dev = double(nnode),
                sdev = double(nnode),
                as.integer(y),
                as.integer(length(y)),
                as.integer(frame$yval),
                as.integer(tree$where),
                as.double(w),
                as.integer(nc), as.double(loss)
                )
        dev <- Z$dev; sdev <- Z$sdev
    } else {
        dev <- tree$frame$dev
        if(!nc) {
            sdev <- .C("VR_dev3",
                       as.integer(nnode),
                       as.integer(nodes),
                       integer(nnode),
                       dev = double(nnode),
                       sdev = double(nnode),
                       as.double(y),
                       as.integer(length(y)),
                       as.double(frame$yval),
                       as.integer(tree$where),
                       as.double(w),
                       )$sdev
        } else  {
            sdev <- -2 * .C("VR_dev2",
                            as.integer(nnode),
                            as.integer(nodes),
                            integer(nnode),
                            dev = double(nnode),
                            sdev = double(nnode),
                            as.integer(y),
                            as.integer(length(y)),
                            as.double(frame$yprob),
                            as.integer(tree$where),
                            as.double(w),
                            )$sdev
        }
    }
    if(missing(newdata) || is.null(newdata)) {
        ndev <- dev
        nsdev <- sdev
    } else {
        if(is.null(attr(newdata, "terms")))
            nd <- model.frame(tree$terms, newdata, na.action=na.pass)
        else nd <- newdata
        y <- model.extract(nd, "response")
        if(missing(nwts)) nwts <- rep(1, length(y))
        where <- pred1.tree(tree, tree.matrix(nd))
        if(method == "misclass") {
            Z <- .C("VR_dev1",
                    as.integer(nnode),
                    as.integer(nodes),
                    integer(nnode),
                    dev = double(nnode),
                    sdev = double(nnode),
                    as.integer(y),
                    as.integer(length(y)),
                    as.integer(frame$yval),
                    as.integer(where),
                    as.double(nwts),
                    as.integer(nc), as.double(loss)
                    )
            ndev <- Z$dev; nsdev <- Z$sdev
        } else {
            if(!nc) {
                Z <- .C("VR_dev3",
                        as.integer(nnode),
                        as.integer(nodes),
                        integer(nnode),
                        dev = double(nnode),
                        sdev = double(nnode),
                        as.double(y),
                        as.integer(length(y)),
                        as.double(frame$yval),
                        as.integer(where),
                        as.double(nwts),
                        )
                ndev <- Z$dev; nsdev <- Z$sdev
            } else {
                yp <- frame$yprob
                yp[yp==0] <- max(0,eps)
                Z <- .C("VR_dev2",
                        as.integer(nnode),
                        as.integer(nodes),
                        integer(nnode),
                        dev = double(nnode),
                        sdev = double(nnode),
                        as.integer(y),
                        as.integer(length(y)),
                        as.double(yp),
                        as.integer(where),
                        as.double(nwts),
                        )
                ndev <- -2 * Z$dev; nsdev <- -2 *Z$sdev
            }
        }
    }
    zp <- .C("VR_prune2",
             n=as.integer(nnode),
             as.integer(nodes),
             as.integer(frame$var == "<leaf>"),
             as.double(dev), as.double(sdev),
             as.double(ndev), as.double(nsdev),
             keep=integer(nnode),
             as.integer(order(nodes)),
             double(nnode),
             integer(nnode),
             double(nnode),
             alpha=double(ndim),
             inode=integer(ndim),
             size=integer(ndim),
             deviance=double(ndim),
             newdev=double(ndim)
             )
    n <- zp$n
    alpha <- zp$alpha[1:n]
    size <- zp$size[1:n]
    index <- 0
    if(missing(k) || is.null(k)) {
        ind <- drop(outer(unique(alpha), alpha, ">=") %*% rep(1, length(alpha)))
        k <- alpha[ind]
        k[1] <- -Inf
        deviance <- zp$newdev[ind]
        size <- size[ind]
        if(!missing(best) && !is.null(best)) {
            index <- ind[sum(best <= size)]
            if(length(index) == 0) {
                warning("best is bigger than tree size")
                index <- 1
            }
        }
    } else {
        if(length(k) == 1) index <- sum(k >= alpha)
        else {
            k <- pmax(k, -1e+100)
            ind <- drop(outer(k, alpha, ">=") %*% rep(1, length(alpha)))
            deviance <- zp$newdev[ind]
            size <- size[ind]
        }
    }
    if(index == 1) return(tree)
    if(index > 1) {
        pnodes <- zp$inode[-1]
        tree <- snip.tree(tree, pnodes[seq(index-1)])
        tree$call$tree <- match.call()$tree
        return(tree)
    }
    obj <- list(size = size, dev = deviance, k = k, method = method)
    class(obj) <- c("prune", "tree.sequence")
    obj
}

predict.tree <-
    function(object, newdata = list(),
             type = c("vector", "tree", "class", "where"),
             split = FALSE, nwts, eps = 1e-3, ...)
{
    which.is.max <- function(x)
    {
        y <- seq(length(x))[x == max(x)]
        if(length(y) > 1) sample(y,1)
        else y
    }

    pred2.tree  <- function(tree, x)
    {
        frame <- tree$frame
        if(!length(frame$yprob)) stop("only for classification trees")
        dimx <- dim(x)
        ypred <- .C("VR_pred2",
                    as.double(x),
                    as.integer(unclass(frame$var) - 1),#0 denotes leaf node
                    as.character(frame$splits[, "cutleft"]),
                    as.character(frame$splits[, "cutright"]),
                    as.integer(sapply(attr(tree, "xlevels"), length)),
                    as.integer(row.names(frame)),
                    as.integer(frame$n),
                    as.integer(nf <- dim(frame)[1]),
                    as.integer(dimx[1]),
                    where = double(nf*dimx[1]),
                    NAOK = TRUE)
        ypred <- matrix(ypred$where, nf)
        dimnames(ypred) <- list(row.names(frame),dimnames(x)[[1]])
        ypred
    }

    if(!inherits(object, "tree") && !inherits(object, "singlenode"))
        stop("Not legitimate tree")
    type <- match.arg(type)
    if(type == "class" && is.null(attr(object, "ylevels")))
        stop("type class only for classification trees")
    if(missing(newdata) || is.null(newdata) & type == "tree")
        return(object)                  #idiot proofing
    if(missing(newdata) || is.null(newdata)) {
        where <- object$where
        newdata <- model.frame.tree(object)
        if(!is.null(w <- object$call$weights))
            nwts <- model.extract(model.frame.tree(object), "weights")
    } else {
        if(is.null(attr(newdata, "terms"))) {
            # newdata is not a model frame.
            Terms <- object$terms
            if(type == "tree") {
                # test if response can be extracted from newdata
                response.vars <- all.vars(formula(Terms)[[2]])
                response.exists <-
                    sapply(response.names, function(nm, newdata)
                           eval(local = newdata,
                                substitute(exists(nm), list(nm=nm))),
                           newdata)
                if(!all(response.exists)) Terms <- delete.response(Terms)
            } else Terms <- delete.response(Terms)
            newdata <- model.frame(Terms, newdata, na.action = na.pass)
        }
        where <- pred1.tree(object, tree.matrix(newdata))
    }
    if(type == "where") return(where)
    frame <- object$frame
    node <- row.names(frame)
    nodes <- as.numeric(node)
    nnode <- length(node)
    if(type != "tree")
        if(is.null(lev <- attr(object, "ylevels"))) {
            if(!split) {
                frame <- frame$yval[where]
                names(frame) <- names(where)
                return(frame)
            } else {
                where <- pred2.tree(object, tree.matrix(newdata))
                leaf <- frame$var=="<leaf>"
                frame <- t(where[leaf, , drop = FALSE]) %*% frame$y[leaf]
                names(frame) <- names(where)
                return(frame)
            }
        } else {
            if(!split) {
                pr <- frame$yprob[where,  , drop = FALSE]
                dimnames(pr)[[1]] <- names(where)
            } else {
                where <- pred2.tree(object, tree.matrix(newdata))
                leaf <- frame$var=="<leaf>"
                pr <- t(where[leaf,,drop = FALSE]) %*% frame$yprob[leaf,,drop=FALSE]
                dimnames(pr) <- list(names(where), lev)
            }
            if(type=="class") {
                cl <- apply(pr, 1, which.is.max)
                return(factor(lev[cl], levels=lev))
            } else return(pr)
        }
    # now must be type = "tree"
    which <- descendants(as.numeric(row.names(frame)))[, where, drop = FALSE]
    if(!all(response.exists)) dev <- rep(NA, nrow(frame))
    else {
        y <- model.extract(newdata, "response")
        if(missing(nwts)) nwts <- rep(1, length(y))
        if(!length(attr(object, "ylevels"))) {
#
#  handle NAs in y separately.
#
            drp <- is.na(y); nwts[drp] <- 0; y[drp] <- 0
            dev <- .C("VR_dev3",
                      as.integer(nnode),
                      as.integer(nodes),
                      integer(nnode),
                      dev = double(nnode),
                      sdev = double(nnode),
                      as.double(y),
                      as.integer(length(y)),
                      as.double(frame$yval),
                      as.integer(where),
                      as.double(nwts)
                      )$dev
            dev[which %*% drp > 0] <- NA
        } else {
            yp <- frame$yprob
            yp[yp==0] <- max(0,eps)
            drp <- is.na(y); nwts[drp] <- 0; y[drp] <- levels(y)[1]
            dev <- -2 * .C("VR_dev2",
                           as.integer(nnode),
                           as.integer(nodes),
                           integer(nnode),
                           dev = double(nnode),
                           sdev = double(nnode),
                           as.integer(y),
                           as.integer(length(y)),
                           as.double(yp),
                           as.integer(where),
                           as.double(nwts)
                           )$dev
            dev[which %*% drp > 0] <- NA
        }
    }
    object$frame$dev <- as.vector(dev)
    object$frame$n <- as.vector(which %*% rep(1, length(where)))
    object$where <- where
    object$call <- match.call()
    object$y <- object$x <- NULL
    object
}

pred1.tree <- function(tree, x)
{
    frame <- tree$frame
    dimx <- dim(x)
    ypred <- .C("VR_pred1",
                as.double(x),
                as.integer(unclass(frame$var) - 1),#0 denotes leaf node
                as.character(frame$splits[, "cutleft"]),
                as.character(frame$splits[, "cutright"]),
                as.integer(sapply(attr(tree, "xlevels"), length)),
                as.integer(row.names(frame)),
                as.integer(frame$n),
                as.integer(nf <- dim(frame)[1]),
                as.integer(dimx[1]),
                as.integer(dimx[2]),
                where = integer(dimx[1]),
                NAOK = TRUE)
    ypred <- ypred$where
    names(ypred) <- dimnames(x)[[1]]
    ypred
}


na.tree.replace <- function(frame)
{
    if(!is.null(j <- attr(attr(frame, "terms"), "response"))) {
        x <- frame[[j]]
        pos <- is.na(x)
        if(any(pos)) {
            frame <- frame[!pos,  , drop = FALSE]
            warning(paste(sum(pos),
                          "observations omitted due to missing values in the response"))
        }
    }
    if(!is.na(j <- match("(weights)", names(frame)))) {
        x <- frame[[j]]
        pos <- is.na(x)
        if(any(pos)) {
            frame <- frame[!pos,  , drop = FALSE]
            warning(paste(sum(pos),
                          "observations omitted due to missing values in the supplied weights"))
        }
    }
    vars <- names(frame)
    names(vars) <- vars
    for(j in names(frame)) {
        x <- frame[[j]]
        pos <- is.na(x)
        if(any(pos)) {
            if(!length(levels(x)))
                stop(paste("continuous variable", j, "contained NAs"))
            else {
                cl <- class(x)
                class(x) <- NULL
                lev <- c(attr(x, "levels"), "NA")
                x[pos] <- length(lev)
                levels(x) <- lev
                class(x) <- cl
            }
            frame[[j]] <- x
        }
    }
    frame
}

prune.misclass <- function(tree, ...)
{
    oc <- match.call()
    oc$method <- "misclass"
    oc[[1]] <- as.name("prune.tree")
    eval(oc, parent.frame())
}
