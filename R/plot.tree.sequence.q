# file treefix/plot.tree.sequence.q copyright (C) 1994-8 B. D. Ripley
#
plot.tree.sequence <- function(x, ..., type = "l", ylim = range(x$dev), 
    order.=c("increasing", "decreasing"))
{
  if(missing(type) && inherits(x, "prune")) type <- "S"
  if(is.null(x$method)) x$method <- "deviance"
  order. <-  match.arg(order.)
  if(order. == "increasing") sign <- +1 else sign <- -1
  plot(sign*x$size, x$dev, axes = FALSE, xlab = "size", ylab = x$method,
       type = type, ylim = ylim, ...)
  box()
  axis(2, ...)
  xaxp <- par("xaxp")
  pos <- sign*seq(xaxp[1], xaxp[2], diff(xaxp[-3])/xaxp[3])
  if(pos[1] == 0) pos[1] <- 1
  n <- length(pos)
  maxsize <- max(x$size)
  if(pos[n] > maxsize) pos[n] <- maxsize
  axis(1, at = sign*pos, lab = pos, ...)
  # axis is broken in R
  xx <- sign*x$size
  lab <- format(signif(x$k, 2))
  ord <- order(xx)
  axis(3, at = xx[ord], lab = lab[ord], ...)
  invisible()
}
