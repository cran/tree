\name{plot.tree.sequence}
\alias{plot.tree.sequence}
\title{
Plot a Tree Sequence
}
\description{
Allows the user to plot a tree sequence.
}
\usage{
plot.tree.sequence(x,  \dots, order=c("increasing", "decreasing"))
}
\arguments{
\item{x}{
object of class \code{tree.sequence}.
This is assumed to be the result of some function that produces
an object with the same named components (\code{size}, \code{deviance}, \code{k})
as that returned by \code{prune.tree()} and \code{shrink.tree()}.
}
\item{order}{
of \code{size} on the plot.  Use \code{"decreasing"} for the natural ordering of \code{k}
and the amount of pruning. Only the first character is needed.
}}
\section{Side Effects}{
plots deviance or number of misclassifications (or total loss)
versus size for a sequence of trees.
}
\details{
This function is a method for the generic function
\code{plot()} for class \code{tree.sequence}.
It can be invoked by calling \code{plot(x)} for an
object \code{x} of the appropriate class, or directly by
calling \code{plot.tree.sequence(x)} regardless of the
class of the object.
}
\seealso{
\code{\link{plot.tree}}.
}
\examples{
zp <- prune.tree(z)
plot(zp, type="s")
zs <- shrink.tree(z)
lines(zs$size, zs$dev)
}
\keyword{tree}
% Converted by Sd2Rd version 0.3-1.
