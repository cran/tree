\name{partition.tree}
\title{Plot the Partitions of a simple Tree Model}
\usage{
partition.tree(tree, label = "yval", add = FALSE, ordvars, ...)
}
\alias{partition.tree}
\description{
  Plot the partitions of a tree involving one or two variables.
}
\arguments{
  \item{tree}{A object of class \code{"tree"}.}
  \item{label}{A character string giving the column of the frame
    component of \code{tree} to be used to label the regions.}
  \item{add}{If true, add to existing plot, otherwise start a new plot.}
  \item{ordvars}{The ordering of the variables to be used in a 2D
    plot. Specify the names in a character string of length 2; the first
    will be used on the \code{x} axis.}
  \item{\dots}{Graphical parameters.}
}
\details{
  This can be used with a regression or classification tree containing
  one or two continuous predictors (only).
  
  If the tree contains one predictor, the predicted value (a regression
  tree) or the probability of the first class (a classification tree) is
  plotted against the predictor over its range in the training set.

  If the tree contains two predictors, a plot is made of the space
  covered by those two predictors and the partition made by the tree is
  superimposed. 
}
\value{
  None.
}
\author{B. D. Ripley}
\seealso{\code{\link{tree}}}

\examples{
data(iris)
ir.tr <- tree(Species ~., iris)
ir.tr
ir.tr1 <- snip.tree(ir.tr, nodes = c(12, 7))
summary(ir.tr1)
par(pty = "s")
plot(iris[, 3],iris[, 4], type="n",
     xlab="petal length", ylab="petal width")
text(iris[, 3], iris[, 4], c("s", "c", "v")[iris[, 5]])
partition.tree(ir.tr1, add = TRUE, cex = 1.5)
}
\keyword{tree}
\keyword{hplot}
