\name{tree}
\title{Fit a Classification or Regression Tree}
\usage{
tree(formula = formula(data), data = parent.frame(),
     weights, subset,
     na.action = na.pass, control = tree.control(nobs, ...),
     method = "recursive.partition",
     split = c("deviance", "gini"),
     model = NULL, x = FALSE, y = TRUE, wts = TRUE, ...)
}
\alias{tree}
\alias{print.tree}
\alias{summary.tree}
\alias{print.summary.tree}
\alias{residuals.tree}
\description{
  A tree is grown by binary recursive partitioning using the response in
  the specified formula and choosing splits from the terms of the
  right-hand-side.
}
\arguments{
  \item{formula}{A formula expression. The left-hand-side (response)
    should be either a numerical vector when a regression tree will be
    fitted or a factor, when a classification tree is produced. The
    right-hand-side should be a series of numeric or factor
    variables separated by \code{+}; there should be no interaction
    terms. Both \code{.} and \code{-} are allowed: regression trees can
    have \code{offset} terms.}
  \item{data}{A data frame in which to preferentially interpret
    \code{formula}, \code{weights} and \code{subset}.}
  \item{weights}{Vector of non-negative observational weights; fractional
    weights are allowed.}
  \item{subset}{An expression specifying the subset of cases to be used.}
  \item{na.action}{A function to filter missing data from the model
    frame. The default is \code{na.pass} (to do nothing) as \code{tree}
    handles missing values (by dropping them down the tree as far
    as possible).}
  \item{control}{A list as returned by \code{tree.control}.}
  \item{method}{character string giving the method to use. The only other
    useful value is \code{"model.frame"}.}
  \item{split}{Splitting criterion to use.}
  \item{model}{If this argument is itself a model frame, then the
    \code{formula} and \code{data} arguments are ignored, and
    \code{model} is used to define the model.}
  \item{x}{logical. If true, the matrix of variables for each case
    is returned.}
  \item{y}{logical. If true, the response variable is returned.}
  \item{wts}{logical. If true, the weights are returned.}
  \item{\dots}{Additional arguments that are passed to
    \code{tree.control}. Normally used for \code{mincut}, \code{minsize}
    or \code{mindev}.}
}
\details{
  A tree is grown by binary recursive partitioning using the response in
  the specified formula and choosing splits from the terms of the
  right-hand-side. Numeric variables are divided into
  \eqn{X < a} and \eqn{X > a}; the levels of an unordered factor
  are divided into
  two non-empty groups. The split which maximizes the reduction in
  impurity is chosen, the data set split and the process
  repeated. Splitting continues until the terminal nodes are too small or
  too few to be split.

  Tree growth is limited to a depth of 31 by the use of integers to
  label nodes.

  Factor predictor variables can have up to 32 levels. This limit is
  imposed for ease of labelling, but since their use in a classification
  tree with three or more levels in a response involves a search over
  \eqn{2^{(k-1)}-1}{2^(k-1) - 1} groupings for \eqn{k} levels,
  the practical limit is much less.
}
\value{
  The value is an object of class \code{"tree"} which has components
  \item{frame}{A data frame with a row for each node, and
    \code{row.names} giving the node numbers. The columns include
    \code{var}, the variable used at the split (or \code{"<leaf>"} for a
    terminal node), \code{n}, the (weighted) number of cases reaching
    that node, \code{dev} the deviance of the node, \code{yval}, the
    fitted value at the node (the mean for regression trees, a majority
    class for classification trees) and \code{split}, a two-column
    matrix of the labels for the left and right splits at the
    node. Classification trees also have \code{yprob}, a matrix of
    fitted probabilities for each response level.}
  \item{where}{An integer vector giving the row number of the frame
    detailing the node to which each case is assigned.}
  \item{terms}{The terms of the formula.}
  \item{call}{The matched call to \code{Tree}.}
  \item{x}{If \code{x = TRUE}, the model matrix.}
  \item{y}{If \code{y = TRUE}, the response.}
  \item{wts}{If \code{wts = TRUE}, the weights.}
  and attributes \code{xlevels} and, for classification trees,
  \code{ylevels}.

  A tree with no splits is of class \code{"singlenode"} which inherits
  from class \code{"tree"}.
}
\references{
  Breiman L., Friedman J. H., Olshen R. A., and  Stone,  C. J. (1984)
  \emph{Classification  and Regression Trees.}  Wadsworth.

  Ripley, B. D. (1996)
  \emph{Pattern Recognition and Neural Networks.}
  Cambridge University Press, Cambridge. Chapter 7.
}
\author{B. D. Ripley}


\seealso{\code{\link{tree.control}}, \code{\link{prune.tree}},
   \code{\link{predict.tree}}, \code{\link{snip.tree}}}

\examples{
library(MASS)
data(cpus)
cpus.ltr <- tree(log10(perf) ~ syct+mmin+mmax+cach+chmin+chmax, cpus)
cpus.ltr
summary(cpus.ltr)
plot(cpus.ltr);  text(cpus.ltr)

data(iris)
ir.tr <- tree(Species ~., iris)
ir.tr
summary(ir.tr)
\testonly{ ## tests of singlenode tree
set.seed(1)
Nvars <- 4;
Nsamples <- 5
N <- Nvars*Nsamples;
dat <- matrix(sample(1:2, N, c(0.5,0.5), replace=TRUE), Nsamples, Nvars)
dat <- as.data.frame(dat)
for (i in 1:Nvars) {
   dat[,i]<- factor(dat[,i], 1:2)
}
names(dat) <- c(paste("X", 1:3, sep=""), "Y")
tr <- tree(Y ~ ., dat)
tr
summary(tr)
deviance(tr)
residuals(tr)
}}
\keyword{tree}
