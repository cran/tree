\name{text.tree}
\alias{text.tree}
\title{Annotate a Tree Plot}
\description{
  Add text to a tree plot.
}
\usage{
\method{text}{tree}(x, splits = TRUE, label = "yval", all = FALSE,
     pretty = NULL, digits = getOption("digits") - 3, 
     adj = par("adj"), xpd = TRUE, ...)
}
\arguments{
 \item{x}{an object of class \code{"tree"}}
 \item{splits}{logical. If \code{TRUE} the splits are labelled}
 \item{label}{The name of column in the \code{frame} component of
   \code{x}, to be used to label the nodes. Can be \code{NULL} to
   suppress node-labelling}
 \item{all}{logical. By default, only the leaves are labelled, but if
   true interior nodes are also labelled.}
 \item{pretty}{the manipulation used for split labels involving
   attributes. See Details.}
 \item{digits}{significant digits for numerical labels.}
 \item{adj, xpd, \dots}{graphical parameters such as \code{cex} and \code{font}.}
}
\details{
  If \code{pretty = 0} then the level names of a factor split attributes
  are used unchanged.  If \code{pretty = NULL}, the levels are presented
  by \code{a}, \code{b}, \dots.  If \code{pretty} is a positive integer,
  \code{\link{abbreviate}} is applied to the labels with that value
  for its argument \code{minlength}.

  If the lettering is vertical (par \code{srt = 90}) and \code{adj} is
  not supplied it is adjusted appropriately.
}
\value{
  None.
}
\author{B. D. Ripley}

\seealso{\code{\link{plot.tree}}}
\examples{
ir.tr <- tree(Species ~., iris)
plot(ir.tr)
text(ir.tr)
}
\keyword{tree}
\keyword{hplot}
