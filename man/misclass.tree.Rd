\name{misclass.tree}
\title{Report Numbers of Misclassifications by a Classification Tree}
\usage{
misclass.tree(tree, detail=F)
}
\alias{misclass.tree}
\arguments{
 \item{tree}{Object of class \code{"tree"}, representing a
   classification tree.}
 \item{detail}{If false, report overall number of
   mis-classifications. If true, report the number of each node.}
}
\description{
  Report the number of mis-classifications made by a classification
  tree, either overall or at each node.
}
\details{
 The quantities returned are weighted by the observational weights if
 these are supplied in the construction of \code{tree}.
}
\value{
  Either the overall number of misclassifications or the number for each node.
}
\author{B.D. Ripley}
\seealso{\code{\link{tree}}}

\examples{
data(iris)
ir.tr <- tree(Species ~., iris)
misclass.tree(ir.tr)
misclass.tree(ir.tr, detail=T)
}
\keyword{tree}
