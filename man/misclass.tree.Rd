% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2007
% It is distributed under the GPL version 2 or 3

\name{misclass.tree}
\title{Misclassifications by a Classification Tree}
\usage{
misclass.tree(tree, detail = FALSE)
}
\alias{misclass.tree}
\arguments{
  \item{tree}{Object of class \code{"tree"}, representing a
    classification tree.}
  \item{detail}{If false, report overall number of
    mis-classifications. If true, report the number at each node.}
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
\author{B. D. Ripley}
\seealso{\code{\link{tree}}}

\examples{
ir.tr <- tree(Species ~., iris)
misclass.tree(ir.tr)
misclass.tree(ir.tr, detail=TRUE)
}
\keyword{tree}
