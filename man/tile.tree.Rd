\name{tile.tree}
\title{Add Class Barplots to a Classification Tree Plot}
\usage{
tile.tree(tree, var, screen.arg = ascr + 1, axes = TRUE)
}
\alias{tile.tree}
\arguments{
  \item{tree}{fitted object of class \code{"tree"}.}
  \item{var}{a factor variable to be displayed: by default it is the
    repsonse factor of the tree.}
  \item{screen.arg}{The screen to be used: default the next after the
    currently active screen.}
  \item{axes}{logical flag for drawing of axes for the barplots.}
}
\description{
  The computes a the frequencies of level of \code{var} for cases
  reaching each leaf of the tree, and plots of barchart of the set of
  frequencies underneath each leaf.
}
\value{
  A matrix of counts of categories (rows) for each leaf (columns). The
  principal effect is the plot.
}
\author{B. D. Ripley}

\seealso{\code{\link{tree.screens}}}

\examples{
library(MASS)
data(fgl)
fgl.tr <- tree(type ~ ., fgl)
summary(fgl.tr)
plot(fgl.tr);  text(fgl.tr, all=TRUE, cex=0.5)
fgl.tr1 <- snip.tree(fgl.tr, node=c(108, 31, 26))
tree.screens()
plot(fgl.tr1)
text(fgl.tr1)
tile.tree(fgl.tr1, fgl$type)
close.screen(all = TRUE)
}
\keyword{tree}
\keyword{hplot}
