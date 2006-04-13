\name{tree.screens}
\title{Split Screen for Plotting Trees}
\usage{
tree.screens(figs, screen.arg = 0, ...)
}
\alias{tree.screens}
\arguments{
  \item{figs}{A specification of the split of the screen. See
    \code{\link{split.screen}} for the allowed forms.}
  \item{screen.arg}{the screen to divide, by default the whole display area.}
  \item{\dots}{plot parameters to be passed to \code{par}.}
}
\description{
  Splits the screen in a way suitable for using \code{tile.tree}.
}
\value{
  A vector of screen numbers for the newly-created screens.
}
\author{B. D. Ripley}

\seealso{\code{\link{tile.tree}}, \code{\link{split.screen}}}

\examples{
data(fgl, package="MASS")
fgl.tr <- tree(type ~ ., fgl)
summary(fgl.tr)
plot(fgl.tr);  text(fgl.tr, all=TRUE, cex=0.5)
fgl.tr1 <- snip.tree(fgl.tr, node=c(108, 31, 26))
tree.screens()
plot(fgl.tr1)
tile.tree(fgl.tr1, fgl$type)
close.screen(all = TRUE)
}
\keyword{tree}
\keyword{hplot}
