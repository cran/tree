\name{snip.tree}
\title{Snip Parts of Tree Objects}
\usage{
snip.tree(tree, nodes, xy.save=F, digits=.Options$digits - 3)
}
\alias{snip.tree}
\arguments{
 \item{tree}{An object of class \code{"tree"}.}
 \item{nodes}{An integer vector giving those nodes that are the roots of
   sub-trees to be snipped off. If missing, the user is invited to
   select a node at which to snip.}
 \item{xy.save}{If true, the `x' and `y' coordinates selected
   interactively are saved in the object \code{.xy} in the global
   environment.}
 \item{digits}{Precision used in printing statistics for selected nodes.}
}
\description{
 \code{snip.tree} has two related functions. If \code{nodes} is
 supplied, it removes those nodes and all their descendants from the
 tree.

 If \code{nodes} is not supplied, the user is invited to select
 nodes interactively; this makes sense only if the tree has already been
 plotted. A node is selected by clicking with the left; its number and
 the deviance of the current tree and that which would remain if that
 node were removed are printed. Selecting the same node again causes it
 to be removed.  Clicking the middle button terminates the selection process.
}
\value{
A tree object  containing  the  nodes  that  remain  after
specified or selected subtrees have been snipped off.
}
\author{B.D. Ripley}
\seealso{\code{\link{tree}}, \code{\link{prune.tree}}}

\keyword{tree}
