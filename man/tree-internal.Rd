\name{tree-internal}
\alias{data.tree}
\alias{descendants}
\alias{labels.tree}
\alias{model.frame.tree}
\alias{node.match}
\alias{pred1.tree}
\alias{tree.depth}
\alias{tree.matrix}
\alias{treeco}
\alias{treepl}
\title{Internal Functions Used by the `tree' Package}
\description{
  Internal functions used in the \code{tree} package.
}
\usage{
data.tree(tree)
descendants(nodes, include=TRUE)
labels.tree(object, pretty = TRUE, collapse = TRUE, ...)
model.frame.tree(formula, ...)
node.match(nodes, nodelist, leaves, print.it = TRUE)
pred1.tree(tree, x)
tree.depth(nodes)
tree.matrix(frame)
treeco(tree, uniform = paste(".Tree.unif", dev.cur(), sep = "."))
treepl(xy, node, erase = FALSE, ...)
}
\note{These are not intended for use by end users.}
\keyword{utilities}
