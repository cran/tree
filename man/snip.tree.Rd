% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2013
% It is distributed under the GPL version 2 or 3

\name{snip.tree}
\title{Snip Parts of Tree Objects}
\usage{
snip.tree(tree, nodes, xy.save = FALSE,
          digits = getOption("digits") - 3)
}
\alias{snip.tree}
\arguments{
  \item{tree}{An object of class \code{"tree"}.}
  \item{nodes}{An integer vector giving those nodes that are the roots of
    sub-trees to be snipped off. If missing, the user is invited to
    select a node at which to snip.}
  \item{xy.save}{If true, the \code{x} and \code{y} coordinates selected
    interactively are saved as attribute \code{.xy} of the returned value.}
  \item{digits}{Precision used in printing statistics for selected nodes.}
}
\description{
  \code{snip.tree} has two related functions. If \code{nodes} is
  supplied, it removes those nodes and all their descendants from the
  tree.

  If \code{nodes} is not supplied, the user is invited to select nodes
  interactively; this makes sense only if the tree has already been
  plotted. A node is selected by clicking with the left mouse button;
  its number and the deviance of the current tree and that which would
  remain if that node were removed are printed. Selecting the same node
  again causes it to be removed (and the lines of its sub-tree erased).
  Clicking any other button terminates the selection process.
}
\value{
  A tree object containing the nodes that remain after
  specified or selected subtrees have been snipped off.
}
\note{
  Prior to version 1.0-34, the saved coordinates were place in object
  \code{.xy} in the workspace.
}

\author{B. D. Ripley}
\seealso{\code{\link{tree}}, \code{\link{prune.tree}}.}

\keyword{tree}
