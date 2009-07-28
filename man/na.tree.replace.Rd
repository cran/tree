% This file is part of the 'tree' package for R
% Copyright (C) B. D. Ripley 1998-2009
% It is distributed under the GPL version 2 or 3

\name{na.tree.replace}
\alias{na.tree.replace}
\title{
  Replace NAs in Predictor Variables
}
\description{
  Adds a new level called \code{"NA"} to any discrete predictor in
  a data frame that contains \code{NA}s.  Stops if any continuous
  predictor contains an \code{NA}.
}
\usage{
na.tree.replace(frame)
}
\arguments{
\item{frame}{
  data frame used to grow a tree.
}}
\value{
  data frame such that a new level named \code{"NA"} is added to
  any discrete predictor in \code{frame} with \code{NA}s.
}
\details{
  This function is used via the \code{na.action} argument to \code{tree}.
}
\seealso{
  \code{\link{tree}}, \code{\link{na.omit}}.
}
\keyword{tree}
