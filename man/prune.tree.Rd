\name{prune.tree}
\alias{prune.tree}
\alias{prune.misclass}
\title{
Cost-complexity Pruning of Tree Object
}
\description{
  Determines a nested sequence of subtrees of the supplied tree by
  recursively "snipping" off the least important splits.
}
\usage{
prune.tree(tree, k = NULL, best = NULL, newdata, nwts,
        method = c("deviance", "misclass"), loss, eps = 1e-3)
prune.misclass(tree, k = NULL, best = NULL, newdata, 
        nwts, loss, eps = 1e-3)
}
\arguments{
  \item{tree}{
    fitted model object of class \code{tree}.  This is assumed to be the result
    of some function that produces an object with the same named
    components as that returned by the \code{tree()} function.
  }
  \item{k}{
    cost-complexity parameter defining either a specific subtree of \code{tree}
    (\code{k} a scalar) or the (optional) sequence of subtrees minimizing the
    cost-complexity measure (\code{k} a vector). If missing, \code{k} is
    determined algorithmically.
  }
  \item{best}{
    integer requesting the size (i.e. number of terminal nodes) of a
    specific subtree in the cost-complexity sequence to be returned. This
    is an alternative way to select a subtree than by supplying a scalar
    cost-complexity parameter \code{k}.  If there is no tree in the sequence of
    the requested size, the next largest is returned.
  }
  \item{newdata}{
    data frame upon which the sequence of cost-complexity subtrees is
    evaluated.  If missing, the data used to grow the tree are used.
  }
  \item{nwts}{
    weights for the \code{newdata} cases.
  }
  \item{method}{
    character string denoting the measure of node heterogeneity used to
    guide cost-complexity pruning.  For regression trees, only the
    default, \code{deviance}, is accepted.  For classification trees, the
    default is \code{deviance} and the alternative is \code{misclass}
    (number of misclassifications or total loss).
  }
  \item{loss}{
    a matrix giving for each true class (row) the numeric loss of
    predicting the class (column).  The classes should be in the order of
    the levels of the response.  It is conventional for a loss matrix to
    have a zero diagonal.
  }
  \item{eps}{
    a lower bound for the probabilities, used to compute deviances if
    events of predicted probability zero occur in \code{newdata}.
  }
}
\details{
  Determines a nested sequence of subtrees of the supplied tree by
  recursively "snipping" off the least important splits, based upon
  the cost-complexity measure.  \code{prune.misclass} is an abbreviation for
  \code{prune.tree(method = "misclass")} for use with \code{cv.tree}.
  
  If \code{k} is supplied, the optimal subtree for that value is returned.

  The response as well as the predictors referred to in the right side
  of the formula in \code{tree} must be present by name in
  \code{newdata}. These data are dropped down each tree in the
  cost-complexity sequence and deviances or losses calculated by
  comparing the supplied response to the prediction.  The function
  \code{\link{cv.tree}()} routinely uses the \code{newdata} argument
  in cross-validating the pruning procedure.  A \code{plot} method
  exists for objects of this class.  It displays the value of the
  deviance, the number of misclassifications or the total loss for
  each subtree in the cost-complexity sequence.  An additional axis
  displays the values of the cost-complexity parameter at each subtree.
}
\value{
  If \code{k} is supplied and is a scalar, a \code{tree} object is
  returned that minimizes the cost-complexity measure for that \code{k}.
  If \code{best} is supplied, a \code{tree} object of size \code{best}
  is returned.  Otherwise, an object of class \code{tree.sequence}
  is returned.  The object contains the following components:

  \item{size}{
    number of terminal nodes in each tree in the cost-complexity pruning
    sequence.
  }
  \item{deviance}{
    total deviance of each tree in the cost-complexity pruning sequence.
  }
  \item{k}{
    the value of the cost-complexity pruning parameter of each tree in 
    the sequence.
  }
}
\examples{
library(MASS)
data(fgl)
fgl.tr <- tree(type ~ ., fgl)
plot(print(fgl.tr))
fgl.cv <- cv.tree(fgl.tr,, prune.tree)
for(i in 2:5)  fgl.cv$dev <- fgl.cv$dev +
   cv.tree(fgl.tr,, prune.tree)$dev
fgl.cv$dev <- fgl.cv$dev/5
plot(fgl.cv)
}
\keyword{tree}
