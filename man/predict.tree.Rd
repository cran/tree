\name{predict.tree}
\alias{predict.tree}
\title{
Predictions from a Fitted Tree Object
}
\description{
    Returns a vector of predicted responses from a fitted tree object.
}
\usage{
predict.tree(object, newdata = list(),
             type = c("vector", "tree", "class", "where"), 
             split = FALSE, nwts, eps = 1e-3)
}
\arguments{
\item{object}{
    fitted model object of class \code{tree}.  This is assumed to be the result
    of some function that produces an object with the same named
    components as that returned by the \code{tree} function.
}
\item{newdata}{
    data frame containing the values at which predictions are required.
    The predictors referred to in the right side
    of \code{formula(object)} must be present by name in \code{newdata}.
    If missing, fitted values are returned.
}
\item{type}{
    character string denoting whether the predictions are
    returned as a vector (default) or as a tree object.
}
\item{split}{
    governs the handling of missing values. If false, cases with missing
    values are dropped down the tree until a leaf is reached or a node
    for which the attribute is missing, and that node is used for
    prediction. If \code{split = T} cases with missing attributes are
    split into fractional cases and dropped down each side of the split.
    The predicted values are averaged over the fractions to give the
    prediction.
}
\item{nwts}{
    weights for the \code{newdata} cases, used when predicting a tree.
}
\item{eps}{
    a lower bound for the probabilities, used if events of predicted
    probability zero occur in \code{newdata} when predicting a tree.
}}
\value{
if \code{type="vector"}:
vector of predicted responses or, if the response is a factor, matrix
of predicted class probabilities.  This new object is obtained by
dropping \code{newdata} down \code{object}.  For factor predictors, if an
observation contains a level not used to grow the tree, it is left at
the deepest possible node and \code{frame$yval} or \code{frame$yprob} at that
node is the prediction.


if \code{type="tree"}:
an object of class \code{"tree"} is returned with new values
for \code{frame$n} and \code{frame$dev}. If
\code{newdata} does not contain a column for the response in the formula
the value of \code{frame$dev} will be {NA}, and if some values in the
response are missing, the some of the deviances will be \code{NA}.


if \code{type="class"}:
for a classification tree, a factor of the  predicted classes (that
with highest posterior probability, with ties split randomly).

if \code{type="where"}:
the nodes the cases reach.
}
\details{
    This function is a method for the generic function
    \code{predict()} for class \code{tree}.
    It can be invoked by calling \code{predict(x)} for an
    object \code{x} of the appropriate class, or directly by
    calling \code{predict.tree(x)} regardless of the
    class of the object.
}
\references{
Ripley, B.D. (1996).
\emph{Pattern Recognition and Neural Networks.}
Cambridge University Press, Cambridge.
}
\details{
    This function is a method for the generic function
    \code{predict()} for class \code{tree}.
    It can be invoked by calling \code{predict(x)} for an
    object \code{x} of the appropriate class, or directly by
    calling \code{predict.tree(x)} regardless of the
    class of the object.
}
\seealso{
\code{\link{predict}}, \code{\link{tree}}.
}
\examples{
data(shuttle)
shuttle.tr <- tree(use ~ ., shuttle, subset=1:253,
                   mindev=1e-6, minsize=2)
shuttle.tr
shuttle1 <- shuttle[254:256, ]  # 3 missing cases
predict(shuttle.tr, shuttle1)
}
\keyword{tree}
