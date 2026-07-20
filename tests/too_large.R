library(tree)
set.seed(3)
n = 200
data = data.frame(X1 = rnorm(n = n, 0, 1), Y1 = as.factor(rbinom(n, 1, 0.5)))
try(tr <- tree( Y1 ~ X1, data = data, minsize = 25, mindev = 0.001 ))
gc() # fail fail if heap is corrupted

## This estimates nmax = 34, and could corrupt R's heap by trying to make a larger tree.

## nmax = 36 suffices
tr <- tree( Y1 ~ X1, data = data, minsize = 25, mindev = 0.001, nmax=40 )
print(tr)
