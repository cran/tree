if(!require(ISLR2)) q("no")
library(tree)
train <- subset(Smarket, Year < 2005)   # 998 obs.
f <- Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5
tree(f, data = train, split = "gini", mincut = 10, minsize = 30)
## was ...: maximum depth reached
tree(f, data = train, split = "gini", mincut = 6, minsize = 30)
## was ... tree is too big
