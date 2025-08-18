# from Duncan Murdoch 2025-07-27
if(!require(ISLR2)) q("no")
options(digits = 3) # too small
library(tree)
High <- factor(ifelse(Carseats$Sales <= 8, "No", "Yes"))
Carseats <- data.frame(Carseats, High)
tree(High ~ . - Sales, Carseats)
