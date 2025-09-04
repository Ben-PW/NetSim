#####################################################################################################

#Data preprocessing script. Import, pre-process, fit ERGMs and check GoFs here. Create objects to be 
#passed into Combine3.R

#####################################################################################################

library(here)
here()
library(statnet)
library(networkdata)

########################################## Florentine Families ######################################
# 16x16 undirected, low density, isolates

data(flo)

flo <- network(flo, directed = F)

flo %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)

flo1 <- ergm(flo ~ edges + absdiff("wealth"))
rm(flo)

###############################  ############################
# 4 x 24 x 24 weighted adjacency matrices

data(acct)

ac1 <- network(acct[[1]], directed = F)
ac2 <- network(acct[[2]], directed = F)

plot(ac1)
plot(ac2)

ac1 %v% "gender"

acct1 <- ergm(ac1 ~ edges + nodematch("gender") + nodematch("job"))
summary(acct1)

acct2 <- ergm(ac2 ~ edges + nodematch("gender") + nodematch("job"))
summary(acct2)

rm(ac1, ac2, acct)

