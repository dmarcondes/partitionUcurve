source("R/utils.R")

#Partitions to test algorithm
n <- 10
part1 <- makePartition(list(c(1,4,5),c(3,6),c(2,7,9),c(8,10)))
part2 <- makePartition(list(c(1),c(4,5),c(3,6),c(2),c(7,9),c(8,10)))
part3 <- makePartition(list(c(1,4,5,8),c(3,6),c(2,7,9,10)))

is.related(part1,part2)
is.related(part2,part1)
is.related(part1,part1)
is.related(part1,part3)
is.related(part3,part2)

#Data
set.seed(1)
x <- sample(x = 1:10,size = 50,replace = T)
y <- as.factor(ifelse(x-5+rnorm(50,0,10) > 0,1,0))
x <- factor(x)
train <- sample(1:50,35,F)
xtrain <- x[train]
ytrain <- y[train]
xval <- x[!(c(1:50) %in% train)]
yval <- y[!(c(1:50) %in% train)]
prop.table(table(xtrain,ytrain))
prop.table(table(xval,yval))

jtrain <- jointDistribution(xtrain,ytrain)
jval <- jointDistribution(xval,yval)

getError(part3,jtrain,jval)
