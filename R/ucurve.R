#' @import plyr
#' @import rlist
#' @export
#' @title U-curve algorith on the Partition Lattice Learning Space
#'
#' @description TBD
#'
#' @details TBD
#'
#' @param xtrain Vector with the training sample of x values.
#' @param ytrain Vector with the training sample of y values.
#' @param xval Vector with the validation sample of x values.
#' @param yval Vector with the validation sample of y values.
#' @param verbose Logical to print a trace of the algorithm.
#' @return \item{hypotheses}{The estimated hypothesis of the global minimums with least VC dimension.}
#' @return \item{partitions}{Partitions of the global minimums with least VC dimension.}
#' @return \item{error}{Error of the global minimums.}
#' @return \item{visited}{Number and percentage of nodes visited.}
#' @examples
#' set.seed(1)
#' x <- sample(x = 1:10,size = 50,replace = T)
#' y <- as.factor(ifelse(x-5+rnorm(50,0,10) > 0,1,0))
#' x <- factor(x)
#' train <- sample(1:50,35,F)
#' xtrain <- x[train]
#' ytrain <- y[train]
#' xval <- x[!(c(1:50) %in% train)]
#' yval <- y[!(c(1:50) %in% train)]
#' ucurve(xtrain,ytrain,xval,yval)

ucurve <- function(xtrain,ytrain,xval,yval,verbose = F){

  #Get sample info
  X <- unique(c(as.character(xtrain),as.character(xval)))
  X <- X[order(X)]
  #if(length(X) > 20)
  #  stop("Sorry, but this algorithm is not scalable! It works only for at most 20 points on X domain.")
  xtrain <- factor(xtrain,X)
  xval <- factor(xval,X)

  Y <- unique(c(as.character(ytrain),as.character(yval)))
  Y <- Y[order(Y)]
  if(length(Y) != 2)
    stop("The algorithm only work for binary classification problems.")
  ytrain <- factor(ytrain,Y)
  yval <- factor(yval,Y)

  #Delete file with visited nodes
  suppressWarnings(system("rm visited.dat"))

  #Get parameters
  n <- length(X)
  search <- 1
  jtrain <- jointDistribution(xtrain,ytrain)
  jval <- jointDistribution(xval,yval)
  part <- makePartition(list(1:n))
  errorPart <- getError(part,jtrain,jval)
  addNode(part,errorPart)
  restart <- findNeighbors(part)
  strongMinimums <- vector()

  while(search){
    #Calculate error
    errorPart <- getError(part,jtrain,jval)

    #Save as visited node
    addNode(part,errorPart)

    #Get neighbours
    N <- findNeighbors(part)

    #Verbose of new step
    if(verbose){
      cat("\n")
      cat(paste("At node",namePartition(part),"with error",round(errorPart,5),"and",length(N),"neighbors"))
      cat("\n")
    }

    #Error of neighboors
    err <- unlist(plyr::alply(rbind(N),2,function(x) getError(getPartition(x),jtrain,jval)))
    if(errorPart <= min(err)){
      if(verbose){
        cat("This node is a Strong Local Minimum!! Storing it and starting algorithm again...")
        cat("\n")
      }
      #Store node as Strong Local Minimum
      strongMinimums <- c(strongMinimums,namePartition(part))

      #Find a new node at second floor to restart algorithm
      restart <- restart[!unlist(lapply(plyr::alply(rbind(restart),2,getPartition),function(x) is.related(x,part)))]

      if(length(restart) > 0){
        #Choose a new node to restart
        part <- getPartition(restart[sample(1:length(restart),1)])
      }
      else{
        #If all the second floor is related to a Strong Local Minimum we are done
        search <- 0
        if(verbose){
          cat("All Strong Local Minimum have been found!! Finishing algorithm...")
        }
      }
    }
    else{
      N <-  N[err < errorPart]
      if(length(strongMinimums) > 0){
        #Get partition of all Strong Local Minimums found so far
        minPart <- plyr::alply(rbind(strongMinimums),2,getPartition)

        #Only neighboors not related to SLM
        r <- unlist(plyr::alply(rbind(N),2,function(y) sum(unlist(lapply(minPart,function(x) is.related(x,y)))) == 0))
      }
      else
        r <- rep(TRUE,length(N))
      N <- N[r & unlist(lapply(plyr::alply(rbind(N),2,getPartition),function(x) length(x) > length(part)))]
      if(length(N) > 0){
        if(verbose){
          cat("Found a greater neighbor with lesser loss. Restarting from it..")
          cat("\n")
        }

        #Neighbor is the new node
        part <- getPartition(sample(N,1))
      }
      else{
        #When the neighboor with less loss is related to a Strong Local Minimum
        if(verbose){
          cat("All neighboors with lesser loss are associated to a Strong Local Minimum or are lesser. Restarting...")
          cat("\n")
        }
        part <- getPartition(restart[sample(1:length(restart),1)])
      }
    }
  }

  #Solution
  e <- apply(rbind(strongMinimums),2,function(x) getError(getPartition(x),jtrain,jval))
  global <- strongMinimums[e == min(e)]
  global <- apply(rbind(global),2,function(x) getPartition(x))
  global <- global[unlist(lapply(global,length)) == min(unlist(lapply(global,length)))]
  solution <- unlist(lapply(global,function(x) namePartition(makePartition(x))))

  #Nodes visited
  visited <- read.table("visited.dat",sep = "e")
  if(verbose){
    cat("\n")
    cat("---------------------------------------------------------------------------------")
    cat("\n")
    cat(paste("Nodes visited: ",nrow(visited)," (",round(100*nrow(visited)/numbers::bell(n),2),"%)",sep = ""))
    cat("\n")
    cat("---------------------------------------------------------------------------------")
  }

  #Solution
  oh <- tapply(cbind(solution),solution,function(x) optimalHyp(part = x,jtrain = jtrain))
  oh <- lapply(oh,function(x) data.frame("X" = X,"fx" = plyr::mapvalues(factor(x[,2]),c("1","0"),Y)))
  error <- unique(tapply(cbind(solution),solution,function(x) getError(part = getPartition(x),jtrain = jtrain,jval = jval)))
  names(error) <- NULL

  return(list("hypotheses" = oh,"partitions" = solution,"error" = error,
              "visited" = paste(nrow(visited)," (",round(100*nrow(visited)/numbers::bell(n),2),"%)",sep = "")))
}
