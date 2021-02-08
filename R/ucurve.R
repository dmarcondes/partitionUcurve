#' @import ggplot2
#' @import scales
#' @import ggpubr
#' @export
#' @title U-curve algorith on the Partition Lattice Learning Space
#'
#' @description TBD
#'
#' @details TBD
#'
#' @param x Vetor.
#' @example "TBD"
#' @return TBD

ucurve <- function(xtrain,ytrain,xval,yval){
  #Get parameters
  X <- unique(c(as.character(xtrain),as.character(xval)))
  X <- X[order(X)]
  xtrain <- factor(xtrain,X)
  xval <- factor(xval,X)
  if(length(X) > 20)
    stop("Sorry, but this algorithm is not scalable! It works only for at most 20 points on X domain.")

  Y <- unique(c(as.character(ytrain),as.character(yval)))
  Y <- Y[order(Y)]
  ytrain <- factor(ytrain,Y)
  yval <- factor(yval,Y)
  if(length(Y) != 2)
    stop("The algorithm only work for binary classification problems.")

  suppressWarnings(system("rm visited.dat"))

  n <- length(X)
  search <- 1
  jtrain <- jointDistribution(xtrain,ytrain)
  jval <- jointDistribution(xval,yval)
  part <- makePartition(list(1:n))
  strongMinimums <- vector()

  while(search){
    #Calculate error
    errorPart <- getError(part,jtrain,jval)

    if(is.null(visitedNode(part)))
      addNode(part,errorPart)

    #Get neighbours
    N <- findNeighbors(part)

    #Verbose of new step
    cat("\n")
    cat(paste("At node",namePartition(part),"with error",round(errorPart,5),"and",length(N),"neighbors"))
    cat("\n")

    #Search neighbours for lesser value
    for(i in 1:length(N)){
      Npart <- getPartition(N[i])
      vis <- visitedNode(Npart)
      if(is.null(vis)){
        errorNeighbor <- getError(Npart,jtrain,jval)
        addNode(Npart,errorNeighbor)
      }
      else
        errorNeighbor <- vis
      if(round(errorNeighbor,7) < round(errorPart,7)){
        part <- Npart
        errorPart <- errorNeighbor
        if(length(strongMinimums) > 0){
          minPart <- apply(rbind(strongMinimums),2,getPartition)
          if(sum(unlist(lapply(minPart,function(x) is.related(x,part)))) == 0 & is.null(vis)){
            i <- length(N) + 1
            cat("Found a neighbor with lesser loss. Restarting from it..")
            cat("\n")
          }
          else
            i <- length(N) + 2
        }
        else{
          if(is.null(vis)){
            i <- length(N) + 1
            cat("Found a neighbor with lesser loss. Restarting from it..")
          }
          else
            i <- length(N) + 2
        }
        break
      }
    }

    if(i == length(N)){
      cat("This node is a Strong Local Minimum!! Storing it and starting algorithm again...")
      cat("\n")
      strongMinimums <- c(strongMinimums,namePartition(part))
      part <- findNeighbors(makePartition(list(1:n)))
      part <- part[sample(1:length(part),length(part),F)]
      minPart <- apply(rbind(strongMinimums),2,getPartition)
      for(j in 1:length(part)){
        if(sum(unlist(lapply(minPart,function(x) is.related(x,getPartition(part[j]))))) == 0 & is.null(visitedNode(getPartition(part[j])))){
          part <- getPartition(part[j])
          j <- length(part) + 1
          break
        }
      }
      if(j == length(part)){
        search <- 0
        cat("All Strong Local Minimum have been found!! Finishing algorithm...")
      }
    }
    else if(i == length(N) + 2){
      cat("This node is associated to Strong Local Minimum or is not a Strong Local Minimum. Starting the algorithm again...")
      cat("\n")
      part <- findNeighbors(makePartition(list(1:n)))
      part <- part[sample(1:length(part),length(part),F)]
      minPart <- apply(rbind(strongMinimums),2,getPartition)
      for(j in 1:length(part)){
        if(sum(unlist(lapply(minPart,function(x) is.related(x,getPartition(part[j]))))) == 0 & is.null(visitedNode(getPartition(part[j])))){
          part <- getPartition(part[j])
          j <- length(part) + 1
          break
        }
        else if(is.null(visitedNode(getPartition(part[j]))))
          addNode(getPartition(part[j]),getError(getPartition(part[j]),jtrain,jval))
      }
      if(j == length(part)){
        search <- 0
        cat("All Strong Local Minimum have been found!! Finishing algorithm...")
      }
    }
  }

  #Solution
  e <- apply(rbind(strongMinimums),2,function(x) getError(getPartition(x),jtrain,jval))
  global <- strongMinimums[e == min(e)]
  global <- apply(rbind(global),2,function(x) getPartition(x))
  global <- global[unlist(lapply(global,length)) == min(unlist(lapply(global,length)))]
  solution <- unlist(lapply(global,function(x) namePartition(makePartition(x))))

  #Nodes visite
  visited <- read.table("visited.dat",sep = "e")
  100*nrow(visited)/numbers::bell(n)
}
