#' @import plyr
#' @import rlist
#' @import fst
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
#' @param stop Number of nodes yet to evaluate to trigger exhaustive search.
#' @param path Path to preprocessed partition files.
#' @param Lh An R object with the partition file if it has already been read.
#' @return \item{hypotheses}{The estimated hypothesis of the global minimums with least VC dimension.}
#' @return \item{partitions}{Partitions of the global minimums with least VC dimension.}
#' @return \item{error}{Validation error of the global minimums.}
#' @return \item{exhausted}{Number of nodes exhausted during algorithm.}
#' @return \item{remain}{Number of nodes remaining after algorithm stopped.}
#' @return \item{finished}{If the algorithm was finished or ended after not finding any Strong Local Minimum,}
#' @return \item{SLMvis}{Number of nodes exhasuted until the last Strong Local Minimum was found.}
#' @return \item{remain_after_prune}{Number of nodes remaining after finding each Strong Local Minimum.}
#' @return \item{exhausted_until_prune}{Number of nodes exhausted until finding each Strong Local Minimum.}
#' @examples
#' set.seed(1)
#' x <- sample(x = c("01","02","03","04","05","06","07","08","09","10"),size = 500,replace = T)
#' y <- as.factor(ifelse(as.numeric(x)-5+rnorm(500,0,10) > 0,1,0))
#' x <- factor(x)
#' train <- sample(1:500,350,F)
#' xtrain <- x[train]
#' ytrain <- y[train]
#' xval <- x[!(c(1:500) %in% train)]
#' yval <- y[!(c(1:500) %in% train)]
#' ucurve(xtrain,ytrain,xval,yval)

ucurve <- function(xtrain,ytrain,xval,yval,verbose = T,stop = 0,path = "~/GDrive/Doutorado/Códigos/Particoes/",
                   Lh = NULL){

  #Get sample info
  X <- unique(c(as.character(xtrain),as.character(xval)))
  X <- X[order(X)]
  if(length(X) > 12)
    stop("Sorry, but this algorithm is not scalable! It works only for at most 12 points on X domain.")
  Y <- unique(c(as.character(ytrain),as.character(yval)))
  Y <- Y[order(Y)]
  if(length(Y) != 2)
    stop("The algorithm only work for binary classification problems.")

  #Turn sample into factor
  xtrain <- factor(xtrain,X)
  xval <- factor(xval,X)
  ytrain <- factor(ytrain,Y)
  yval <- factor(yval,Y)

  #Delete file with visited nodes
  suppressWarnings(system("rm visited.dat"))

  #Get parameters
  n <- length(X) #Length of domain
  search <- 1 #Start search
  jtrain <- jointDistribution(xtrain,ytrain) #Train joint distribution
  jval <- jointDistribution(xval,yval) #Validation joint distribution
  part <- makePartition(list(1:n)) #First part to start algorithm
  strongMinimums <- vector() #Declare vector of strong minimums
  vis <- 0 #Start number of exhasted nodes
  SLMvis <- 0 #Exhasuted nodes until last Strong Local Minimum found
  remain_after_prune <- vector() #Nodes remaining after prunig when found each Strong Local Minimum
  exhausted_until_prune <- vector() #Nodes exhausted until each Strong Local Minimum was found

  #Get partitions
  if(is.null(Lh)){
    cat("\n")
    cat("Reading partitions...")
    cat("\n")
    Lh <- fst::read.fst(paste(path,n,"_part.fst",sep = ""))
    nm <- colnames(Lh)
    Lh <- data.frame(Lh)
    colnames(Lh) <- nm
  }
  else{
    if(!is.data.frame(Lh)){
      nm <- colnames(Lh)
      Lh <- data.frame(Lh)
      colnames(Lh) <- nm
    }
  }
  restart <- data.frame("r" = colnames(Lh),"f" = 0) #How many times the algorithm started from each partition

  if(verbose){
    cat("------------------------------------------------------------------------------\n")
    cat("Starting U-curve algorithm...\n")
    cat("------------------------------------------------------------------------------\n")
  }

  #Start algorithm
  while(search){
    #If there is no more nodes to restart from end algorithm
    if(nrow(restart) == 0)
      break

    #If the current partition is not on the space choose another
    if(!(namePartition(part) %in% colnames(Lh))){
      part <- sample(restart$r,1) #Sample new part
      restart$f[restart$r == part] <- restart$f[restart$r == part] + 1 #Add another start to the sampled part
      restart <- restart[restart$f == 0,] #Only parts from which the algorithm has not started
      part <- getPartition(part) #Get partition
    }

    #Calculate error
    errorPart <- getError(part,jtrain,jval)

    #Add to nodes exhasuted
    vis <- vis + 1

    #Save as visited node
    addNode(part,errorPart)

    #Get neighbours
    N <- findNeighbors(part)

    #Error of neighboors
    err <- unlist(lapply(data.frame(rbind(N)),function(x) getError(getPartition(x),jtrain,jval)))

    #If the node is a Strong Local Minimum
    if(errorPart <= min(err)){
      #Nodes exhasuted until found the Strong Local Minimum
      SLMvis <- vis

      #Store node as Strong Local Minimum
      strongMinimums <- c(strongMinimums,namePartition(part))

      #Erase nodes associated to it
      part_num <- Lh[,colnames(Lh) == namePartition(part)]
      Lh <- Lh[,unlist(lapply(Lh,function(x) !is.related(x,part_num)))]
      restart <- restart[restart$r %in% colnames(Lh),]

      if(verbose){
        cat("\n")
        cat("Found a SLM. There are ",nrow(restart)," (",round(100*nrow(restart)/numbers::bell(n),2),"%) nodes remaining after ",
            vis," (",round(100*vis/numbers::bell(n),2),"%) exhausted...",sep = "")
        cat("\n")
      }

      #Nodes remaining after pruning
      remain_after_prune <- c(remain_after_prune,nrow(restart))

      #Nodes exhasuted until found Strong Local Minimum
      exhausted_until_prune <- c(exhausted_until_prune,vis)

      if(nrow(restart) <= stop){
        #If less nodes remaining than stop, stop algorithm
        search <- 0
      }
      else{
        #Choose a new node to restart
        part <- sample(restart$r,1) #Sample new part
        restart$f[restart$r == part] <- restart$f[restart$r == part] + 1 #Add another start to the sampled part
        restart <- restart[restart$f == 0,] #Only parts from which the algorithm has not started
        part <- getPartition(part) #Get partition
      }
    }
    else{
      #Erase last node
      if(!is.null(ncol(Lh)))
        Lh <- Lh[,colnames(Lh) != namePartition(part)]

      #Neighboors with lesser loss still on space
      N1 <-  N[err <= errorPart] #Get neighboors with less loss
      N1 <- N1[N1 %in% colnames(Lh)] #Getneighboors still on the space
      if(length(N1) > 0)
        part <- getPartition(sample(N1,1)) #If some neighboor is available to continue, sample next node from it
      else{
        if(nrow(restart) > 0){
          #Choose a new node to restart
          part <- sample(restart$r,1) #Sample new part
          restart$f[restart$r == part] <- restart$f[restart$r == part] + 1 #Add another start to the sampled part
          restart <- restart[restart$f == 0,] #Only parts from which the algorithm has not started
          part <- getPartition(part) #Get partition
        }
        else
          break #End algorithm
      }
    }
    #After 1000 points have been exhausted, end algorithm if it takes more than twice the exhasuted nodes to find next Strong Local Minimum
    if(SLMvis > 1000)
      if(vis >= 2*SLMvis)
        break
  }

  remain <- nrow(restart) #Number of nodes remaining

  #If nodes remaining start exhaustion
  if(remain > 0 & remain <= stop){
    finished <- 1
    if(verbose){
      cat("\n")
      cat(paste("Ending U-curve algorithm and starting exhaustion of ",nrow(restart)," (",round(100*nrow(restart)/numbers::bell(n),2),
                "%) nodes remaining...",sep = ""))
      cat("\n")
    }

    #Get error of all remaining
    err <- unlist(lapply(data.frame(rbind(restart$r)),function(x) getError(getPartition(x),jtrain,jval)))

    #Get minimum error of strong local minimums
    m <- min(unlist(lapply(data.frame(rbind(strongMinimums)),function(x) getError(getPartition(x),jtrain,jval))))

    #Add nodes with less error to strong local minimums
    strongMinimums <- c(strongMinimums,colnames(Lh)[err <= m])

    if(verbose){
      cat("\n")
      cat("---------------------------------------------------------------------------------")
      cat("\n")
      cat(paste("There was ",remain," (",round(100*remain/numbers::bell(n),2),"%)"," nodes to visit after exhausting ",
                vis," (",round(100*vis/numbers::bell(n),2),"%)",sep = ""))
      cat("\n")
      cat("---------------------------------------------------------------------------------")
    }
  }
  else if(remain > stop){
    finished <- 0
    if(verbose){
      cat("\n")
      cat(paste("Ending U-curve algorithm with ",nrow(restart)," (",round(100*nrow(restart)/numbers::bell(n),2),
                "%) nodes remaining. No Strong Local Minimum was found after the exhaustion of ",SLMvis,
                " nodes until the exhaustion of ",vis," nodes...",sep = ""))
      cat("\n")
    }
  }
  else if(remain == 0){
    finished <- 1
    if(verbose){
      cat("\n")
      cat("---------------------------------------------------------------------------------")
      cat("\n")
      cat(paste("There was ",remain," (",round(100*remain/numbers::bell(n),2),"%)"," nodes to visit after exhausting ",
                vis," (",round(100*vis/numbers::bell(n),2),"%)",sep = ""))
      cat("\n")
      cat("---------------------------------------------------------------------------------")
    }
  }

  #Solution
  e <- apply(rbind(strongMinimums),2,function(x) getError(getPartition(x),jtrain,jval))
  global <- strongMinimums[e == min(e)]
  global <- apply(rbind(global),2,function(x) getPartition(x))
  global <- global[unlist(lapply(global,length)) == min(unlist(lapply(global,length)))]
  solution <- unique(unlist(lapply(global,function(x) namePartition(makePartition(x)))))

  #Solution
  oh <- tapply(cbind(solution),solution,function(x) optimalHyp(part = x,jtrain = jtrain))
  oh <- lapply(oh,function(x) data.frame("X" = X,"fx" = plyr::mapvalues(factor(x[,2]),c("1","0"),Y)))
  error <- unique(tapply(cbind(solution),solution,function(x) getError(part = getPartition(x),jtrain = jtrain,jval = jval)))
  names(error) <- NULL

  return(list("hypotheses" = oh,"partitions" = solution,"error" = error,
              "exhausted" = vis,"remain" = remain,"finished" = finished,"SLMvis" = ifelse(!finished,SLMvis,NA),
              "remain_after_prune" = remain_after_prune,"exhausted_until_prune" = exhausted_until_prune))
}
