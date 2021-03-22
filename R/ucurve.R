#' @import plyr
#' @import rlist
#' @import fst
#' @import parallel
#' @import numbers
#' @import stringi
#' @import ggplot2
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
#' @param optimal Logical indicating if the algorithm should return an optimal solution.
#' @param exhaust Number of points to exhaust before stopping the algorithm.
#' @param sampleNeigh Either false to consider all neighboors, or the maximum number of neighboors to sample at each exhaustion. If a number,then optimal should be false.
#' @param verbose Logical to print a trace of the algorithm.
#' @param stop Number of nodes yet to evaluate to trigger exhaustive search.
#' @param path Path to preprocessed partition files.
#' @param Lh A data frame with the partition lattice.
#' @param cores Number of cores for parallel computing.
#' @return \item{hypotheses}{The estimated hypothesis of the global minimums with least VC dimension.}
#' @return \item{partitions}{Partitions of the global minimums with least VC dimension.}
#' @return \item{error}{Validation error of the global minimums.}
#' @return \item{exhausted}{Number of nodes exhausted during algorithm.}
#' @return \item{remain}{Number of nodes remaining after algorithm stopped.}
#' @return \item{finished}{If the algorithm was finished or ended after not finding any Strong Local Minimum.}
#' @return \item{SLMvis}{Number of nodes exhasuted until the last Strong Local Minimum was found.}
#' @return \item{remain_after_prune}{Number of nodes remaining after finding each Strong Local Minimum.}
#' @return \item{exhausted_until_prune}{Number of nodes exhausted until finding each Strong Local Minimum.}
#' @return \item{optimal}{Wheter an optimal solution was returned.}
#' @examples
#' set.seed(1)
#' x <- sample(x = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20"),size = 500,replace = T)
#' y <- as.factor(ifelse(as.numeric(x)-10+rnorm(500,0,20/3) > 0,1,0))
#' x <- factor(x)
#' train <- sample(1:500,350,F)
#' xtrain <- x[train]
#' ytrain <- y[train]
#' xval <- x[!(c(1:500) %in% train)]
#' yval <- y[!(c(1:500) %in% train)]
#' u <- ucurve(xtrain,ytrain,xval,yval,optimal = F,sampleNeigh = 5000)

ucurve <- function(xtrain,ytrain,xval,yval,optimal = T,exhaust = 1000,sampleNeigh = F,
                   verbose = T,stop = 0,path = "~/GDrive/Doutorado/Códigos/Particoes/",Lh = NULL,cores = 4){

  #Test if parameters correct
  if(!!sampleNeigh){
    if(optimal)
      stop("Can only sample neighboors for suboptimal algorithm...")
  }

  #Timing
  start_time <- Sys.time()

  #Get sample info
  X <- unique(c(as.character(xtrain),as.character(xval)))
  X <- X[order(X)]
  n <- length(X) #Length of domain
  if(optimal)
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

  if(verbose){
    cat("------------------------------------------------------------------------------\n")
    cat(paste(length(unique(X)),"points on the domain\n"))
    cat(paste(ifelse(n <= 218,format(bell(n),big.mark = ",",scientific = F),"???"),
                            "nodes on the Partition Lattice Learning Space\n"))
    cat(paste("Training sample size:",length(xtrain),"\n"))
    cat(paste("Validation sample size:",length(xval),"\n"))
    cat(paste("Returning",ifelse(optimal,"optimal","suboptimal"),"hypotheses\n"))
    cat("------------------------------------------------------------------------------\n")
    cat("\n")
  }

  #Delete file with visited nodes
  suppressWarnings(system("rm visited.dat"))

  #Get parameters
  search <- 1 #Start search
  jtrain <- jointDistribution(xtrain,ytrain) #Train joint distribution
  jval <- jointDistribution(xval,yval) #Validation joint distribution
  part <- makePartition(list(1:n)) #First part to start algorithm
  strongMinimums <- vector() #Declare vector of strong minimums
  vis <- 0 #Start number of exhasted nodes
  SLMvis <- 0 #Exhasuted nodes until last Strong Local Minimum found
  exhausted_until_prune <- vector() #Nodes exhausted until each Strong Local Minimum was found
  min_error <- NULL
  if(optimal){
    remain_after_prune <- vector() #Nodes remaining after prunig when found each Strong Local Minimum
  }

  #Get partitions
  if(optimal){
    if(is.null(Lh)){
      cat("\n")
      cat("Reading partitions...")
      cat("\n")
      Lh <- read.fst(paste(path,"part_",n,".fst",sep = ""))
    }
  }

  if(verbose){
    cat("------------------------------------------------------------------------------\n")
    cat("Starting U-curve algorithm...\n")
    cat("------------------------------------------------------------------------------\n")
  }

  if(verbose)
    cat("*")

  #Start algorithm
  while(search){
    if(optimal){
      #If there is no more nodes to restart from end algorithm
      if(length(colnames(Lh)) == 0)
        break
    }

    #Calculate error
    errorPart <- getError(part,jtrain,jval)

    #Add to nodes exhasuted
    vis <- vis + 1

    #Save as visited node
    addNode(part,errorPart)

    #Get neighbours
    if(!sampleNeigh)
      N <- findNeighbors(part,n,cores)
    else
      N <- as.vector(sampleNeighbors(sample_size = sampleNeigh,part = part,n = n,cores = cores))

    #Error of neighboors
    err <- unlist(mclapply(data.frame(rbind(N)),function(x) getError(getPartition(x),jtrain,jval),mc.cores = cores))
    if(optimal)
      Lh <- Lh[,!(colnames(Lh) %in% N[err > errorPart])]

    #If the node is a Strong Local Minimum
    if(errorPart <= min(err)){
      #Plot
      if(is.null(min_error))
        min_error <- errorPart
      else
        min_error <- c(min_error,errorPart)#ifelse(errorPart < min(min_error),errorPart,min(min_error)))
      p <- data.frame("x" = 1:(length(strongMinimums) + 1),y = min_error)
      if(length(strongMinimums) > 0){
        p <- ggplot(p,aes(x = x,y = y,group = 1)) + theme_linedraw() + geom_point(color = "salmon") + geom_line(color = "salmon") +
          xlab("Strong Local Minimum") + ylab("Validation Error") +
          scale_x_continuous(breaks = 1:1e4)
        print(p)
      }


      #Nodes exhasuted until found the Strong Local Minimum
      SLMvis <- vis

      #Store node as Strong Local Minimum
      strongMinimums <- c(strongMinimums,namePartition(part))

      if(optimal){
        #Erase nodes associated to it
        num_part <- Lh[,colnames(Lh) == namePartition(part)]
        tmp <- Lh[,!unlist(mclapply(Lh,function(x) is.related(x,num_part),mc.cores = cores))]
      }

      if(verbose){
        if(optimal){
          cat("\n")
          cat("Found a SLM. There are ",ncol(Lh)," (",ifelse(n <= 218,round(100*ncol(Lh)/bell(n),5),0),"%) nodes remaining after ",
              vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%) exhausted",sep = "")
          cat("\n")
        }
        else{
          cat("\n")
          cat("Found a SLM after exhausting ",vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = "")
          cat("\n")
        }
      }

      if(optimal){
        #Nodes remaining after pruning
        remain_after_prune <- c(remain_after_prune,length(colnames(Lh)))
      }

      #Nodes exhasuted until found Strong Local Minimum
      exhausted_until_prune <- c(exhausted_until_prune,vis)

      if(optimal){
        if(length(colnames(Lh)) <= stop){
          #If less nodes remaining than stop, stop algorithm
          search <- 0
        }
        else{
          #Choose a new node to restart
          part <- sample(colnames(Lh),1) #Sample new part
          part <- getPartition(part) #Get partition
        }
      }
      else{
        #Sample new partition to restart
        part <- samplePartition(n)
      }
      #If exhausted more points than desired
      if(!optimal){
        if(vis > exhaust)
          break
      }
      if(verbose)
        cat("*")
    }
    else{
      N1 <-  N[err <= errorPart] #Get neighboors with less loss

      if(optimal){
        #Erase last node
        Lh <- Lh[,colnames(Lh) != namePartition(part)]

        #Neighboors with lesser loss still on space
        N1 <- N1[N1 %in% colnames(Lh)] #Get neighboors still on the space
      }

      #If some neighboor is available to continue, sample next node from it
      if(length(N1) > 0){
          part <- getPartition(sample(N1,1))
          if(verbose)
            cat("-")
      }
      else{
        if(optimal){
          if(length(colnames(Lh)) > 0){
            #Choose a new node to restart
            part <- sample(colnames(Lh),1) #Sample new part
            part <- getPartition(part) #Get partition
          }
          else
            break #End algorithm
        }
        else{
          #Sample new partition to restart
          part <- samplePartition(n)
        }

        #If exhausted more points than desired
        if(!optimal){
          if(vis > exhausted)
            break
        }
        if(verbose)
          cat("*")
      }
    }
  }

  if(optimal){
    remain <- length(colnames(Lh)) #Number of nodes remaining

  #If nodes remaining start exhaustion
    if(remain > 0 & remain <= stop){
      finished <- 1
      if(verbose){
        cat("\n")
        cat(paste("Ending U-curve algorithm and starting exhaustion of ",ncol(Lh)," (",
                  ifelse(n <= 218,round(100*ncol(Lh)/bell(n),5),0),
                  "%) nodes remaining...",sep = ""))
        cat("\n")
      }

      #Get error of all remaining
      err <- unlist(mclapply(data.frame(rbind(colnames(Lh))),function(x) getError(getPartition(x),jtrain,jval),mc.cores = cores))

      #Get minimum error of strong local minimums
      m <- min(unlist(mclapply(data.frame(rbind(strongMinimums)),function(x) getError(getPartition(x),jtrain,jval),mc.cores = cores)))

      #Add nodes with less error to strong local minimums
      strongMinimums <- c(strongMinimums,colnames(Lh)[err <= m])

      if(verbose){
        cat("\n")
        cat("---------------------------------------------------------------------------------")
        cat("\n")
        cat(paste("There was ",remain," (",ifelse(n <= 218,round(100*remain/bell(n),5),0),"%)"," nodes to visit after exhausting ",
                  vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = ""))
        cat("\n")
        cat("---------------------------------------------------------------------------------")
      }
    }
    else if(remain == 0){
      finished <- 1
      if(verbose){
        cat("\n")
        cat("---------------------------------------------------------------------------------")
        cat("\n")
        cat(paste("There was ",remain," (",ifelse(n <= 218,round(100*remain/bell(n),5),0),"%)"," nodes to visit after exhausting ",
                  vis," (",ifelse(n <= 218,round(100*vis/bell(n),5),0),"%)",sep = ""))
        cat("\n")
        cat("---------------------------------------------------------------------------------")
      }
    }
  }
  else{
    finished <- 1
    if(verbose){
      cat("\n")
      cat(paste("U-curve algorithm returned a suboptimal solution after the exhaustion of ",vis," (",
                ifelse(n <= 218,round(100*vis/bell(n),5),0),
                "%) nodes...",sep = ""))
      cat("\n")
    }
  }

  #Solution
  e <- apply(rbind(strongMinimums),2,function(x) getError(getPartition(x),jtrain,jval))
  global <- strongMinimums[e == min(e)]
  global <- apply(rbind(global),2,function(x) getPartition(x))
  #global <- global[unlist(lapply(global,length)) == min(unlist(lapply(global,length)))]
  solution <- unique(unlist(lapply(global,function(x) namePartition(makePartition(x)))))

  #Solution
  oh <- tapply(cbind(solution),solution,function(x) optimalHyp(part = x,jtrain = jtrain))
  oh <- lapply(oh,function(x) data.frame("X" = X,"fx" = plyr::mapvalues(factor(x[,2]),c("1","0"),Y)))
  error <- unique(tapply(cbind(solution),solution,function(x) getError(part = getPartition(x),jtrain = jtrain,jval = jval)))
  names(error) <- NULL

  #Timing
  end_time <- Sys.time()

  #Classifier
  classifier <- function(x){
    y <- unlist(lapply(unique(oh),function(y) y[y[,1] == x,2]))
    y <- table(y)
    y <- names(y)[y == max(y)]
    if(length(y) > 1)
      y <- sample(y,1)
    return(y)
  }

  return(list("classifier" = classifier,"hypotheses" = unique(oh),"partitions" = solution,"error" = error,
              "exhausted" = vis,"remain" = ifelse(optimal,remain,NA),"finished" = finished,"SLMvis" = ifelse(!finished,SLMvis,NA),
              "remain_after_prune" = ifelse(optimal,remain_after_prune,NA),"exhausted_until_prune" = exhausted_until_prune,
              "time" = end_time - start_time,"optimal" = optimal))
}
