####Utility functions for partitionUcurve package#####

#Test if partion and order parts by minimum
makePartition <- function(part){
  #Is partition
  v <- unlist(part)
  if(!(min(v) == 1 & max(v) == length(v) & length(unique(v)) == length(v)))
    stop("It is not a partition. Please review the code.")

  #Order partition
  part <- part[order(unlist(lapply(part,min)))]
  part <- lapply(part,function(x) x[order(x,decreasing = F)])

  #Name partition
  names(part) <- unlist(lapply(part,function(x) paste(x,collapse = ",")))

  return(part)
}

#Get name of partition
namePartition <- function(part){
  return(paste(names(part),collapse = "-"))
}

#Get name from index
nameFromIndex <- function(index){
  n <- length(index)
  name <- namePartition(makePartition(tapply(c(1:n),index,function(x) x)))

  return(name)
}

#Get partition from name
getPartition <- function(name){
  part <- sapply(unlist(strsplit(x = name,split = "-")),function(x) lapply(strsplit(x,","),as.numeric))
  return(part)
}

#Teste if two partitions are related
is.related <- function(part1,part2){
  return(sum(part1 <= part2) == length(part1) | sum(part1 >= part2) == length(part1))
}

#Find neighbors of a partition
findNeighbors <- function(part){

  #All breaks of partition
  n <- max(unlist(part))
  m <- matrix(ncol = n)
  for(i in 1:length(part))
    m[1,part[[i]]] <- i
  size <- max(m)
  nsize <- size + 1

  if(max(m) < n){
    for(i in 1:max(m)){
      if(length(m[1,m[1,] == i]) > 1){
        mtemp <- expand.grid(replicate(length(m[1,m[1,] == i]), c(0,1), simplify=FALSE))
        mtemp <- mtemp[2:(nrow(mtemp)/2),]
        mtemp <- as.matrix(nsize*mtemp)
        mtemp[mtemp == 0] <- i
        m2 <- matrix(0,nrow = nrow(mtemp),ncol = n)
        m2[,m[1,] == i] <- as.matrix(mtemp)
        m3 <- t(replicate(nrow(m2),m[1,]))
        m3[,m[1,] == i] <- 0
        m2 <- m2 + m3
        m <- rbind(m,m2)
      }
    }
    breaks <- apply(m,1,nameFromIndex)[-1]
  }
  else
    breaks <- NULL

  #All unions of partition
  if(length(part) > 1){
    nameUnion <- combn(1:length(part),2)
    if(ncol(nameUnion) > 1){
      union <- plyr::alply(nameUnion,2,function(x) c(x[[1]],x[[2]],part[[x[1]]],part[[x[2]]]))
      union <- lapply(union,function(x) rlist::list.append(list(x[-c(1,2)]),part[-c(x[1:2])]))
      union <- lapply(union,function(y) lapply(rapply(y, enquote, how="unlist"), eval))
      union <- lapply(union,makePartition)
      union <- unlist(lapply(union,namePartition))
      names(union) <- NULL
    }
    else
      union <- namePartition(makePartition(list(c(part[[1]],part[[2]]))))
    if(is.null(breaks))
        return(union)
    else
      return(c(breaks,union))
  }
  else
    return(breaks)
}

#Joint distribution
jointDistribution <- function(x,y){
  return(prop.table(table(x,y)))
}

#Calculate error of a partition
getError <- function(part,jtrain,jval){

  #Estimated hypothesis
  parJoint <- lapply(part,function(x) colSums(rbind(jtrain[x,],c(0,0))))
  optimalComplement <- unlist(lapply(parJoint,function(x) c(1:2)[x == min(x)][1]))

  #Error
  parVal <- unlist(lapply(part,function(x) colSums(rbind(jval[x,],c(0,0)))))[optimalComplement + seq(0,2*(length(part)-1),2)]
  error <- sum(parVal)

  return(error)
}

#Add a visited node to file
addNode <- function(part,error){
  name <- paste(namePartition(part),"e",error,sep = "")
  write(name,file = "visited.dat",append = T)
}

#Test if node has been visited
visitedNode <- function(part){
  name <- namePartition(part)
  vis <- suppressWarnings(system(paste("grep ",name," visited.dat"),intern = T))
  if(length(vis) > 0)
    vis <- as.numeric(unlist(strsplit(vis,"e"))[2])
  else
    vis <- NULL

  return(vis)
}

#Get optimal hypothesis of partition
optimalHyp <- function(part,jtrain){
  part <- getPartition(part)

  #Estimated hypothesis
  parJoint <- lapply(part,function(x) colSums(rbind(jtrain[x,],c(0,0))))
  optimalComplement <- unlist(lapply(parJoint,function(x) c(0:1)[x == min(x)][1]))

  #Hypothesis
  domain <- unlist(strsplit(names(optimalComplement),","))
  domain <- domain[order(as.numeric(domain))]
  f <- vector()
  for(i in domain)
    f[domain == i] <- optimalComplement[unlist(lapply(strsplit(names(optimalComplement),","),function(x) i %in% x))]
  f <- data.frame("X" = domain,"fx" = f)

  return(f)
}
