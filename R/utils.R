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

  #Name partition
  names(part) <- unlist(lapply(part,function(x) paste(x,collapse = ",")))
  return(part)
}

#Teste if two partitions are related
is.related <- function(part1,part2){
  if(all.equal(part1,part2) == T)
    return(T)
  m1 <- max(part1)
  m2 <- max(part2)
  if(m1 == m2)
    return(F)
  if(m1 > m2)
    return(sum(!unlist(lapply(tapply(part2,part1,function(x) x),function(x) length(unique(x)) == 1))) == 0)
  if(m1 < m2)
    return(sum(!unlist(lapply(tapply(part1,part2,function(x) x),function(x) length(unique(x)) == 1))) == 0)
}

numberNeighboors <- function(part){
  #Number of blocks unions
  numb_union <- (length(part)*(length(part) - 1))/2

  #Number of blocks breaks
  part <- lapply(part,length)
  part <- part[part > 1]
  if(length(part) > 0)
    numb_break <- sum(unlist(lapply(part,function(x) ifelse(x %% 2 == 0,sum(choose(x,0:(x/2 - 1))) + choose(x,x/2)/2,
                                  sum(choose(x,0:floor(x/2)))))))
  else
    numb_break <- 0
  return(list("total" = numb_union+numb_break,"union" = numb_union,"breaks" = numb_break))
}

#Find neighbors of a partition
findNeighbors <- function(part,n,cores){

  #All breaks of partition
  m <- matrix(ncol = n)
  for(i in 1:length(part))
    m[1,part[[i]]] <- i
  size <- max(m)
  nsize <- size + 1

  if(max(m) < n){
    for(i in 1:max(m)){
      if(length(m[1,m[1,] == i]) > 1){
        mtemp <- expand.grid(rlist::list.append(c(1),replicate(length(m[1,m[1,] == i])-1, c(0,1), simplify=FALSE)))
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
    breaks <- unlist(mclapply(data.frame(t(m)),nameFromIndex,mc.cores = cores))[-1]
  }
  else
    breaks <- NULL

  #All unions of partition
  if(length(part) > 1){
    nameUnion <- combn(1:length(part),2)
    if(ncol(nameUnion) > 1){
      union <- plyr::alply(nameUnion,2,function(x) c(x[[1]],x[[2]],part[[x[1]]],part[[x[2]]]))
      union <- mclapply(union,function(x) rlist::list.append(list(x[-c(1,2)]),part[-c(x[1:2])]),mc.cores = cores)
      union <- mclapply(union,function(y) lapply(rapply(y, enquote, how="unlist"), eval),mc.cores = cores)
      union <- mclapply(union,makePartition,mc.cores = cores)
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

#Sample neighbors of a partition
sampleNeighbors <- function(sample_size,part,n,cores){
  numb <- numberNeighboors(part)
  if(numb$total < sample_size)
    return(findNeighbors(part,n,cores))

  #Sample sizes
  if(numb$breaks < numb$union){
    sbreak <- floor(sample_size * numb$breaks/numb$total)
    sbreak <- ifelse(sbreak <= numb$breaks,sbreak,numb$breaks)
    sunion <- sample_size - sbreak
  }
  else{
    sunion <- floor(sample_size * numb$union/numb$total)
    sunion <- ifelse(sunion <= numb$union,sunion,numb$union)
    sbreak <- sample_size - sunion
  }

  #All breaks of partition
  m <- matrix(ncol = n)
  for(i in 1:length(part))
    m[1,part[[i]]] <- i
  size <- max(m)
  nsize <- size + 1
  avai_blocks <- names(table(m[1,]))[table(m[1,]) > 1]
  wheight <- table(m[1,])[table(m[1,]) > 1]

  mtemp <- mclapply(data.frame(rbind(c(1:sbreak))),function(x){
    #Sample block
    block <- sample(x = avai_blocks,size = 1,prob = wheight)

    #Size new
    size_new <- sample(x = 1:(sum(m[1,] == block)-1),size = 1)

    #Which to turn
    turn <- sample(x = c(1:ncol(m))[m[1,] == block],size = size_new)

    #New block
    mtemp <- m[1,]
    mtemp[turn] <- nsize
    return(as.list(mtemp))},mc.cores = cores)
  m <- as.matrix(data.table::rbindlist(mtemp))

  if(nrow(m) > 1)
    breaks <- unlist(mclapply(data.frame(t(m)),nameFromIndex,mc.cores = cores))
  else
    breaks <- NULL

  #All unions of partition
  if(length(part) > 1 & sunion > 1){
    nameUnion <- combn(1:length(part),2)
    nameUnion <- nameUnion[,sample(x = 1:ncol(nameUnion),size = sunion)]
    if(ncol(nameUnion) > 1){
      union <- plyr::alply(nameUnion,2,function(x) c(x[[1]],x[[2]],part[[x[1]]],part[[x[2]]]))
      union <- mclapply(union,function(x) rlist::list.append(list(x[-c(1,2)]),part[-c(x[1:2])]),mc.cores = cores)
      union <- mclapply(union,function(y) lapply(rapply(y, enquote, how="unlist"), eval),mc.cores = cores)
      union <- mclapply(union,makePartition,mc.cores = cores)
      union <- unlist(lapply(union,namePartition))
      names(union) <- NULL
    }
    else{
      union <- namePartition(makePartition(list(c(part[[1]],part[[2]]))))
    }

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

#Sample a partition
samplePartition <- function(n){
  #Size of new partion
  size <- rbinom(1,n,1/2)

  if(size == 1)
    return(makePartition(list(1:n)))
  if(size == n)
    return(makePartition(as.list(1:n)))

  #Sample an ordenation of the domain
  o <- sample(1:n,n,F)

  #Sample last position of each partition
  pos <- sample(1:(n-1),size-1,F)
  pos <- pos[order(pos,decreasing = F)]

  #Creating part
  part <- list()
  pos[size] <- n
  for(i in 1:size)
    part[[i]] <- o[ifelse(i == 1,1,pos[i-1]+1):pos[i]]
  part <- makePartition(part)

  return(part)
}

