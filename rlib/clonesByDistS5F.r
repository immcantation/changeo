#' Generates clones by distance method with S5F mutability model
#' 
#' @author     Gur Yaari, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2013.12.09

load("S5F_Substitution.RData")

dist_seq_fast<-function(seq1,seq2){  
  retVal = 0
  sapply(1:length(seq1),function(x){
    i = seq1[x]
    j = seq2[x]
    if(nchar(i)==1) {
      retVal <<- retVal + mean(S5F_Substitution[i,j],S5F_Substitution[j,i],na.rm=T)
    } else {
      retVal <<- retVal + mean( S5F_Substitution[substr(j,3,3),i], S5F_Substitution[substr(i,3,3),j], na.rm=T )
    }    
  })
  #sum(symmetric_distance_array[paste(seq1,seq2)])  
  return(retVal)
}

switch_row <- function(Mat, i, j){
	if (i != j) {
		origNames <- rownames(Mat)[c(i,j)]
		tmp<-Mat[j,]
		Mat[j,] <- Mat[i,]
		Mat[i,] <- tmp
		rownames(Mat)[c(j,i)] <- origNames 
	}
	return(Mat)
}

switch_col <- function(Mat, i, j){
	if (i != j) {
		origNames <- colnames(Mat)[c(i,j)]
		tmp<-Mat[,j]
		Mat[,j] <- Mat[,i]
		Mat[,i] <- tmp
		colnames(Mat)[c(j,i)] <- origNames
	}
	return(Mat)
}

ComparStrings<-function(String1,String2,Threshold=2){
  tmp<-0
  for(i in 1:nchar(String1)){
    if(substr(String1,i,i)!=substr(String2,i,i)){
      tmp<-tmp+1;
      if(tmp>=Threshold)break
    }
  }
  if(tmp<Threshold)return(1)
  else return(0)
}

s2c<-function(String){
  N<-nchar(String)
  return(substring(String,1:N,1:N))
}

#' Generates clones by distance method with S5F mutability model
#' 
#' @param   Strings   a vector of junction sequences as strings
#' @param   Thresh    a numerical distance threshold
#' @return  a list of junction string vectors defining clones
getClones <- function(Strings, Thresh) {

  #   arg <- commandArgs(TRUE) 
  #   arg<-c("5","AAAAAAAAA|TTTTTTAAA|ccccccccA|ttttttttA|AAcAgAAAA|ttttttttt|ggggggggg|gggAagggg|tAAAAAAAA|AAtAAAAAA")
  #   arg<-c("10","TGTGCGAAAGAGGGATATTGTAGTAGTACCAGCTGTTTATATCGGGAACCTTTTGATATCTG|TGTGCGAGGGATCCGGGGCTATATTGTAGTGGTGGTGGCTGCGCGAATGCTTTTGATGTTTG|")
  #   Thresh <- as.numeric(arg[1])
  #   Strings <- toupper(strsplit(arg[2],"\\|")[[1]])
  Strings <- toupper(Strings)
  
  N<-length(Strings)
  Mat<-diag(N)
  
  Clone1 <- sapply(Strings,function(x){  
                            lenString <- nchar(x)
                            
                            #Pos 1
                            pos1 =  substr(x,1,1)
                            
                            #Pos 2
                            pos2 =  substr(x,2,2)  
                            
                            #Middle
                            posMid <- 3:(lenString-2)
                            posMiddle <-  substr(rep(x,lenString-5),(posMid-2),(posMid+2))  
                            
                            #Pos N-1
                            posN_1 <- substr(x,lenString-2,lenString-2)  
                            
                            #Pos N
                            posN <- substr(x,lenString-1,lenString-1)  
                            
                            return( c(pos1, pos2, posMiddle, posN_1, posN) )
                            })
  
  nameStrings = paste(rep("A",N),1:N,sep="")
  t<-sapply(1:N, function(i)c(rep.int(0,i-1),sapply(i:N,function(j){
                                                          dist_seq_fast(Clone1[,i],Clone1[,j])
                                                    })))
  colnames(t)<-Strings
  rownames(t)<-Strings
  BinaryDist<-t
  BinaryDist[BinaryDist<Thresh & BinaryDist>0]<-1
  BinaryDist[BinaryDist>=Thresh]<-0
  diag(BinaryDist) <- rep(1,N)
  tmp <- sapply(1:nrow(BinaryDist),function(i)BinaryDist[1:i,i]<<-BinaryDist[i,1:i])
  Mat <- BinaryDist
  
  #colnames(Mat) = Strings
  #rownames(Mat) = Strings
  if(N>2) {
    Blocks<-NULL
    Current<-2
    for(i in 1:(N-1)) {  
    while(Mat[min(Current,N),i]==1 & Current<N)Current=Current+1;
      indeces<-which(Mat[min(Current,N):N,i]==1)
      #print(indeces)
      for(Ind in indeces) {
        if(Current>N) break    
        Mat <- switch_row(Mat, Current+Ind-match(Ind,indeces), Current)
        Mat <- switch_col(Mat, Current+Ind-match(Ind,indeces), Current)
        Current<-Current+1
      }
      if(Current==i+1) {
        Blocks<-c(Blocks,i)
        Current<-i+2
      }
      if(Current>=N)break
    }
      if(sum(Mat[,N]==0)==N-1)  Blocks <- c(Blocks,N-1)
    Blocks <- c(Blocks,N)
  
  } else if(N==2) {
    if(Mat[2,1]==0)Blocks=c(1,2)
    else Blocks=2 
  } else {
    Blocks=1
  }
             
  #retString <- NULL
  
  #startInd=c(1,1+Blocks[-length(Blocks)])
  #for(i in 1:length(startInd)){
  #  retString <- c(retString, "|", rownames(Mat)[(startInd[i]:Blocks[i])])
  #}
  #retString <- paste(retString,collapse=",")
  #retString <- gsub("|,","|",retString,fixed=T)
  #retString <- gsub(",|","|",retString,fixed=T)
  #cat("$$$",retString<-substring(retString,2,nchar(retString)),"\n",sep="")

	returnBlocks <- list()
	startInd=c(1,1+Blocks[-length(Blocks)])
	for(i in 1:length(startInd)) {
    		returnBlocks[[i]] <- rownames(Mat)[(startInd[i]:Blocks[i])]
	}
	return(returnBlocks)

}
