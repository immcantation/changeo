#' Generates distance to nearest neighbor by using S5F mutability model
#' 
#' @author     Namita Gupta, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2013.12.09
library(stringr)
load("S5F_Substitution.RData")

dist_seq_fast<-function(seq1,seq2){  
  retVal = 0
  sapply(1:length(seq1),function(x){
    i = seq1[x]
    j = seq2[x]
    if(nchar(i)==1) {
      retVal <- retVal + mean(S5F_Substitution[i,j],S5F_Substitution[j,i],na.rm=T)
    } else {
      print(i)
      print(j)
      retVal <- retVal + mean( S5F_Substitution[substr(j,3,3),i], S5F_Substitution[substr(i,3,3),j], na.rm=T )
    }    
  })
  #sum(symmetric_distance_array[paste(seq1,seq2)])  
  return(retVal)
}

distToNearest <- function(file, genotyped=F, grouping='first') {
	
	# Read input file into dataframe
	clip <- read.delim(file, as.is=T)
	# Create new column for distance to nearest neighbor
	clip$DIST_NEIGHBOR = rep(NA, nrow(clip))
	
	if(genotyped) { 
		v_col <- "V_CALL_GENOTYPED"
	} else {
		v_col <- "V_CALL"
	}
	j_col <- "J_CALL"

	# Recreate V column based on grouping method
	if(grouping == 'first') {
		V <- sapply(clip[,v_col], function(x){
					return(unique(na.omit(str_extract(strsplit(x,',')[[1]][1], 
													perl('IG[HLK][VDJ]\\d+[-/\\w]*'))))) })
		J <- sapply(clip[,j_col], function(x){
					return(unique(na.omit(str_extract(strsplit(x,',')[[1]][1], 
													perl('IG[HLK][VDJ]\\d+[-/\\w]*'))))) })
	} else if(grouping == 'set') {
		stop("Error: Set grouping is not yet implemented\n")
	} else {
		stop("Error: Unrecognized grouping method\n")
	}
	
	junc_len <- "JUNCTION_GAP_LENGTH"
	dist <- "DIST_NEIGHBOR"
	
	# Iterate over rows
	for(i in 1:nrow(clip)) {
		if(is.na(clip[i,dist])) {
			indices = which(V==V[i] & J==J[i] & clip[,junc_len]==clip[i,junc_len])
      if(length(indices)>1) {
  			Strings <- clip[indices,"JUNCTION"]
  			Strings <- toupper(Strings)
  			N<-length(Strings)
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
  			Mat<-sapply(1:N, function(i)c(rep.int(0,i-1),sapply(i:N,function(j){
  											dist_seq_fast(Clone1[,i],Clone1[,j])
  										})))
  			Mat <- Mat + t(Mat)
  			colnames(Mat)<-indices
  			rownames(Mat)<-indices
  			clip[indices,dist] <- sapply(1:N, function(i){ min(Mat[-i,i]) })
      }
		}
	}
	return(clip)
}