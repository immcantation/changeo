#' Generates distance to nearest neighbor by using S5F mutability model
#' 
#' @author     Namita Gupta, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2013.12.09
setwd('/home/ng265/workspace/changeo_v4/rlib')
library(stringr)
load("S5F_Substitution.RData")
load("Old_Symmetric_Distance.RData")

#' String to Character
#'
#' @param   String    a string
#' @return  a vector of each character
s2c<-function(String){
  N<-nchar(String)
  return(substring(String,1:N,1:N))
}

#' Get symmetric distance
#'
#' @param   seq1   the first nucleotide sequence
#' @param   seq2   the second nucleotide sequence
#' @return  distance between two sequences based on symmetric distance matrix
old_dist_seq_fast<-function(seq1,seq2){
  sum(symmetric_distance_array[paste(seq1,seq2)])
}

#' Get S5F distance
#'
#' @param   seq1   the first nucleotide sequence
#' @param   seq2   the second nucleotide sequence
#' @return  distance between two sequences based on S5F model
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

#' Get distance to nearest sequence sharing same V/J/JuncLen
#'
#' @param   file        CLIP-DB file to read in
#' @param   clip        CLIP dataframe if it has already been read in
#' @param   genotyped   whether file is genotyped
#' @param   grouping    'first' means first gene call is used
#' @param   dist        'S5F' uses S5F distance function, anything else uses old function
#' @return  CLIP dataframe with DIST_NEAREST column added
distToNearest <- function(file="", clip=NA, genotyped=F, grouping='first', dist="S5F") {
	
	if(file != "") {
    # Read input file into dataframe
	  clip <- read.delim(file, as.is=T)
	} else if(!is.data.frame(clip)) {
    stop('Must submit either Clip file or data frame')
	}
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
		V <- as.character(sapply(clip[,v_col], function(x){
					return(unique(na.omit(str_extract(strsplit(x,',')[[1]][1], 
													perl('IG[HLK][VDJ]\\d+[-/\\w]*'))))) }))
		J <- as.character(sapply(clip[,j_col], function(x){
					return(unique(na.omit(str_extract(strsplit(x,',')[[1]][1], 
													perl('IG[HLK][VDJ]\\d+[-/\\w]*'))))) }))
	} else if(grouping == 'set') {
		stop("Error: Set grouping is not yet implemented\n")
	} else {
		stop("Error: Unrecognized grouping method\n")
	}
	
	junc_len <- "JUNCTION_GAP_LENGTH"
	dist <- "DIST_NEIGHBOR"
	
	# Iterate over rows
	for(k in 1:nrow(clip)) {
    if(k %% 500 == 0) { cat(100*k/nrow(clip),"%\n") }
		if(is.na(clip[k,dist])) {
			indices = which(V==V[k] & J==J[k] & clip[,junc_len]==clip[k,junc_len])
      if(length(indices)>1) {
  			Strings <- clip[indices,"JUNCTION"]
  			Strings <- toupper(Strings)
  			N<-length(Strings)
        if(dist=="S5F") {
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
        } else {
          Clone1 <- sapply(Strings,s2c)
          nameStrings = paste(rep("A",N),1:N,sep="")
          Mat<-sapply(1:N, function(i)c(rep.int(0,i-1),sapply(i:N,function(j){
                          old_dist_seq_fast(Clone1[,i],Clone1[,j])
                        })))
          
        }
  			Mat <- Mat + t(Mat)
  			colnames(Mat)<-indices
  			rownames(Mat)<-indices
  			clip[indices,dist] <- sapply(1:N, function(i){ min(Mat[-i,i]) })
      }
		}
	}
	return(clip)
}

