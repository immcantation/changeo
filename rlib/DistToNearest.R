#' Generates distance to nearest neighbor by using S5F mutability model
#' 
#' @author     Namita Gupta, Gur Yaari, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2013.12.09
setwd('/home/ng265/workspace/changeo_v4/rlib')
library(stringr)
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


#' Get distance to nearest sequence sharing same V/J/JuncLen
#'
#' @param   file        DB file to read in
#' @param   db          dataframe if file has already been read in
#' @param   genotyped   whether file is genotyped
#' @param   grouping    'first' means first gene call is used
#' @param   dist        'S5F' uses S5F distance function, anything else uses old function
#' @return  DB dataframe with DIST_NEAREST column added
distToNearest <- function(file="", db=NA, genotyped=F, grouping='first') {
	
	if(file != "") {
    # Read input file into dataframe
	  db <- read.delim(file, as.is=T)
	} else if(!is.data.frame(db)) {
    stop('Must submit either DB file or data frame')
	}
	# Create new column for distance to nearest neighbor
	db$DIST_NEIGHBOR = rep(NA, nrow(db))
	
	if(genotyped) { 
		v_col <- "V_CALL_GENOTYPED"
	} else {
		v_col <- "V_CALL"
	}
	j_col <- "J_CALL"

	# Recreate V column based on grouping method
	if(grouping == 'first') {
		V <- as.character(sapply(db[,v_col], function(x){
					return(unique(na.omit(str_extract(strsplit(x,',')[[1]][1], 
													perl('IG[HLK][VDJ]\\d+[-/\\w]*'))))) }))
		J <- as.character(sapply(db[,j_col], function(x){
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
	for(k in 1:nrow(db)) {
		if(k %% 500 == 0) { cat(100*k/nrow(db),"%\n") }
		if(is.na(db[k,dist])) {
			indices = which(V==V[k] & J==J[k] & db[,junc_len]==db[k,junc_len])
			if(length(indices)>1) {
				Strings <- db[indices,"JUNCTION"]
				Strings <- toupper(Strings)
				N<-length(Strings)
				Clone1 <- sapply(Strings,s2c)
				nameStrings = paste(rep("A",N),1:N,sep="")
				Mat<-sapply(1:N, function(i)c(rep.int(0,i-1),sapply(i:N,function(j){
												old_dist_seq_fast(Clone1[,i],Clone1[,j])
											})))
				Mat <- Mat + t(Mat)
				colnames(Mat)<-indices
				rownames(Mat)<-indices
				db[indices,dist] <- sapply(1:N, function(i){ min(Mat[-i,i]) })
			}
		}
	}
	return(db)
}

