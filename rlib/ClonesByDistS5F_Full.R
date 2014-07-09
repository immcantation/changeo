#' Generates clones by distance method with S5F mutability model
#' 
#' @author     Gur Yaari, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2014.7.8

load("S5F_Targeting.RData")
S5F_Substitution <- S5F_Targeting[["Substitution"]]
S5F_Mutability <- S5F_Targeting[["Mutability"]]  
S5F_Substitution_Array <- S5F_Targeting[["S5F_Substitution_Array"]]


#' Get S5F distance between two sequences of same length broken down into 5-mers
#'
#' @param   seq1   the first nucleotide sequence
#' @param   seq2   the second nucleotide sequence
#' 
#' @return  distance between two sequences based on S5F model
dist_seq_fast<-function(seq1,seq2) {
	# Compute distance only on fivemers that have mutations
	fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
	seq1 <- seq1[fivemersWithMu]
	seq2 <- seq2[fivemersWithMu]
	fivemersWithMu <- (substr(seq1,3,3)!="N" & substr(seq2,3,3)!="N")
	seq1 <- seq1[fivemersWithMu]
	seq2 <- seq2[fivemersWithMu]  
	a <- tryCatch({
				if(length(seq1)==1) {
					seq1_to_seq2 <- S5F_Substitution[substr(seq2,3,3),seq1] * S5F_Mutability[seq1]
					seq2_to_seq1 <- S5F_Substitution[substr(seq1,3,3),seq2] * S5F_Mutability[seq2]
				} else {
					seq1_to_seq2 <- sum( diag(S5F_Substitution[substr(seq2,3,3),seq1]) *  S5F_Mutability[seq1] )
					seq2_to_seq1 <- sum( diag(S5F_Substitution[substr(seq1,3,3),seq2]) *  S5F_Mutability[seq2] )
				}
				return( mean(c(seq1_to_seq2, seq2_to_seq1)) )
			},error = function(e){
				return(NA)
			})
}


#' Switch two named rows in a matrix
#'
#' @param   Mat  the matrix whose rows are to be switched
#' @param   i    the first row to switch
#' @param   j    the second row to switch
#' 
#' @return  matrix with rows i and j switched and corrected rownames
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


#' Switch two named columns in a matrix
#'
#' @param   Mat  the matrix whose columns are to be switched
#' @param   i    the first column to switch
#' @param   j    the second column to switch
#' 
#' @return  matrix with columns i and j switched and corrected columnnames
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


#' Generates clones by distance method with S5F mutability model
#' 
#' @param   Strings   a vector of junction sequences as strings
#' @param   Thresh    a numerical distance threshold
#' 
#' @return  a list of junction string vectors defining clones
getClones <- function(Strings, Thresh) {

	#   arg <- commandArgs(TRUE) 
	#   arg<-c("5","AAAAAAAAA|TTTTTTAAA|ccccccccA|ttttttttA|AAcAgAAAA|ttttttttt|ggggggggg|gggAagggg|tAAAAAAAA|AAtAAAAAA")
	#   arg<-c("10","ACGNACGT|ACGTACGG|ACGTACGG|ACGTACGT|ATGTACGT|ACGTACGT|ACGTACGT|ACGTACGT")
	#   Thresh <- as.numeric(arg[1])
	#   Strings <- toupper(strsplit(arg[2],"\\|")[[1]])
	
	Strings <- toupper(Strings)
	  
	# Add "NN" to the start and end of each sequence (junction)
	StringsOrig <- Strings
	Strings <- as.vector(sapply(Strings,function(x){paste("NN",x,"NN",sep="")}))
	    
	N<-length(Strings)
	Mat<-diag(N)

	Clone1 <- sapply(Strings,function(x){  
                              lenString <- nchar(x)
                              fivemersPos <- 3:(lenString-2)
                              fivemers <-  substr(rep(x,lenString-4),(fivemersPos-2),(fivemersPos+2))
                              return(fivemers)
                             },simplify="matrix")
  
	t<-sapply(1:N, function(i) c(rep.int(0,i-1),
					             sapply(i:N,function(j){
                                              dist_seq_fast(Clone1[,i],Clone1[,j])
                                            })))
	  
	# Create binary matrix based on whether values are below threshold
	BinaryDist <- t
	colnames(BinaryDist) <- StringsOrig
	rownames(BinaryDist) <- StringsOrig
	BinaryDist[BinaryDist<Thresh & BinaryDist>0] <- 1
	BinaryDist[BinaryDist>=Thresh] <- 0
	diag(BinaryDist) <- rep(1,N)
	tmp <- sapply(1:nrow(BinaryDist),function(i)BinaryDist[1:i,i] <<- BinaryDist[i,1:i])
	Mat <- BinaryDist
	  
	# Rearrange rows/columns (essentially single linkage clustering) to form blocks of junctions
	if(N>2) {
	 	Blocks<-NULL
		Current<-2
		for(i in 1:(N-1)) {  
			while(Mat[min(Current,N),i]==1 & Current<N) Current=Current+1;
			indices <- which(Mat[min(Current,N):N,i]==1)
			#print(indices)
			for(Ind in indices) {
				if(Current>N) break    
				Mat <- switch_row(Mat, Current+Ind-match(Ind,indices), Current)
				Mat <- switch_col(Mat, Current+Ind-match(Ind,indices), Current)
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
		if(Mat[2,1]==0) Blocks=c(1,2)
		else Blocks=2 
	} else {
		Blocks=1
	}
	  
	# Create list of junction string vectors that are single clones
	returnBlocks <- list()
	startInd=c(1,1+Blocks[-length(Blocks)])
	for(i in 1:length(startInd)) {
		returnBlocks[[i]] <- rownames(Mat)[(startInd[i]:Blocks[i])]
	}
	return(returnBlocks)
}
