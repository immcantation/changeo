#' Generates clones by distance method with mutability model
#' 
#' @author     Gur Yaari, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2013.11.10

#dyn.load("/home/gur/SHMDistance.so")
#arg<-c("5","AAAAAAAAA|TTTTTTAAA|ccccccccA|ttttttttA|AAcAgAAAA|ttttttttt|ggggggggg|gggAagggg|tAAAAAAAA|AAtAAAAAA")
#arg <- commandArgs(TRUE) 
#J_Length=50
# N<-5000
# arg<-c("5",paste(sapply(1:N,function(i)paste(sample(NUCLEOTIDES,J_Length,replace=TRUE),collapse="")),collapse="|"))
# 
# arg<-c("5","tgtgcgagaactggtacggtggtaacgtcagggtactactacggaatggacgtctgg|tgtgcggtgactacggtggagactccgatgttccagtcctacggtatgaacgtctgg")
#source("http://selection.med.yale.edu/baseline/Baseline_Functions.r")

dyn.load("SHMDistance.so")

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


#' Generates clones by old distance method
#' 
#' @param   Strings   a vector of junction sequences as strings
#' @param   Thresh    a numerical distance threshold
#' 
#' @return  a list of junction string vectors defining clones
getClones <- function(Strings, Thresh) {
	Strings <- toupper(Strings)
	StringsNumeric<-chartr(c("ACGTN.-"),"1234566",Strings)
	
	N<-length(Strings)
	Mat<-diag(N)
	
	tmp<-.C("SHMDistance", Nstrings=as.integer(length(StringsNumeric[1:N])), Nnucs=as.integer(nchar(StringsNumeric[1])),S=StringsNumeric[1:N],MATRIX=matrix(0,N,N))
	
	t<-tmp[[4]]
	
	# Create binary matrix based on whether values are below threshold
	colnames(t)<-Strings
	rownames(t)<-Strings
	BinaryDist<-t
	BinaryDist[BinaryDist>=Thresh]<-0
	BinaryDist[BinaryDist<Thresh & BinaryDist>0]<-1
	diag(BinaryDist) <- rep(1,N)
	tmp <- sapply(1:nrow(BinaryDist),function(i)BinaryDist[1:i,i]<-BinaryDist[i,1:i])
	Mat <- BinaryDist
	
	# Rearrange rows/columns (essentially single linkage clustering) to form blocks of junctions
	if(N>2) {
		Blocks<-NULL
		Current<-2
		for(i in 1:(N-1)){  
			while(Mat[min(Current,N),i]==1 & Current<N) Current=Current+1;
			indices<-which(Mat[min(Current,N):N,i]==1)
			#print(indices)
			for(Ind in indices){
				if(Current>N) break    
				Mat <- switch_row(Mat, Current+Ind-match(Ind,indices), Current)
				Mat <- switch_col(Mat, Current+Ind-match(Ind,indices), Current)
				Current<-Current+1
			}
			if(Current==i+1){
				Blocks<-c(Blocks,i)
				Current<-i+2
			}
			if(Current>=N)break
		}
		if(sum(Mat[,N]==0)==N-1) Blocks <- c(Blocks,N-1)
		Blocks <- c(Blocks,N)	
	} else if (N==2) {
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
