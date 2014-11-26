#' Generates clones by distance method with S5F mutability model
#' 
#' @author     Gur Yaari, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2014.11.24

# Imports
suppressPackageStartupMessages(require(shm))
suppressPackageStartupMessages(require(alakazam))


#' Switch two rows in a matrix with named dimensions
#'
#' @param   Mat   input matrix
#' @param   i     row index
#' @param   j     row index
#' @return  matrix with rows i and j switched
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


#' Switch two columns in a matrix with named dimensions
#'
#' @param   Mat   input matrix
#' @param   i     column index
#' @param   j     column index
#' @return  matrix with columns i and j switched
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
#' @param   model     string defining name of the model to load.
#'                    One of "hs5f" or "m3n".
#' @return  a list of junction string vectors defining clones
getClones <- function(Strings, Thresh, model="hs5f") {

    model_data <- loadModel(model)

	Strings <- toupper(Strings)
	
	StringsOrig <- Strings
	# Change '.' gaps to '-'
	Strings <- gsub('.', '-', Strings, fixed=T)
	# Add "NN" to the start and end of each sequence (junction)
	Strings <- as.vector(sapply(Strings,function(x){paste("NN",x,"NN",sep="")}))
	
	N<-length(Strings)
	Mat<-diag(N)
	
	Clone1 <- sapply(Strings,function(x){  
				lenString <- nchar(x)
				fivemersPos <- 3:(lenString-2)
				fivemers <-  substr(rep(x,lenString-4),(fivemersPos-2),(fivemersPos+2))
				return(fivemers)
			},simplify="matrix")
	
	t<-sapply(1:N, function(i)c(rep.int(0,i-1),sapply(i:N,function(j){
									dist_seq_fast(Clone1[,i],Clone1[,j],
									              model_data[["sub"]], model_data[["mut"]])
								})))
	BinaryDist<-t
	colnames(BinaryDist)<-StringsOrig
	rownames(BinaryDist)<-StringsOrig
	BinaryDist[BinaryDist>=Thresh]<-0
	BinaryDist[BinaryDist<Thresh & BinaryDist>0]<-1
	diag(BinaryDist) <- rep(1,N)
	tmp <- sapply(1:nrow(BinaryDist),function(i)BinaryDist[1:i,i]<<-BinaryDist[i,1:i])
	Mat <- BinaryDist
	
	if(N>2) {
		Blocks<-NULL
		Current<-2
		for(i in 1:(N-1)) {  
			while(Mat[min(Current,N),i]==1 & Current<N)Current=Current+1;
			indices<-which(Mat[min(Current,N):N,i]==1)
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
		if(Mat[2,1]==0)Blocks=c(1,2)
		else Blocks=2 
	} else {
		Blocks=1
	}
	
	returnBlocks <- list()
	startInd=c(1,1+Blocks[-length(Blocks)])
	for(i in 1:length(startInd)) {
		returnBlocks[[i]] <- rownames(Mat)[(startInd[i]:Blocks[i])]
	}
	return(returnBlocks)
}