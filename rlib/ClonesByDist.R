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

	Strings <- toupper(Strings) #Convert junctions to upper-case
  StringsOrig <- Strings #Save original Strings

  # Change '.' gaps to 'N' - IGNORES GAPS FOR NOW TO ELIMINATE NA SITUATION
	Strings <- gsub('.', 'N', Strings, fixed=T)

  # Add "NN" to the start and end of each sequence (junction)
  # This helps determine targeting for the first 2 and the last 2 nucletides
	Strings <- as.vector(sapply(Strings,function(x){paste("NN",x,"NN",sep="")}))

	N<-length(Strings)
	Mat<-diag(N)

  # Break the junctions into 5-mers and create a sliding window matrix
  # (each column is a sequence)
	matSeqSlidingFiveMer <- sapply(Strings,function(x){
	  slidingArrayOf5mers(x)
			},simplify="matrix")

  # Compute pairwise distance between all sequences' fivemers (by column)
  BinaryDist <- sapply(1:N, function(i)c(rep.int(0,i-1),sapply(i:N,function(j){
									dist_seq_fast(matSeqSlidingFiveMer[,i], matSeqSlidingFiveMer[,j],
									              model_data[["subs"]],
                                model_data[["mut"]])
								})))
	colnames(BinaryDist)<-StringsOrig
	rownames(BinaryDist)<-StringsOrig

  #colnames(BinaryDist)<- 1:ncol(BinaryDist)
  #rownames(BinaryDist)<- 1:nrow(BinaryDist)

  # Perform single linkage clustering to determine clones
  hc <- hclust(as.dist(BinaryDist), method="single")
  #plot(hc)
  clones <- cutree(hc,h=Thresh)
  listClones <- list()
  for(i in unique(sort(clones))){
    listClones[[i]] <- names(clones)[clones==i]
  }

	return(listClones)
}
