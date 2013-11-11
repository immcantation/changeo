#' Generates clones by distance method with mutability model
#' 
#' @author Gur Yaari, Namita Gupta, Jason Vander Heiden
#' @date 2013.11.10
#' 
#' @param Strings   a named vector of junction sequences
#' @param Thresh    mutation threshold
#' 
#' @return a string of grouped junctions

#dyn.load("/home/gur/SHMDistance.so")
#arg<-c("5","AAAAAAAAA|TTTTTTAAA|ccccccccA|ttttttttA|AAcAgAAAA|ttttttttt|ggggggggg|gggAagggg|tAAAAAAAA|AAtAAAAAA")
#arg <- commandArgs(TRUE) 
#J_Length=50
# N<-5000
# arg<-c("5",paste(sapply(1:N,function(i)paste(sample(NUCLEOTIDES,J_Length,replace=TRUE),collapse="")),collapse="|"))
# 
# arg<-c("5","tgtgcgagaactggtacggtggtaacgtcagggtactactacggaatggacgtctgg|tgtgcggtgactacggtggagactccgatgttccagtcctacggtatgaacgtctgg")
#source("http://selection.med.yale.edu/baseline/Baseline_Functions.r")
source("Baseline_Functions.R")
dyn.load("SHMDistance.so")


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



getClones <- function(Strings, Thresh) {
	#Thresh <- as.numeric(arg[1])
	#Strings <- toupper(strsplit(arg[2],"\\|")[[1]])
	Strings <- toupper(strsplit(Strings,"\\|")[[1]])
	StringsNumeric<-chartr(c("ACGTN-"),"123456",Strings)
	
	N<-length(Strings)
	Mat<-diag(N)
	
	tmp<-.C("SHMDistance", Nstrings=as.integer(length(StringsNumeric[1:N])), Nnucs=as.integer(nchar(StringsNumeric[1])),S=StringsNumeric[1:N],MATRIX=matrix(0,N,N))
	
	t<-tmp[[4]]
	
	colnames(t)<-Strings
	rownames(t)<-Strings
	BinaryDist<-t
	BinaryDist[BinaryDist<Thresh & BinaryDist>0]<-1
	BinaryDist[BinaryDist>=Thresh]<-0
	diag(BinaryDist) <- rep(1,N)
	tmp <- sapply(1:nrow(BinaryDist),function(i)BinaryDist[1:i,i]<<-BinaryDist[i,1:i])
	Mat <- BinaryDist
	
	if(N>2) {
		Blocks<-NULL
		Current<-2
		for(i in 1:(N-1)){  
			while(Mat[min(Current,N),i]==1 & Current<N)Current=Current+1;
			indeces<-which(Mat[min(Current,N):N,i]==1)
			#print(indeces)
			for(Ind in indeces){
				if(Current>N) break    
				Mat <- switch_row(Mat, Current+Ind-match(Ind,indeces), Current)
				Mat <- switch_col(Mat, Current+Ind-match(Ind,indeces), Current)
				Current<-Current+1
			}
			if(Current==i+1){
				Blocks<-c(Blocks,i)
				Current<-i+2
			}
			if(Current>=N)break
		}
		if(sum(Mat[,N]==0)==N-1)  Blocks <- c(Blocks,N-1)
		Blocks <- c(Blocks,N)
		
	} else if (N==2) {
		if(Mat[2,1]==0) Blocks=c(1,2)
		else Blocks=2 
	} else {
		Blocks=1
	}
	
	retString <- NULL
	
	startInd=c(1,1+Blocks[-length(Blocks)])
	for(i in 1:length(startInd)){
		retString <- c(retString, "|", rownames(Mat)[(startInd[i]:Blocks[i])])
	}
	retString <- paste(retString,collapse=",")
	retString <- gsub("|,","|",retString,fixed=T)
	retString <- gsub(",|","|",retString,fixed=T)
	retString <- substring(retString,2,nchar(retString))
	#cat("$$$",retString<-substring(retString,2,nchar(retString)),"\n",sep="")

	return(retString)
}
