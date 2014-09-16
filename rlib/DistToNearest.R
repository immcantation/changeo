#' Generates distance to nearest neighbor by using S5F mutability model
#' 
#' @author     Namita Gupta, Gur Yaari, Mohamed Uduman
#' @copyright  Copyright 2013 Kleinstein Lab, Yale University. All rights reserved
#' @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
#' @date       2014.09.16

library(stringr)
library(utils)
library(parallel)
library(doMC) 
library(plyr)
registerDoMC(cores=detectCores())


#' Get S5F distance between two sequences of same length broken down into fivemers
#'
#' @param   seq1   the first nucleotide sequence
#' @param   seq2   the second nucleotide sequence
#' @return  distance between two sequences based on S5F model
dist_seq_fast<-function(seq1,seq2){  
  #Compute distance only on fivemers that have mutations
  fivemersWithMu <- substr(seq1,3,3)!=substr(seq2,3,3)
  fivemersWithNonNuc <- ( !is.na(match(substr(seq1,3,3),c("A","C","G","T"))) & !is.na(match(substr(seq2,3,3),c("A","C","G","T"))) )
  seq1 <- seq1[fivemersWithMu & fivemersWithNonNuc]
  seq2 <- seq2[fivemersWithMu & fivemersWithNonNuc] 
  a <- tryCatch({
  if(length(seq1)==1){
    #seq1_to_seq2 <- S5F_Substitution[substr(seq2,3,3),seq1] * S5F_Mutability[seq1]
    seq1_to_seq2 <- Targeting[["Targeting"]][substr(seq2,3,3),seq1]
    #seq2_to_seq1 <- S5F_Substitution[substr(seq1,3,3),seq2] * S5F_Mutability[seq2]
    seq2_to_seq1 <- Targeting[["Targeting"]][substr(seq1,3,3),seq2]
  }else{
    seq1_to_seq2 <- sum( diag(Targeting[["Targeting"]][substr(seq2,3,3),seq1]) )
    seq2_to_seq1 <- sum( diag(Targeting[["Targeting"]][substr(seq1,3,3),seq2]) )
  }
  return( mean(c(seq1_to_seq2, seq2_to_seq1)) )
  },error = function(e){
    return(NA)
  })
}


#' Given an array of junction sequenes, find the distance to the closest sequence
#'
#' @param   arrJunctions   
#' @return  distances to the closest sequence
getDistanceToClosest <- function(arrJunctions){ 
  
  #Initialize array of distances
  arrJunctionsDist <- rep(NA,length(arrJunctions))
  
  #Filter unique junctions
  arrJunctionsUnique <- unique(arrJunctions)
  
  #Map indexes of unique to its non-unique in the original arrJunctions
  indexJunctions <- match(arrJunctions, arrJunctionsUnique)
  
  #Identify junctions with multiple non-unique sequences and set its distances to 0 
  indexJunctionsCounts <- table(indexJunctions)
  indexRepeated <- as.numeric(names(indexJunctionsCounts)[indexJunctionsCounts>1])
  indexRepeated <- indexJunctions%in%indexRepeated
  arrJunctionsDist[ indexRepeated ] <- rep(0,sum(indexRepeated))
  names(arrJunctionsDist) <- arrJunctions
  
  #Compute distances between junctions
  numbOfUniqueJuctions <- length(arrJunctionsUnique)
  arrUniqueJunctionsDist <- rep(NA,numbOfUniqueJuctions)
  if(numbOfUniqueJuctions>1){
    arrJunctionsUnique <- toupper(arrJunctionsUnique)
    arrJunctionsUnique <- gsub('.', '-', arrJunctionsUnique, fixed=T)
    arrJunctionsUnique <- as.vector(sapply(arrJunctionsUnique,function(x){paste("NN",x,"NN",sep="")})) 
    matSequenceFivemers <- sapply(arrJunctionsUnique,
                                  function(x){  
                                    lenString <- nchar(x)
                                    fivemersPos <- 3:(lenString-2)
                                    fivemers <-  substr(rep(x,lenString-4),(fivemersPos-2),(fivemersPos+2))
                                    return(fivemers)
                                  }
                                  , simplify="matrix"
    )
    matDistance <-sapply(1:numbOfUniqueJuctions, function(i)c(rep.int(0,i-1),sapply(i:numbOfUniqueJuctions,function(j){
      dist_seq_fast(matSequenceFivemers[,i],matSequenceFivemers[,j])
    })))
    matDistance <- matDistance + t(matDistance)
    arrUniqueJunctionsDist <- sapply(1:numbOfUniqueJuctions, function(i){ min(matDistance[-i,i]) })    
    names(arrUniqueJunctionsDist) <- arrJunctionsUnique
  }
  
  #Fill the distances for the sequences that are unique
  arrJunctionsDist[is.na(arrJunctionsDist)] <- arrUniqueJunctionsDist[indexJunctionsCounts==1]
  return(round(arrJunctionsDist,4))
}


#' Get distance to nearest sequence sharing same V/J/JuncLen
#'
#' @param   db          dataframe
#' @param   genotyped   whether file is genotyped
#' @param   grouping    'first' means first gene call is used
#' @param   model       SHM targeting model 
#' @return  DB dataframe with DIST_NEAREST column added
distToNearest <- function(db, genotyped=F, grouping='first', model="HS5F") {
  	
	if(!is.data.frame(db))
		stop('Must submit either DB file or data frame')
	
	if(genotyped) { 
		v_col <- "V_CALL_GENOTYPED"
	} else {
		v_col <- "V_CALL"
	}
	j_col <- "J_CALL"
	
	# Parse V and J Column to get gene
	cat("V+J Column parsing\n")
	
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
	
	# Load targeting model
	cat("Loading Targeting Model\n")
	
	if(model=="HS5F"){
		load("HS5F_Targeting.RData", envir=.GlobalEnv)
	} else if(model==""){
		load("MRS5NF_Targeting.RData", envir=.GlobalEnv)
	} else if(model=="MTri"){
		load("MTri_Targeting.RData", envir=.GlobalEnv)
	} else{
		stop("Error: Unrecognized targeting model\n")
	}
	# Create new column for distance to nearest neighbor
	db$DIST_NEIGHBOR = rep(NA, nrow(db))
	
	# Create new column for parsed V and J
	junc_len <- "JUNCTION_GAP_LENGTH"
	db[,"V"] <- V
	db[,"J"] <- J
	db[,"ROW_ID"] <- 1:nrow(db)
	
	cat("Calculating distance to nearest neighbor\n")
	db <- ddply( ddply(db,c("V","J","JUNCTION_GAP_LENGTH"), transform, 
					"DIST_NEIGHBOR"=getDistanceToClosest(JUNCTION), .parallel=TRUE), 
			.(ROW_ID) , .parallel=TRUE)
	
	db <- db[,!(names(db) %in% c("V","J","ROW_ID"))]
	return(db)
}