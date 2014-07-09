#########################################################################################
# License Agreement
# 
# THIS WORK IS PROVIDED UNDER THE TERMS OF THIS CREATIVE COMMONS PUBLIC LICENSE 
# ("CCPL" OR "LICENSE"). THE WORK IS PROTECTED BY COPYRIGHT AND/OR OTHER 
# APPLICABLE LAW. ANY USE OF THE WORK OTHER THAN AS AUTHORIZED UNDER THIS LICENSE 
# OR COPYRIGHT LAW IS PROHIBITED.
# 
# BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE TO BE 
# BOUND BY THE TERMS OF THIS LICENSE. TO THE EXTENT THIS LICENSE MAY BE CONSIDERED 
# TO BE A CONTRACT, THE LICENSOR GRANTS YOU THE RIGHTS CONTAINED HERE IN 
# CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND CONDITIONS.
#
# BASELIne: Bayesian Estimation of Antigen-Driven Selection in Immunoglobulin Sequences
# Coded by: Mohamed Uduman & Gur Yaari
# Copyright 2012 Kleinstein Lab
# Version: 1.1 (7.03.2012)
#########################################################################################

# Global variables  
  
  FILTER_BY_MUTATIONS = 1000

  # Nucleotides
  NUCLEOTIDES = c("A","C","G","T")
  
  # Amino Acids
  AMINO_ACIDS <- c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G")
  names(AMINO_ACIDS) <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
  names(AMINO_ACIDS) <- names(AMINO_ACIDS)

  #Amino Acid Traits
  #"*" "A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y"
  #B = "Hydrophobic/Burried"  N = "Intermediate/Neutral"  S="Hydrophilic/Surface") 
  TRAITS_AMINO_ACIDS_CHOTHIA98 <- c("*","N","B","S","S","B","N","N","B","S","B","B","S","N","S","S","N","N","B","B","N")
  names(TRAITS_AMINO_ACIDS_CHOTHIA98) <- sort(unique(AMINO_ACIDS))
  TRAITS_AMINO_ACIDS <- array(NA,21)
  
  # Codon Table
  CODON_TABLE <- as.data.frame(matrix(NA,ncol=64,nrow=12))

  # Substitution Model: Smith DS et al. 1996
  substitution_Literature_Mouse <- matrix(c(0, 0.156222928, 0.601501588, 0.242275484, 0.172506739, 0, 0.241239892, 0.586253369, 0.54636291, 0.255795364, 0, 0.197841727, 0.290240811, 0.467680608, 0.24207858, 0),nrow=4,byrow=T,dimnames=list(NUCLEOTIDES,NUCLEOTIDES))
  substitution_Flu_Human <- matrix(c(0,0.2795596,0.5026927,0.2177477,0.1693210,0,0.3264723,0.5042067,0.4983549,0.3328321,0,0.1688130,0.2021079,0.4696077,0.3282844,0),4,4,byrow=T,dimnames=list(NUCLEOTIDES,NUCLEOTIDES))
  substitution_Flu25_Human <- matrix(c(0,0.2580641,0.5163685,0.2255674,0.1541125,0,0.3210224,0.5248651,0.5239281,0.3101292,0,0.1659427,0.1997207,0.4579444,0.3423350,0),4,4,byrow=T,dimnames=list(NUCLEOTIDES,NUCLEOTIDES))
  load("FiveS_Substitution.RData")

  # Mutability Models: Shapiro GS et al. 2002
  triMutability_Literature_Human <- matrix(c(0.24, 1.2, 0.96, 0.43, 2.14, 2, 1.11, 1.9, 0.85, 1.83, 2.36, 1.31, 0.82, 0.52, 0.89, 1.33, 1.4, 0.82, 1.83, 0.73, 1.83, 1.62, 1.53, 0.57, 0.92, 0.42, 0.42, 1.47, 3.44, 2.58, 1.18, 0.47, 0.39, 1.12, 1.8, 0.68, 0.47, 2.19, 2.35, 2.19, 1.05, 1.84, 1.26, 0.28, 0.98, 2.37, 0.66, 1.58, 0.67, 0.92, 1.76, 0.83, 0.97, 0.56, 0.75, 0.62, 2.26, 0.62, 0.74, 1.11, 1.16, 0.61, 0.88, 0.67, 0.37, 0.07, 1.08, 0.46, 0.31, 0.94, 0.62, 0.57, 0.29, NA, 1.44, 0.46, 0.69, 0.57, 0.24, 0.37, 1.1, 0.99, 1.39, 0.6, 2.26, 1.24, 1.36, 0.52, 0.33, 0.26, 1.25, 0.37, 0.58, 1.03, 1.2, 0.34, 0.49, 0.33, 2.62, 0.16, 0.4, 0.16, 0.35, 0.75, 1.85, 0.94, 1.61, 0.85, 2.09, 1.39, 0.3, 0.52, 1.33, 0.29, 0.51, 0.26, 0.51, 3.83, 2.01, 0.71, 0.58, 0.62, 1.07, 0.28, 1.2, 0.74, 0.25, 0.59, 1.09, 0.91, 1.36, 0.45, 2.89, 1.27, 3.7, 0.69, 0.28, 0.41, 1.17, 0.56, 0.93, 3.41, 1, 1, NA, 5.9, 0.74, 2.51, 2.24, 2.24, 1.95, 3.32, 2.34, 1.3, 2.3, 1, 0.66, 0.73, 0.93, 0.41, 0.65, 0.89, 0.65, 0.32, NA, 0.43, 0.85, 0.43, 0.31, 0.31, 0.23, 0.29, 0.57, 0.71, 0.48, 0.44, 0.76, 0.51, 1.7, 0.85, 0.74, 2.23, 2.08, 1.16, 0.51, 0.51, 1, 0.5, NA, NA, 0.71, 2.14), nrow=64,byrow=T)
  triMutability_Literature_Mouse <- matrix(c(1.31, 1.35, 1.42, 1.18, 2.02, 2.02, 1.02, 1.61, 1.99, 1.42, 2.01, 1.03, 2.02, 0.97, 0.53, 0.71, 1.19, 0.83, 0.96, 0.96, 0, 1.7, 2.22, 0.59, 1.24, 1.07, 0.51, 1.68, 3.36, 3.36, 1.14, 0.29, 0.33, 0.9, 1.11, 0.63, 1.08, 2.07, 2.27, 1.74, 0.22, 1.19, 2.37, 1.15, 1.15, 1.56, 0.81, 0.34, 0.87, 0.79, 2.13, 0.49, 0.85, 0.97, 0.36, 0.82, 0.66, 0.63, 1.15, 0.94, 0.85, 0.25, 0.93, 1.19, 0.4, 0.2, 0.44, 0.44, 0.88, 1.06, 0.77, 0.39, 0, 0, 0, 0, 0, 0, 0.43, 0.43, 0.86, 0.59, 0.59, 0, 1.18, 0.86, 2.9, 1.66, 0.4, 0.2, 1.54, 0.43, 0.69, 1.71, 0.68, 0.55, 0.91, 0.7, 1.71, 0.09, 0.27, 0.63, 0.2, 0.45, 1.01, 1.63, 0.96, 1.48, 2.18, 1.2, 1.31, 0.66, 2.13, 0.49, 0, 0, 0, 2.97, 2.8, 0.79, 0.4, 0.5, 0.4, 0.11, 1.68, 0.42, 0.13, 0.44, 0.93, 0.71, 1.11, 1.19, 2.71, 1.08, 3.43, 0.4, 0.67, 0.47, 1.02, 0.14, 1.56, 1.98, 0.53, 0.33, 0.63, 2.06, 1.77, 1.46, 3.74, 2.93, 2.1, 2.18, 0.78, 0.73, 2.93, 0.63, 0.57, 0.17, 0.85, 0.52, 0.31, 0.31, 0, 0, 0.51, 0.29, 0.83, 0.54, 0.28, 0.47, 0.9, 0.99, 1.24, 2.47, 0.73, 0.23, 1.13, 0.24, 2.12, 0.24, 0.33, 0.83, 1.41, 0.62, 0.28, 0.35, 0.77, 0.17, 0.72, 0.58, 0.45, 0.41), nrow=64,byrow=T)
  triMutability_Names <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")
  load("FiveS_Mutability.RData")

# Functions
  
  # Translate codon to amino acid
  translateCodonToAminoAcid<-function(Codon){
     return(AMINO_ACIDS[Codon])
  }

  # Translate amino acid to trait change
  translateAminoAcidToTraitChange<-function(AminoAcid){
     return(TRAITS_AMINO_ACIDS[AminoAcid])
  }
    
  # Initialize Amino Acid Trait Changes
  initializeTraitChange <- function(traitChangeModel=1,species=1,traitChangeFileName=NULL){
    if(!is.null(traitChangeFileName)){
      tryCatch(
          traitChange <- read.delim(traitChangeFileName,sep="\t",header=T)
          , error = function(ex){
            cat("Error|Error reading trait changes. Please check file name/path and format.\n")
            q()
          }
        )
    }else{
      traitChange <- TRAITS_AMINO_ACIDS_CHOTHIA98
    }
    TRAITS_AMINO_ACIDS <<- traitChange
 } 
  
  # Read in formatted nucleotide substitution matrix
  initializeSubstitutionMatrix <- function(substitutionModel,species,subsMatFileName=NULL){
    if(!is.null(subsMatFileName)){
      tryCatch(
          subsMat <- read.delim(subsMatFileName,sep="\t",header=T)
          , error = function(ex){
            cat("Error|Error reading substitution matrix. Please check file name/path and format.\n")
            q()
          }
        )
      if(sum(apply(subsMat,1,sum)==1)!=4) subsMat = t(apply(subsMat,1,function(x)x/sum(x)))
    }else{
      if(substitutionModel==1)subsMat <- substitution_Literature_Mouse
      if(substitutionModel==2)subsMat <- substitution_Flu_Human      
      if(substitutionModel==3)subsMat <- substitution_Flu25_Human      
       
    }

    if(substitutionModel==0){
      subsMat <- matrix(1,4,4)
      subsMat[,] = 1/3
      subsMat[1,1] = 0
      subsMat[2,2] = 0
      subsMat[3,3] = 0
      subsMat[4,4] = 0
    }
    
    
    NUCLEOTIDESN = c(NUCLEOTIDES,"N", "-")
    if(substitutionModel==5){
      subsMat <- FiveS_Substitution
      return(subsMat)
    }else{
      subsMat <- rbind(subsMat,rep(NA,4),rep(NA,4))
      return( matrix(data.matrix(subsMat),6,4,dimnames=list(NUCLEOTIDESN,NUCLEOTIDES) ) )
    }
  }

   
  # Read in formatted Mutability file
  initializeMutabilityMatrix <- function(mutabilityModel=1, species=1,mutabilityMatFileName=NULL){
    if(!is.null(mutabilityMatFileName)){
        tryCatch(
            mutabilityMat <- read.delim(mutabilityMatFileName,sep="\t",header=T)
            , error = function(ex){
              cat("Error|Error reading mutability matrix. Please check file name/path and format.\n")
              q()
            }
          )
    }else{
      mutabilityMat <- triMutability_Literature_Human
      if(species==2) mutabilityMat <- triMutability_Literature_Mouse
    }

  if(mutabilityModel==0){ mutabilityMat <- matrix(1,64,3)}
  
    if(mutabilityModel==5){
      mutabilityMat <- FiveS_Mutability
      return(mutabilityMat)
    }else{
      return( matrix( data.matrix(mutabilityMat), 64, 3, dimnames=list(triMutability_Names,1:3)) )
    }
  }

  # Read FASTA file formats
  # Modified from read.fasta from the seqinR package
  baseline.read.fasta <-
  function (file = system.file("sequences/sample.fasta", package = "seqinr"), 
      seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE, 
      set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, 
      strip.desc = FALSE,  sizeof.longlong = .Machine$sizeof.longlong, 
      endian = .Platform$endian, apply.mask = TRUE) 
  {
      seqtype <- match.arg(seqtype)
  
          lines <- readLines(file)
          
          if (legacy.mode) {
              comments <- grep("^;", lines)
              if (length(comments) > 0) 
                  lines <- lines[-comments]
          }
          
          
          ind_groups<-which(substr(lines, 1L, 3L) == ">>>")
          lines_mod<-lines
  
          if(!length(ind_groups)){
              lines_mod<-c(">>>All sequences combined",lines)            
          }
          
          ind_groups<-which(substr(lines_mod, 1L, 3L) == ">>>")
  
          lines <- array("BLA",dim=(length(ind_groups)+length(lines_mod)))
          id<-sapply(1:length(ind_groups),function(i)ind_groups[i]+i-1)+1
          lines[id] <- "THIS IS A FAKE SEQUENCE"
          lines[-id] <- lines_mod
          rm(lines_mod)
  
  		ind <- which(substr(lines, 1L, 1L) == ">")
          nseq <- length(ind)
          if (nseq == 0) {
               stop("no line starting with a > character found")
          }        
          start <- ind + 1
          end <- ind - 1
  
          while( any(which(ind%in%end)) ){
            ind=ind[-which(ind%in%end)]
            nseq <- length(ind)
            if (nseq == 0) {
                stop("no line starting with a > character found")
            }        
            start <- ind + 1
            end <- ind - 1        
          }
          
          end <- c(end[-1], length(lines))
          sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], collapse = ""))
          if (seqonly) 
              return(sequences)
          nomseq <- lapply(seq_len(nseq), function(i) {
          
              #firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
              substr(lines[ind[i]], 2, nchar(lines[ind[i]]))
          
          })
          if (seqtype == "DNA") {
              if (forceDNAtolower) {
                  sequences <- as.list(tolower(chartr(".","-",sequences)))
              }else{
                  sequences <- as.list(toupper(chartr(".","-",sequences)))
              }
          }
          if (as.string == FALSE) 
              sequences <- lapply(sequences, s2c)
          if (set.attributes) {
              for (i in seq_len(nseq)) {
                  Annot <- lines[ind[i]]
                  if (strip.desc) 
                    Annot <- substr(Annot, 2L, nchar(Annot))
                  attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                    Annot = Annot, class = switch(seqtype, AA = "SeqFastaAA", 
                      DNA = "SeqFastadna"))
              }
          }
          names(sequences) <- nomseq
          return(sequences)
  }

  
  # Replaces non FASTA characters in input files with N  
  replaceNonFASTAChars <-function(inSeq="ACGTN-AApA"){
    gsub('[^ACGTNacgt[:punct:]-[:punct:].]','N',inSeq,perl=TRUE)
  }    
  
  # Find the germlines in the FASTA list
  germlinesInFile <- function(seqIDs){
    firstChar = sapply(seqIDs,function(x){substr(x,1,1)})
    secondChar = sapply(seqIDs,function(x){substr(x,2,2)})
    return(firstChar==">" & secondChar!=">")
  }
  
  # Find the groups in the FASTA list
  groupsInFile <- function(seqIDs){
    sapply(seqIDs,function(x){substr(x,1,2)})==">>"
  }

  # In the process of finding germlines/groups, expand from the start to end of the group
  expandTillNext <- function(vecPosToID){    
    IDs = names(vecPosToID)
    posOfInterests =  which(vecPosToID)
  
    expandedID = rep(NA,length(IDs))
    expandedIDNames = gsub(">","",IDs[posOfInterests])
    startIndexes = c(1,posOfInterests[-1])
    stopIndexes = c(posOfInterests[-1]-1,length(IDs))
    expandedID  = unlist(sapply(1:length(startIndexes),function(i){
                                    rep(i,stopIndexes[i]-startIndexes[i]+1)
                                  }))
    names(expandedID) = unlist(sapply(1:length(startIndexes),function(i){
                                    rep(expandedIDNames[i],stopIndexes[i]-startIndexes[i]+1)
                                  }))  
    return(expandedID)                                                                                                  
  }
    
  # Process FASTA (list) to return a matrix[input, germline)
  processInputAdvanced <- function(inputFASTA){
  
    seqIDs = names(inputFASTA)
    numbSeqs = length(seqIDs)
    posGermlines1 = germlinesInFile(seqIDs)
    numbGermlines = sum(posGermlines1)
    posGroups1 = groupsInFile(seqIDs)
    numbGroups = sum(posGroups1)
    consDef = NA
    
    if(numbGermlines==0){
      posGermlines = 2
      numbGermlines = 1  
    }
  
      glPositionsSum = cumsum(posGermlines1)
      glPositions = table(glPositionsSum)
      #Find the position of the conservation row
      consDefPos = as.numeric(names(glPositions[names(glPositions)!=0 & glPositions==1]))+1  
    if( length(consDefPos)> 0 ){
      consDefID =  match(consDefPos, glPositionsSum) 
      #The coservation rows need to be pulled out and stores seperately 
      consDef =  inputFASTA[consDefID]
      inputFASTA =  inputFASTA[-consDefID]
  
      seqIDs = names(inputFASTA)
      numbSeqs = length(seqIDs)
      posGermlines1 = germlinesInFile(seqIDs)
      numbGermlines = sum(posGermlines1)
      posGroups1 = groupsInFile(seqIDs)
      numbGroups = sum(posGroups1)
      if(numbGermlines==0){
        posGermlines = 2
        numbGermlines = 1  
      }    
    }
    
    posGroups <- expandTillNext(posGroups1)
    posGermlines <- expandTillNext(posGermlines1)
    posGermlines[posGroups1] = 0
    names(posGermlines)[posGroups1] = names(posGroups)[posGroups1]
    posInput = rep(TRUE,numbSeqs)
    posInput[posGroups1 | posGermlines1] = FALSE
    
    matInput = matrix(NA, nrow=sum(posInput), ncol=2)
    rownames(matInput) = seqIDs[posInput]
    colnames(matInput) = c("Input","Germline")
    
    vecInputFASTA = unlist(inputFASTA)  
    matInput[,1] = vecInputFASTA[posInput]
    matInput[,2] = vecInputFASTA[ which( names(inputFASTA)%in%paste(">",names(posGermlines)[posInput],sep="") )[ posGermlines[posInput]] ]
    
    germlines = posGermlines[posInput]
    groups = posGroups[posInput]
    
    return( list("matInput"=matInput, "germlines"=germlines, "groups"=groups, "conservationDefinition"=consDef ))      
  }


  # Replace leading and trailing dashes in the sequence
  replaceLeadingTrailingDashes <- function(x,readEnd){
    iiGap = unlist(gregexpr("-",x[1]))
    ggGap = unlist(gregexpr("-",x[2]))  
    #posToChange = intersect(iiGap,ggGap)
    
    
    seqIn = replaceLeadingTrailingDashesHelper(x[1])
    seqGL = replaceLeadingTrailingDashesHelper(x[2])
    seqTemplate = rep('N',readEnd)
    seqIn <- c(seqIn,seqTemplate[(length(seqIn)+1):readEnd])
    seqGL <- c(seqGL,seqTemplate[(length(seqGL)+1):readEnd])
#    if(posToChange!=-1){
#      seqIn[posToChange] = "-"
#      seqGL[posToChange] = "-"
#    }
  
    seqIn = c2s(seqIn[1:readEnd])
    seqGL = c2s(seqGL[1:readEnd])
  
    lenGL = nchar(seqGL)
    if(lenGL<readEnd){
      seqGL = paste(seqGL,c2s(rep("N",readEnd-lenGL)),sep="")
    }
  
    lenInput = nchar(seqIn)
    if(lenInput<readEnd){
      seqIn = paste(seqIn,c2s(rep("N",readEnd-lenInput)),sep="")
    }    
    return( c(seqIn,seqGL) )
  }  

  replaceLeadingTrailingDashesHelper <- function(x){
    grepResults = gregexpr("-*",x)
    grepResultsPos = unlist(grepResults)
    grepResultsLen =  attr(grepResults[[1]],"match.length")   
    x = s2c(x)
    if(x[1]=="-"){
      x[1:grepResultsLen[1]] = "N"      
    }
    if(x[length(x)]=="-"){
      x[(length(x)-grepResultsLen[length(grepResultsLen)]+1):length(x)] = "N"      
    }
    return(x)
  }



  
  # Check sequences for indels
  checkForInDels <- function(matInputP){
    insPos <- checkInsertion(matInputP)
    delPos <- checkDeletions(matInputP)
    return(list("Insertions"=insPos, "Deletions"=delPos))
  }

  # Check sequences for insertions
  checkInsertion <- function(matInputP){
    insertionCheck = apply( matInputP,1, function(x){
                                          inputGaps <- as.vector( gregexpr("-",x[1])[[1]] )
                                          glGaps <- as.vector( gregexpr("-",x[2])[[1]] )                                          
                                          return( is.finite( match(FALSE, glGaps%in%inputGaps ) ) )
                                        })   
    return(as.vector(insertionCheck))
  }
  # Fix inserstions
  fixInsertions <- function(matInputP){
    insPos <- checkInsertion(matInputP)
    sapply((1:nrow(matInputP))[insPos],function(rowIndex){
                                                x <- matInputP[rowIndex,]
                                                inputGaps <- gregexpr("-",x[1])[[1]]
                                                glGaps <- gregexpr("-",x[2])[[1]]
                                                posInsertions <- glGaps[!(glGaps%in%inputGaps)]
                                                inputInsertionToN <- s2c(x[2])
                                                inputInsertionToN[posInsertions]!="-"
                                                inputInsertionToN[posInsertions] <- "N"
                                                inputInsertionToN <- c2s(inputInsertionToN)
                                                matInput[rowIndex,2] <<- inputInsertionToN 
                                              })                                                               
    return(insPos)
  } 
    
  # Check sequences for deletions
  checkDeletions <-function(matInputP){
    deletionCheck = apply( matInputP,1, function(x){
                                          inputGaps <- as.vector( gregexpr("-",x[1])[[1]] )
                                          glGaps <- as.vector( gregexpr("-",x[2])[[1]] )
                                          return( is.finite( match(FALSE, inputGaps%in%glGaps ) ) )
                                      })
    return(as.vector(deletionCheck))                                      
  }
  # Fix sequences with deletions
  fixDeletions <- function(matInputP){
    delPos <- checkDeletions(matInputP)    
    sapply((1:nrow(matInputP))[delPos],function(rowIndex){
                                                x <- matInputP[rowIndex,]
                                                inputGaps <- gregexpr("-",x[1])[[1]]
                                                glGaps <- gregexpr("-",x[2])[[1]]
                                                posDeletions <- inputGaps[!(inputGaps%in%glGaps)]
                                                inputDeletionToN <- s2c(x[1])
                                                inputDeletionToN[posDeletions] <- "N"
                                                inputDeletionToN <- c2s(inputDeletionToN)
                                                matInput[rowIndex,1] <<- inputDeletionToN 
                                              })                                                                   
    return(delPos)
  }  
    

  # Trim DNA sequence to the last codon
  trimToLastCodon <- function(seqToTrim){
    seqLen = nchar(seqToTrim)  
    trimmedSeq = s2c(seqToTrim)
    poi = seqLen
    tailLen = 0
    
    while(trimmedSeq[poi]=="-" || trimmedSeq[poi]=="."){
      tailLen = tailLen + 1
      poi = poi - 1   
    }
    
    trimmedSeq = c2s(trimmedSeq[1:(seqLen-tailLen)])
    seqLen = nchar(trimmedSeq)
    # Trim sequence to last codon
  	if( getCodonPos(seqLen)[3] > seqLen )
  	  trimmedSeq = substr(seqToTrim,1, ( (getCodonPos(seqLen)[1])-1 ) )
    
    return(trimmedSeq)
  }
  
  # Given a nuclotide position, returns the pos of the 3 nucs that made the codon
  # e.g. nuc 86 is part of nucs 85,86,87
  getCodonPos <- function(nucPos){
    codonNum =  (ceiling(nucPos/3))*3
    return( (codonNum-2):codonNum)
  }
  
  # Given a nuclotide position, returns the codon number
  # e.g. nuc 86  = codon 29
  getCodonNumb <- function(nucPos){
    return( ceiling(nucPos/3) )
  }
  
  # Given a codon, returns all the nuc positions that make the codon
  getCodonNucs <- function(codonNumb){
    getCodonPos(codonNumb*3)
  }  

  computeCodonTable <- function(testID=1){
                  
    if(testID<=4){    
      # Pre-compute every codons
      intCounter = 1
      for(pOne in NUCLEOTIDES){
        for(pTwo in NUCLEOTIDES){
          for(pThree in NUCLEOTIDES){
            codon = paste(pOne,pTwo,pThree,sep="")
            colnames(CODON_TABLE)[intCounter] =  codon
            intCounter = intCounter + 1
            CODON_TABLE[,codon] = mutationTypeOptimized(cbind(permutateAllCodon(codon),rep(codon,12)))
          }  
        }
      }
      chars = c("N","A","C","G","T", "-")
      for(a in chars){
        for(b in chars){
          for(c in chars){
            if(a=="N" | b=="N" | c=="N"){ 
              #cat(paste(a,b,c),sep="","\n") 
              CODON_TABLE[,paste(a,b,c,sep="")] = rep(NA,12)
            }
          }  
        }
      }
      
      chars = c("-","A","C","G","T")
      for(a in chars){
        for(b in chars){
          for(c in chars){
            if(a=="-" | b=="-" | c=="-"){ 
              #cat(paste(a,b,c),sep="","\n") 
              CODON_TABLE[,paste(a,b,c,sep="")] = rep(NA,12)
            }
          }  
        }
      }
      CODON_TABLE <<- as.matrix(CODON_TABLE)
    }
  }
  
  collapseClone <- function(vecInputSeqs,glSeq,readEnd,nonTerminalOnly=0){
  #print(length(vecInputSeqs))
    vecInputSeqs = unique(vecInputSeqs) 
    if(length(vecInputSeqs)==1){
      return( list( c(vecInputSeqs,glSeq), F) )
    }else{
      charInputSeqs <- sapply(vecInputSeqs, function(x){
                                              s2c(x)[1:readEnd]
                                            })
      charGLSeq <- s2c(glSeq)
      matClone <- sapply(1:readEnd, function(i){
                                            posNucs = unique(charInputSeqs[i,])
                                            posGL = charGLSeq[i]
                                            error = FALSE                                            
                                            if(posGL=="-" & sum(!(posNucs%in%c("-","N")))==0 ){
                                              return(c("-",error))
                                            }
                                            if(length(posNucs)==1)
                                              return(c(posNucs[1],error))
                                            else{
                                              if("N"%in%posNucs){
                                                error=TRUE
                                              }
                                              if(sum(!posNucs[posNucs!="N"]%in%posGL)==0){
                                                return( c(posGL,error) )  
                                              }else{
                                                #return( c(sample(posNucs[posNucs!="N"],1),error) )  
                                                if(nonTerminalOnly==0){
                                                  return( c(sample(charInputSeqs[i,charInputSeqs[i,]!="N" & charInputSeqs[i,]!=posGL],1),error) )  
                                                }else{
                                                  posNucs = charInputSeqs[i,charInputSeqs[i,]!="N" & charInputSeqs[i,]!=posGL]
                                                  posNucsTable = table(posNucs)
                                                  if(sum(posNucsTable>1)==0){
                                                    return( c(posGL,error) )
                                                  }else{
                                                    return( c(sample( posNucs[posNucs%in%names(posNucsTable)[posNucsTable>1]],1),error) )
                                                  }
                                                }
                                                
                                              }
                                            } 
                                          })
      
                                          
      #print(length(vecInputSeqs))                                        
      return(list(c(c2s(matClone[1,]),glSeq),"TRUE"%in%matClone[2,]))
    }
  }

  # Compute the expected for each sequence-germline pair
  getExpectedIndividual <- function(matInput){    
    facGL <- factor(matInput[,2])
    facLevels = levels(facGL)
    LisGLs_MutabilityU = mclapply(1:length(facLevels),  function(x){
                                                      computeMutabilities(facLevels[x])
                                                    })
    facIndex = match(facGL,facLevels)
    
    LisGLs_Mutability = mclapply(1:nrow(matInput),  function(x){
                                                      cInput = rep(NA,nchar(matInput[x,1]))
                                                      cInput[s2c(matInput[x,1])!="N"] = 1
                                                      LisGLs_MutabilityU[[facIndex[x]]] * cInput                                                   
                                                    })
                                                    
    LisGLs_Targeting =  mclapply(1:dim(matInput)[1],  function(x){
                                                      computeTargeting(matInput[x,2],LisGLs_Mutability[[x]])
                                                    })
                                                    
    LisGLs_MutationTypes  = mclapply(1:length(matInput[,2]),function(x){
                                                    #print(x)
                                                    computeMutationTypes(matInput[x,2])
                                                })
    
    LisGLs_Exp = mclapply(1:dim(matInput)[1],  function(x){
                                                  computeExpected(LisGLs_Targeting[[x]],LisGLs_MutationTypes[[x]])
                                                })
    
    ul_LisGLs_Exp =  unlist(LisGLs_Exp)                                            
    return(matrix(ul_LisGLs_Exp,ncol=4,nrow=(length(ul_LisGLs_Exp)/4),byrow=T))
  }

  # Compute mutabilities of sequence based on the tri-nucleotide model
  computeMutabilities <- function(paramSeq){
    seqLen = nchar(paramSeq)
    seqMutabilites = rep(NA,seqLen)
  
    gaplessSeq = gsub("-", "", paramSeq)
    gaplessSeqLen = nchar(gaplessSeq)
    gaplessSeqMutabilites = rep(NA,gaplessSeqLen)
    
    if(mutabilityModel!=5){
      pos<- 3:(gaplessSeqLen)
      subSeq =  substr(rep(gaplessSeq,gaplessSeqLen-2),(pos-2),(pos+2))    
      gaplessSeqMutabilites[pos] =      
        tapply( c(
                                        getMutability( substr(subSeq,1,3), 3) , 
                                        getMutability( substr(subSeq,2,4), 2), 
                                        getMutability( substr(subSeq,3,5), 1) 
                                        ),rep(1:(gaplessSeqLen-2),3),mean,na.rm=TRUE
                                      )
      #Pos 1
      subSeq =  substr(gaplessSeq,1,3)
      gaplessSeqMutabilites[1] =  getMutability(subSeq , 1)
      #Pos 2
      subSeq =  substr(gaplessSeq,1,4)
      gaplessSeqMutabilites[2] =  mean( c(
                                            getMutability( substr(subSeq,1,3), 2) , 
                                            getMutability( substr(subSeq,2,4), 1) 
                                          ),na.rm=T
                                      ) 
      seqMutabilites[which(s2c(paramSeq)!="-")]<- gaplessSeqMutabilites
      return(seqMutabilites)
    }else{
      
      pos<- 3:(gaplessSeqLen)
      subSeq =  substr(rep(gaplessSeq,gaplessSeqLen-2),(pos-2),(pos+2))    
      gaplessSeqMutabilites[pos] = sapply(subSeq,function(x){ getMutability5(x) }, simplify=T)
      seqMutabilites[which(s2c(paramSeq)!="-")]<- gaplessSeqMutabilites
      return(seqMutabilites)
    }

  }

  # Returns the mutability of a triplet at a given position
  getMutability <- function(codon, pos=1:3){
    triplets <- rownames(mutability)
    mutability[  match(codon,triplets) ,pos]
  }

  getMutability5 <- function(fivemer){
    return(mutability[fivemer])
  }

  # Returns the substitution probabilty
  getTransistionProb <- function(nuc){
    substitution[nuc,]
  }

  getTransistionProb5 <- function(fivemer){    
    if(any(which(fivemer==colnames(substitution)))){
      return(substitution[,fivemer])
    }else{
      return(array(NA,4))
    }
  }

  # Given a nuc, returns the other 3 nucs it can mutate to
  canMutateTo <- function(nuc){
    NUCLEOTIDES[- which(NUCLEOTIDES==nuc)]
  }
  
  # Given a nucleotide, returns the probabilty of other nucleotide it can mutate to 
  canMutateToProb <- function(nuc){
    substitution[nuc,canMutateTo(nuc)]
  }

  # Compute targeting, based on precomputed mutatbility & substitution  
  computeTargeting <- function(param_strSeq,param_vecMutabilities){

    if(substitutionModel!=5){
      vecSeq = s2c(param_strSeq)
      matTargeting = sapply( 1:length(vecSeq), function(x) { param_vecMutabilities[x] * getTransistionProb(vecSeq[x]) } )  
      #matTargeting = apply( rbind(vecSeq,param_vecMutabilities),2, function(x) { as.vector(as.numeric(x[2]) * getTransistionProb(x[1])) } )
      dimnames( matTargeting ) =  list(NUCLEOTIDES,1:(length(vecSeq))) 
      return (matTargeting)
    }else{
      
      seqLen = nchar(param_strSeq)
      seqsubstitution = matrix(NA,ncol=seqLen,nrow=4)
      paramSeq <- param_strSeq
      gaplessSeq = gsub("-", "", paramSeq)
      gaplessSeqLen = nchar(gaplessSeq)
      gaplessSeqSubstitution  = matrix(NA,ncol=gaplessSeqLen,nrow=4) 
      
      pos<- 3:(gaplessSeqLen)
      subSeq =  substr(rep(gaplessSeq,gaplessSeqLen-2),(pos-2),(pos+2))    
      gaplessSeqSubstitution[,pos] = sapply(subSeq,function(x){ getTransistionProb5(x) }, simplify=T)
      seqsubstitution[,which(s2c(paramSeq)!="-")]<- gaplessSeqSubstitution
      matTargeting <- param_vecMutabilities * seqsubstitution
      dimnames( matTargeting ) =  list(NUCLEOTIDES,1:(seqLen)) 
      return (matTargeting)      
    }
  }  

  # Compute the mutations types   
  computeMutationTypes <- function(param_strSeq){
  #cat(param_strSeq,"\n")
    #vecSeq = trimToLastCodon(param_strSeq)
    lenSeq = nchar(param_strSeq)
    vecCodons = sapply({1:(lenSeq/3)}*3-2,function(x){substr(param_strSeq,x,x+2)})
    matMutationTypes = matrix( unlist(CODON_TABLE[,vecCodons]) ,ncol=lenSeq,nrow=4, byrow=F)
    dimnames( matMutationTypes ) =  list(NUCLEOTIDES,1:(ncol(matMutationTypes)))
    return(matMutationTypes)   
  }  
  computeMutationTypesFast <- function(param_strSeq){
    matMutationTypes = matrix( CODON_TABLE[,param_strSeq] ,ncol=3,nrow=4, byrow=F)
    #dimnames( matMutationTypes ) =  list(NUCLEOTIDES,1:(length(vecSeq)))
    return(matMutationTypes)   
  }  
  mutationTypeOptimized <- function( matOfCodons ){
   apply( matOfCodons,1,function(x){ mutationType(x[2],x[1]) } ) 
  }  

  # Returns a vector of codons 1 mutation away from the given codon
  permutateAllCodon <- function(codon){
    cCodon = s2c(codon)
    matCodons = t(array(cCodon,dim=c(3,12)))
    matCodons[1:4,1] = NUCLEOTIDES
    matCodons[5:8,2] = NUCLEOTIDES
    matCodons[9:12,3] = NUCLEOTIDES
    apply(matCodons,1,c2s)
  }

  # Given two codons, tells you if the mutation is R or S (based on your definition)
  mutationType <- function(codonFrom,codonTo){
    if(testID==4){
      if( is.na(codonFrom) | is.na(codonTo) | is.na(translateCodonToAminoAcid(codonFrom)) | is.na(translateCodonToAminoAcid(codonTo)) ){
        return(NA)
      }else{
        mutationType = "S"
        if( translateAminoAcidToTraitChange(translateCodonToAminoAcid(codonFrom)) != translateAminoAcidToTraitChange(translateCodonToAminoAcid(codonTo)) ){
          mutationType = "R"                                                              
        }
        if(translateCodonToAminoAcid(codonTo)=="*" | translateCodonToAminoAcid(codonFrom)=="*"){
          mutationType = "Stop"
        }
        return(mutationType)
      }  
    }else if(testID==5){  
      if( is.na(codonFrom) | is.na(codonTo) | is.na(translateCodonToAminoAcid(codonFrom)) | is.na(translateCodonToAminoAcid(codonTo)) ){
        return(NA)
      }else{
        if(codonFrom==codonTo){
          mutationType = "S"
        }else{
          codonFrom = s2c(codonFrom)
          codonTo = s2c(codonTo)  
          mutationType = "Stop"
          nucOfI = codonFrom[which(codonTo!=codonFrom)]
          if(nucOfI=="C"){
            mutationType = "R"  
          }else if(nucOfI=="G"){
            mutationType = "S"
          }
        }
        return(mutationType)
      }
    }else{
      if( is.na(codonFrom) | is.na(codonTo) | is.na(translateCodonToAminoAcid(codonFrom)) | is.na(translateCodonToAminoAcid(codonTo)) ){
        return(NA)
      }else{
        mutationType = "S"
        if( translateCodonToAminoAcid(codonFrom) != translateCodonToAminoAcid(codonTo) ){
          mutationType = "R"                                                              
        }
        if(translateCodonToAminoAcid(codonTo)=="*" | translateCodonToAminoAcid(codonFrom)=="*"){
          mutationType = "Stop"
        }
        return(mutationType)
      }  
    }    
  }

  
  #given a mat of targeting & it's corresponding mutationtypes returns 
  #a vector of Exp_RCDR,Exp_SCDR,Exp_RFWR,Exp_RFWR
  computeExpected <- function(paramTargeting,paramMutationTypes){
    # Replacements
    RPos = which(paramMutationTypes=="R")  
      #FWR
      Exp_R_FWR = sum(paramTargeting[ RPos[which(FWR_Nuc_Mat[RPos]==T)] ],na.rm=T)
      #CDR
      Exp_R_CDR = sum(paramTargeting[ RPos[which(CDR_Nuc_Mat[RPos]==T)] ],na.rm=T)
    # Silents
    SPos = which(paramMutationTypes=="S")  
      #FWR
      Exp_S_FWR = sum(paramTargeting[ SPos[which(FWR_Nuc_Mat[SPos]==T)] ],na.rm=T)
      #CDR
      Exp_S_CDR = sum(paramTargeting[ SPos[which(CDR_Nuc_Mat[SPos]==T)] ],na.rm=T)
  
      return(c(Exp_R_CDR,Exp_S_CDR,Exp_R_FWR,Exp_S_FWR))
  }
  
  # Count the mutations in a sequence
  # each mutation is treated independently 
  analyzeMutations2NucUri_website <- function( rev_in_matrix ){
    paramGL = rev_in_matrix[2,]
    paramSeq = rev_in_matrix[1,]  
    
    #Fill seq with GL seq if gapped
    #if( any(paramSeq=="-") ){
    #  gapPos_Seq =  which(paramSeq=="-")
    #  gapPos_Seq_ToReplace = gapPos_Seq[paramGL[gapPos_Seq] != "-"]
    #  paramSeq[gapPos_Seq_ToReplace] =  paramGL[gapPos_Seq_ToReplace]
    #}
  
  
    #if( any(paramSeq=="N") ){
    #  gapPos_Seq =  which(paramSeq=="N")
    #  gapPos_Seq_ToReplace = gapPos_Seq[paramGL[gapPos_Seq] != "N"]
    #  paramSeq[gapPos_Seq_ToReplace] =  paramGL[gapPos_Seq_ToReplace]
    #}  
      
    analyzeMutations2NucUri(  matrix(c( paramGL, paramSeq  ),2,length(paramGL),byrow=T)  )
    
  }

  #1 = GL 
  #2 = Seq
  analyzeMutations2NucUri <- function( in_matrix=matrix(c(c("A","A","A","C","C","C"),c("A","G","G","C","C","A")),2,6,byrow=T) ){
    paramGL = in_matrix[2,]
    paramSeq = in_matrix[1,]
    paramSeqUri = paramGL
    #mutations = apply(rbind(paramGL,paramSeq), 2, function(x){!x[1]==x[2]})
    mutations_val = paramGL != paramSeq   
    if(any(mutations_val)){
      mutationPos = {1:length(mutations_val)}[mutations_val]  
      mutationPos = mutationPos[sapply(mutationPos, function(x){!any(paramSeq[getCodonPos(x)]=="N")})]
      length_mutations =length(mutationPos)
      mutationInfo = rep(NA,length_mutations)
      if(any(mutationPos)){  

        pos<- mutationPos
        pos_array<-array(sapply(pos,getCodonPos))
        codonGL =  paramGL[pos_array]
        
        codonSeq = sapply(pos,function(x){
                                  seqP = paramGL[getCodonPos(x)]
                                  muCodonPos = {x-1}%%3+1 
                                  seqP[muCodonPos] = paramSeq[x]
                                  return(seqP)
                                })      
        GLcodons =  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
        Seqcodons =   apply(codonSeq,2,c2s)
        mutationInfo = apply(rbind(GLcodons , Seqcodons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})     
        names(mutationInfo) = mutationPos
    }
    if(any(!is.na(mutationInfo))){
      return(mutationInfo[!is.na(mutationInfo)])    
    }else{
      return(NA)
    }
    
    
    }else{
      return (NA)
    }
  }
  
  processNucMutations2 <- function(mu){
    if(!is.na(mu)){
      #R
      if(any(mu=="R")){
        Rs = mu[mu=="R"]
        nucNumbs = as.numeric(names(Rs))
        R_CDR = sum(as.integer(CDR_Nuc[nucNumbs]),na.rm=T)
        R_FWR = sum(as.integer(FWR_Nuc[nucNumbs]),na.rm=T)      
      }else{
        R_CDR = 0
        R_FWR = 0
      }    
      
      #S
      if(any(mu=="S")){
        Ss = mu[mu=="S"]
        nucNumbs = as.numeric(names(Ss))
        S_CDR = sum(as.integer(CDR_Nuc[nucNumbs]),na.rm=T)
        S_FWR = sum(as.integer(FWR_Nuc[nucNumbs]),na.rm=T)      
      }else{
        S_CDR = 0
        S_FWR = 0
      }    
      
      
      retVec = c(R_CDR,S_CDR,R_FWR,S_FWR)
      retVec[is.na(retVec)]=0
      return(retVec)
    }else{
      return(rep(0,4))
    }
  }        
  
  
  ## Z-score Test
  computeZScore <- function(mat, test="Focused"){
    matRes <- matrix(NA,ncol=2,nrow=(nrow(mat)))
    if(test=="Focused"){
      #Z_Focused_CDR
      #P_Denom = sum( mat[1,c(5,6,8)], na.rm=T )
      P = apply(mat[,c(5,6,8)],1,function(x){(x[1]/sum(x))})
      R_mean = apply(cbind(mat[,c(1,2,4)],P),1,function(x){x[4]*(sum(x[1:3]))})
      R_sd=sqrt(R_mean*(1-P))
      matRes[,1] = (mat[,1]-R_mean)/R_sd
    
      #Z_Focused_FWR
      #P_Denom = sum( mat[1,c(7,6,8)], na.rm=T )
      P = apply(mat[,c(7,6,8)],1,function(x){(x[1]/sum(x))})
      R_mean = apply(cbind(mat[,c(3,2,4)],P),1,function(x){x[4]*(sum(x[1:3]))})
      R_sd=sqrt(R_mean*(1-P))
      matRes[,2] = (mat[,3]-R_mean)/R_sd
    }
  
    if(test=="Local"){
      #Z_Focused_CDR
      #P_Denom = sum( mat[1,c(5,6,8)], na.rm=T )
      P = apply(mat[,c(5,6)],1,function(x){(x[1]/sum(x))})
      R_mean = apply(cbind(mat[,c(1,2)],P),1,function(x){x[3]*(sum(x[1:2]))})
      R_sd=sqrt(R_mean*(1-P))
      matRes[,1] = (mat[,1]-R_mean)/R_sd
    
      #Z_Focused_FWR
      #P_Denom = sum( mat[1,c(7,6,8)], na.rm=T )
      P = apply(mat[,c(7,8)],1,function(x){(x[1]/sum(x))})
      R_mean = apply(cbind(mat[,c(3,4)],P),1,function(x){x[3]*(sum(x[1:2]))})
      R_sd=sqrt(R_mean*(1-P))
      matRes[,2] = (mat[,3]-R_mean)/R_sd
    }
    
    if(test=="Imbalanced"){
      #Z_Focused_CDR
      #P_Denom = sum( mat[1,c(5,6,8)], na.rm=T )
      P = apply(mat[,5:8],1,function(x){((x[1]+x[2])/sum(x))})
      R_mean = apply(cbind(mat[,1:4],P),1,function(x){x[5]*(sum(x[1:4]))})
      R_sd=sqrt(R_mean*(1-P))
      matRes[,1] = (mat[,1]-R_mean)/R_sd
    
      #Z_Focused_FWR
      #P_Denom = sum( mat[1,c(7,6,8)], na.rm=T )
      P = apply(mat[,5:8],1,function(x){((x[3]+x[4])/sum(x))})
      R_mean = apply(cbind(mat[,1:4],P),1,function(x){x[5]*(sum(x[1:4]))})
      R_sd=sqrt(R_mean*(1-P))
      matRes[,2] = (mat[,3]-R_mean)/R_sd
    }    
      
    matRes[is.nan(matRes)] = NA
    return(matRes)
  }

  # Return a p-value for a z-score
  z2p <- function(z){
    p=NA
    if( !is.nan(z) && !is.na(z)){   
      if(z>0){
        p = (1 - pnorm(z,0,1))
      } else if(z<0){
        p = (-1 * pnorm(z,0,1))
      } else{
        p = 0.5
      }
    }else{
      p = NA
    }
    return(p)
  }    
  
  
  ## Bayesian  Test

  # Fitted parameter for the bayesian framework
BAYESIAN_FITTED<-c(0.407277142798302, 0.554007336744485, 0.63777155771234, 0.693989162719009, 0.735450014674917, 0.767972534429806, 0.794557287143399, 0.816906816601605, 0.83606796225341, 0.852729446430296, 0.867370424541641, 0.880339760590323, 0.891900995024999, 0.902259181289864, 0.911577919359,0.919990301665853, 0.927606458124537, 0.934518806350661, 0.940805863754375, 0.946534836475715, 0.951763691199255, 0.95654428191308, 0.960920179487397, 0.964930893680829, 0.968611312149038, 0.971992459313836, 0.975102110004818, 0.977964943023096, 0.980603428208439, 0.983037660179428, 0.985285800977406, 0.987364285326685, 0.989288037855441, 0.991070478823525, 0.992723699729969, 0.994259575477392, 0.995687688867975, 0.997017365051493, 0.998257085153047, 0.999414558305388, 1.00049681357804, 1.00151036237481, 1.00246080204981, 1.00335370751909, 1.0041939329768, 1.0049859393417, 1.00573382091263, 1.00644127217376, 1.00711179729107, 1.00774845526417, 1.00835412715854, 1.00893143010366, 1.00948275846309, 1.01001030293661, 1.01051606798079, 1.01100188771288, 1.01146944044216, 1.01192026195449, 1.01235575766094, 1.01277721370986)
  CONST_i <- sort(c(((2^(seq(-39,0,length.out=201)))/2)[1:200],(c(0:11,13:99)+0.5)/100,1-(2^(seq(-39,0,length.out=201)))/2))
  
  # Given x, M & p, returns a pdf 
  calculate_bayes <- function ( x=3, N=10, p=0.33,
                                i=CONST_i,
                                max_sigma=20,length_sigma=4001
                              ){
    if(!0%in%N){
      G <- max(length(x),length(N),length(p))
      x=array(x,dim=G)
      N=array(N,dim=G)
      p=array(p,dim=G)
      sigma_s<-seq(-max_sigma,max_sigma,length.out=length_sigma)
      sigma_1<-log({i/{1-i}}/{p/{1-p}})
      index<-min(N,60)
      y<-dbeta(i,x+BAYESIAN_FITTED[index],N+BAYESIAN_FITTED[index]-x)*(1-p)*p*exp(sigma_1)/({1-p}^2+2*p*{1-p}*exp(sigma_1)+{p^2}*exp(2*sigma_1))
      if(!sum(is.na(y))){
        tmp<-approx(sigma_1,y,sigma_s)$y
        tmp/sum(tmp)/{2*max_sigma/{length_sigma-1}}
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }  
  # Given a mat of observed & expected, return a list of CDR & FWR pdf for selection
  computeBayesianScore <- function(mat, test="Focused", max_sigma=20,length_sigma=4001){
    flagOneSeq = F
    if(nrow(mat)==1){
      mat=rbind(mat,mat)
      flagOneSeq = T
    }
    if(test=="Focused"){
      #CDR
      P = c(apply(mat[,c(5,6,8)],1,function(x){(x[1]/sum(x))}),0.5)
      N = c(apply(mat[,c(1,2,4)],1,function(x){(sum(x))}),0)
      X = c(mat[,1],0)
      bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesCDR = bayesCDR[-length(bayesCDR)]
  
      #FWR
      P = c(apply(mat[,c(7,6,8)],1,function(x){(x[1]/sum(x))}),0.5)
      N = c(apply(mat[,c(3,2,4)],1,function(x){(sum(x))}),0)
      X = c(mat[,3],0)
      bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesFWR = bayesFWR[-length(bayesFWR)]     
    }
    
    if(test=="Local"){
      #CDR
      P = c(apply(mat[,c(5,6)],1,function(x){(x[1]/sum(x))}),0.5)
      N = c(apply(mat[,c(1,2)],1,function(x){(sum(x))}),0)
      X = c(mat[,1],0)
      bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesCDR = bayesCDR[-length(bayesCDR)]
  
      #FWR
      P = c(apply(mat[,c(7,8)],1,function(x){(x[1]/sum(x))}),0.5)
      N = c(apply(mat[,c(3,4)],1,function(x){(sum(x))}),0)
      X = c(mat[,3],0)
      bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesFWR = bayesFWR[-length(bayesFWR)]     
    } 
     
    if(test=="Imbalanced"){
      #CDR
      P = c(apply(mat[,c(5:8)],1,function(x){((x[1]+x[2])/sum(x))}),0.5)
      N = c(apply(mat[,c(1:4)],1,function(x){(sum(x))}),0)
      X = c(apply(mat[,c(1:2)],1,function(x){(sum(x))}),0)
      bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesCDR = bayesCDR[-length(bayesCDR)]
  
      #FWR
      P = c(apply(mat[,c(5:8)],1,function(x){((x[3]+x[4])/sum(x))}),0.5)
      N = c(apply(mat[,c(1:4)],1,function(x){(sum(x))}),0)
      X = c(apply(mat[,c(3:4)],1,function(x){(sum(x))}),0)
      bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesFWR = bayesFWR[-length(bayesFWR)]     
    }

    if(test=="ImbalancedSilent"){
      #CDR
      P = c(apply(mat[,c(6,8)],1,function(x){((x[1])/sum(x))}),0.5)
      N = c(apply(mat[,c(2,4)],1,function(x){(sum(x))}),0)
      X = c(apply(mat[,c(2,4)],1,function(x){(x[1])}),0)
      bayesCDR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesCDR = bayesCDR[-length(bayesCDR)]
  
      #FWR
      P = c(apply(mat[,c(6,8)],1,function(x){((x[2])/sum(x))}),0.5)
      N = c(apply(mat[,c(2,4)],1,function(x){(sum(x))}),0)
      X = c(apply(mat[,c(2,4)],1,function(x){(x[2])}),0)
      bayesFWR = apply(cbind(X,N,P),1,function(x){calculate_bayes(x=x[1],N=x[2],p=x[3],max_sigma=max_sigma,length_sigma=length_sigma)})    
      bayesFWR = bayesFWR[-length(bayesFWR)]     
    }
        
    if(flagOneSeq==T){
      bayesCDR = bayesCDR[1]  
      bayesFWR = bayesFWR[1]
    }
    return( list("CDR"=bayesCDR, "FWR"=bayesFWR) )
  }
  
  ##Covolution
  break2chunks<-function(G=1000){
  base<-2^round(log(sqrt(G),2),0)
  return(c(rep(base,floor(G/base)-1),base+G-(floor(G/base)*base)))
  }  
  
  PowersOfTwo <- function(G=100){
    exponents <- array()
    i = 0
    while(G > 0){
      i=i+1
      exponents[i] <- floor( log2(G) )
      G <- G-2^exponents[i]
    }
    return(exponents)
  }
  
  convolutionPowersOfTwo <- function( cons, length_sigma=4001 ){
    G = ncol(cons)
    if(G>1){
      for(gen in log(G,2):1){
        ll<-seq(from=2,to=2^gen,by=2)
        sapply(ll,function(l){cons[,l/2]<<-weighted_conv(cons[,l],cons[,l-1],length_sigma=length_sigma)})
      }
    }
    return( cons[,1] )
  }
  
  convolutionPowersOfTwoByTwos <- function( cons, length_sigma=4001,G=1 ){
    if(length(ncol(cons))) G<-ncol(cons)
    groups <- PowersOfTwo(G)
    matG <- matrix(NA, ncol=length(groups), nrow=length(cons)/G )
    startIndex = 1
    for( i in 1:length(groups) ){
      stopIndex <- 2^groups[i] + startIndex - 1
      if(stopIndex!=startIndex){
        matG[,i] <- convolutionPowersOfTwo( cons[,startIndex:stopIndex], length_sigma=length_sigma )
        startIndex = stopIndex + 1
      }
      else {
        if(G>1) matG[,i] <- cons[,startIndex:stopIndex]
        else matG[,i] <- cons
        #startIndex = stopIndex + 1
      }
    }
    return( list( matG, groups ) )
  }
  
  weighted_conv<-function(x,y,w=1,m=100,length_sigma=4001){
    lx<-length(x)
    ly<-length(y)
    if({lx<m}| {{lx*w}<m}| {{ly}<m}| {{ly*w}<m}){
      if(w<1){
        y1<-approx(1:ly,y,seq(1,ly,length.out=m))$y
        x1<-approx(1:lx,x,seq(1,lx,length.out=m/w))$y
        lx<-length(x1)
        ly<-length(y1)
      }
      else {
        y1<-approx(1:ly,y,seq(1,ly,length.out=m*w))$y
        x1<-approx(1:lx,x,seq(1,lx,length.out=m))$y
        lx<-length(x1)
        ly<-length(y1)
      }
    }
    else{
      x1<-x
      y1<-approx(1:ly,y,seq(1,ly,length.out=floor(lx*w)))$y
      ly<-length(y1)
    }
    tmp<-approx(x=1:(lx+ly-1),y=convolve(x1,rev(y1),type="open"),xout=seq(1,lx+ly-1,length.out=length_sigma))$y
    tmp[tmp<=0] = 0
    return(tmp/sum(tmp))
  }
  
  calculate_bayesGHelper <- function( listMatG,length_sigma=4001 ){
    matG <- listMatG[[1]]
    groups <- listMatG[[2]]
    i = 1
    resConv <- matG[,i]
    denom <- 2^groups[i]
    if(length(groups)>1){
      while( i<length(groups) ){
        i = i + 1
        resConv <- weighted_conv(resConv, matG[,i], w= {{2^groups[i]}/denom} ,length_sigma=length_sigma)
        #cat({{2^groups[i]}/denom},"\n")
        denom <- denom + 2^groups[i]
      }
    }
    return(resConv)
  }
  
  # Given a list of PDFs, returns a convoluted PDF    
  groupPosteriors <- function( listPosteriors, max_sigma=20, length_sigma=4001 ,Threshold=2 ){  
    listPosteriors = listPosteriors[ !is.na(listPosteriors) ]
    Length_Postrior<-length(listPosteriors)
    if(Length_Postrior>1 & Length_Postrior<=Threshold){
      cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
      listMatG <- convolutionPowersOfTwoByTwos(cons,length_sigma=length_sigma)
      y<-calculate_bayesGHelper(listMatG,length_sigma=length_sigma)
      return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }else if(Length_Postrior==1) return(listPosteriors[[1]])
    else  if(Length_Postrior==0) return(NA)
    else {
      cons = matrix(unlist(listPosteriors),length(listPosteriors[[1]]),length(listPosteriors))
      y = fastConv(cons,max_sigma=max_sigma, length_sigma=length_sigma )
      return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }
  }

  fastConv<-function(cons, max_sigma=20, length_sigma=4001){
    chunks<-break2chunks(G=ncol(cons))
    if(ncol(cons)==3) chunks<-2:1
    index_chunks_end <- cumsum(chunks)
    index_chunks_start <- c(1,index_chunks_end[-length(index_chunks_end)]+1)
    index_chunks <- cbind(index_chunks_start,index_chunks_end)
    
    case <- sum(chunks!=chunks[1])
    if(case==1) End <- max(1,((length(index_chunks)/2)-1))
    else End <- max(1,((length(index_chunks)/2)))
    
    firsts <- sapply(1:End,function(i){
          	    indexes<-index_chunks[i,1]:index_chunks[i,2]
          	    convolutionPowersOfTwoByTwos(cons[ ,indexes])[[1]]
          	  })
    if(case==0){
    	result<-calculate_bayesGHelper( convolutionPowersOfTwoByTwos(firsts) )
    }else if(case==1){
      last<-list(calculate_bayesGHelper(
      convolutionPowersOfTwoByTwos( cons[ ,index_chunks[length(index_chunks)/2,1]:index_chunks[length(index_chunks)/2,2]] )
                                      ),0)
      result_first<-calculate_bayesGHelper(convolutionPowersOfTwoByTwos(firsts))
      result<-calculate_bayesGHelper(
        list(
          cbind(
          result_first,last[[1]]),
          c(log(index_chunks_end[length(index_chunks)/2-1],2),log(index_chunks[length(index_chunks)/2,2]-index_chunks[length(index_chunks)/2,1]+1,2))
        )
      )
    }
    return(as.vector(result))
  }
    
  # Computes the 95% CI for a pdf
  calcBayesCI <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
    cdf = cumsum(Pdf)
    cdf = cdf/cdf[length(cdf)]  
    return( c(sigma_s[findInterval(low,cdf)-1] , sigma_s[findInterval(up,cdf)]) ) 
  }
  
  # Computes a mean for a pdf
  calcBayesMean <- function(Pdf,max_sigma=20,length_sigma=4001){
    if(length(Pdf)!=length_sigma) return(NA)
    sigma_s=seq(-max_sigma,max_sigma,length.out=length_sigma)
    norm = {length_sigma-1}/2/max_sigma
    return( (Pdf%*%sigma_s/norm)  ) 
  }
  
  # Returns the mean, and the 95% CI for a pdf
  calcBayesOutputInfo <- function(Pdf,low=0.025,up=0.975,max_sigma=20, length_sigma=4001){
    if(is.na(Pdf)) 
     return(rep(NA,3))  
    bCI = calcBayesCI(Pdf=Pdf,low=low,up=up,max_sigma=max_sigma,length_sigma=length_sigma)
    bMean = calcBayesMean(Pdf=Pdf,max_sigma=max_sigma,length_sigma=length_sigma)
    return(c(bMean, bCI))
  }   

  # Computes the p-value of a pdf
  computeSigmaP <- function(Pdf, length_sigma=4001, max_sigma=20){
    if(length(Pdf)>1){
      norm = {length_sigma-1}/2/max_sigma
      pVal = {sum(Pdf[1:{{length_sigma-1}/2}]) + Pdf[{{length_sigma+1}/2}]/2}/norm
      if(pVal>0.5){
        pVal = pVal-1
      }
      return(pVal)
    }else{
      return(NA)
    }
  }    
  
  # Compute p-value of two distributions
  compareTwoDistsFaster <-function(sigma_S=seq(-20,20,length.out=4001), N=10000, dens1=runif(4001,0,1), dens2=runif(4001,0,1)){
  #print(c(length(dens1),length(dens2)))
  if(length(dens1)>1 & length(dens2)>1 ){
    dens1<-dens1/sum(dens1)
    dens2<-dens2/sum(dens2)
    cum2 <- cumsum(dens2)-dens2/2
    tmp<- sum(sapply(1:length(dens1),function(i)return(dens1[i]*cum2[i])))
    #print(tmp)
    if(tmp>0.5)tmp<-tmp-1
    return( tmp )
    }
    else {
    return(NA)
    }
    #return (sum(sapply(1:N,function(i)(sample(sigma_S,1,prob=dens1)>sample(sigma_S,1,prob=dens2))))/N)
  }  
  
  # get number of seqeunces contributing to the sigma (i.e. seqeunces with mutations)
  numberOfSeqsWithMutations <- function(matMutations,test=1){
    cdrSeqs <- 0
    fwrSeqs <- 0    
    if(test==1){#focused
      cdrMutations <- apply(matMutations, 1, function(x){ sum(x[c(1,2,4)]) })
      fwrMutations <- apply(matMutations, 1, function(x){ sum(x[c(3,4,2)]) })
      if( any(which(cdrMutations>0)) ) cdrSeqs <- sum(cdrMutations>0)
      if( any(which(fwrMutations>0)) ) fwrSeqs <- sum(fwrMutations>0) 
    }
    if(test==2){#local
      cdrMutations <- apply(matMutations, 1, function(x){ sum(x[c(1,2)]) })
      fwrMutations <- apply(matMutations, 1, function(x){ sum(x[c(3,4)]) })
      if( any(which(cdrMutations>0)) ) cdrSeqs <- sum(cdrMutations>0)
      if( any(which(fwrMutations>0)) ) fwrSeqs <- sum(fwrMutations>0) 
    }
  return(c("CDR"=cdrSeqs, "FWR"=fwrSeqs))
}  



shadeColor <- function(sigmaVal=NA,pVal=NA){
  if(is.na(sigmaVal) & is.na(pVal)) return(NA)
  if(is.na(sigmaVal) & !is.na(pVal)) sigmaVal=sign(pVal)
  if(is.na(pVal) || pVal==1 || pVal==0){
    returnColor = "#FFFFFF";
  }else{
    colVal=abs(pVal);
    
    if(sigmaVal<0){      
        if(colVal>0.1)
          returnColor = "#CCFFCC";
        if(colVal<=0.1)
          returnColor = "#99FF99";
        if(colVal<=0.050)
          returnColor = "#66FF66";
        if(colVal<=0.010)
          returnColor = "#33FF33";
        if(colVal<=0.005)
          returnColor = "#00FF00";
      
    }else{
      if(colVal>0.1)
        returnColor = "#FFCCCC";
      if(colVal<=0.1)
        returnColor = "#FF9999";
      if(colVal<=0.05)
        returnColor = "#FF6666";
      if(colVal<=0.01)
        returnColor = "#FF3333";
      if(colVal<0.005)
        returnColor = "#FF0000";
    }
  }
  
  return(returnColor)
}



plotHelp <- function(xfrac=0.05,yfrac=0.05,log=FALSE){
  if(!log){
    x = par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac
    y = par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac
  }else {
    if(log==2){
      x = par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac
      y = 10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac)
    }
    if(log==1){
      x = 10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac)
      y = par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac
    }
    if(log==3){
      x = 10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac)
      y = 10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac)
    }
  }
  return(c("x"=x,"y"=y))
}

# SHMulation

  # Based on targeting, introduce a single mutation & then update the targeting 
  oneMutation <- function(){
    # Pick a postion + mutation
    posMutation = sample(1:(seqGermlineLen*4),1,replace=F,prob=as.vector(seqTargeting))
    posNucNumb = ceiling(posMutation/4)                    # Nucleotide number
    posNucKind = 4 - ( (posNucNumb*4) - posMutation )   # Nuc the position mutates to
  
    #mutate the simulation sequence
    seqSimVec <-  s2c(seqSim)
    seqSimVec[posNucNumb] <- NUCLEOTIDES[posNucKind]
    seqSim <<-  c2s(seqSimVec)
    
    #update Mutability, Targeting & MutationsTypes
    updateMutabilityNTargeting(posNucNumb)
  
    #return(c(posNucNumb,NUCLEOTIDES[posNucKind])) 
    return(posNucNumb)
  }  
  
  updateMutabilityNTargeting <- function(position){
    min_i<-max((position-2),1)
    max_i<-min((position+2),nchar(seqSim))
    min_ii<-min(min_i,3)
    
    #mutability - update locally
    seqMutability[(min_i):(max_i)] <<- computeMutabilities(substr(seqSim,position-4,position+4))[(min_ii):(max_i-min_i+min_ii)]
    
    
    #targeting - compute locally
    seqTargeting[,min_i:max_i] <<- computeTargeting(substr(seqSim,min_i,max_i),seqMutability[min_i:max_i])                 
    seqTargeting[is.na(seqTargeting)] <<- 0
    #mutCodonPos = getCodonPos(position) 
    mutCodonPos = seq(getCodonPos(min_i)[1],getCodonPos(max_i)[3])
    #cat(mutCodonPos,"\n")                                                  
    mutTypeCodon = getCodonPos(position)
    seqMutationTypes[,mutTypeCodon] <<- computeMutationTypesFast( substr(seqSim,mutTypeCodon[1],mutTypeCodon[3]) ) 
    # Stop = 0
    if(any(seqMutationTypes[,mutCodonPos]=="Stop",na.rm=T )){
      seqTargeting[,mutCodonPos][seqMutationTypes[,mutCodonPos]=="Stop"] <<- 0
    }
    
  
    #Selection
    selectedPos = (min_i*4-4)+(which(seqMutationTypes[,min_i:max_i]=="R"))  
    # CDR
    selectedCDR = selectedPos[which(matCDR[selectedPos]==T)]
    seqTargeting[selectedCDR] <<-  seqTargeting[selectedCDR] *  exp(selCDR)
    seqTargeting[selectedCDR] <<- seqTargeting[selectedCDR]/baseLineCDR_K
        
    # FWR
    selectedFWR = selectedPos[which(matFWR[selectedPos]==T)]
    seqTargeting[selectedFWR] <<-  seqTargeting[selectedFWR] *  exp(selFWR)
    seqTargeting[selectedFWR] <<- seqTargeting[selectedFWR]/baseLineFWR_K      
    
  }  
  


  # Validate the mutation: if the mutation has not been sampled before validate it, else discard it.   
  validateMutation <- function(){  
    if( !(mutatedPos%in%mutatedPositions) ){ # if it's a new mutation
      uniqueMutationsIntroduced <<- uniqueMutationsIntroduced + 1
      mutatedPositions[uniqueMutationsIntroduced] <<-  mutatedPos  
    }else{
      if(substr(seqSim,mutatedPos,mutatedPos)==substr(seqGermline,mutatedPos,mutatedPos)){ # back to germline mutation
        mutatedPositions <<-  mutatedPositions[-which(mutatedPositions==mutatedPos)]
        uniqueMutationsIntroduced <<-  uniqueMutationsIntroduced - 1
      }      
    }
  }  
  
  
  
  # Places text (labels) at normalized coordinates 
  myaxis <- function(xfrac=0.05,yfrac=0.05,log=FALSE,w="text",cex=1,adj=1,thecol="black"){
    par(xpd=TRUE)
    if(!log)
      text(par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac,par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac,w,cex=cex,adj=adj,col=thecol)
    else {
    if(log==2)
    text(
      par()$usr[1]-(par()$usr[2]-par()$usr[1])*xfrac,
      10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac),
      w,cex=cex,adj=adj,col=thecol)
    if(log==1)
      text(
      10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac),
      par()$usr[4]+(par()$usr[4]-par()$usr[3])*yfrac,
      w,cex=cex,adj=adj,col=thecol)
    if(log==3)
      text(
      10^((par()$usr[1])-((par()$usr[2])-(par()$usr[1]))*xfrac),
      10^((par()$usr[4])+((par()$usr[4])-(par()$usr[3]))*yfrac),
      w,cex=cex,adj=adj,col=thecol)
    }
    par(xpd=FALSE)
  }
  
  
  
  # Count the mutations in a sequence
  analyzeMutations <- function( inputMatrixIndex, model = 0 , multipleMutation=0, seqWithStops=0){

    paramGL = s2c(matInput[inputMatrixIndex,2])
    paramSeq = s2c(matInput[inputMatrixIndex,1])            
    
    #if( any(paramSeq=="N") ){
    #  gapPos_Seq =  which(paramSeq=="N")
    #  gapPos_Seq_ToReplace = gapPos_Seq[paramGL[gapPos_Seq] != "N"]
    #  paramSeq[gapPos_Seq_ToReplace] =  paramGL[gapPos_Seq_ToReplace]
    #}        
    mutations_val = paramGL != paramSeq   
    
    if(any(mutations_val)){
      mutationPos = which(mutations_val)#{1:length(mutations_val)}[mutations_val]  
      length_mutations =length(mutationPos)
      mutationInfo = rep(NA,length_mutations)
                          
      pos<- mutationPos
      pos_array<-array(sapply(pos,getCodonPos))
      codonGL =  paramGL[pos_array]
      codonSeqWhole =  paramSeq[pos_array]
      codonSeq = sapply(pos,function(x){
                                seqP = paramGL[getCodonPos(x)]
                                muCodonPos = {x-1}%%3+1 
                                seqP[muCodonPos] = paramSeq[x]
                                return(seqP)
                              })
      GLcodons =  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
      SeqcodonsWhole =  apply(matrix(codonSeqWhole,length_mutations,3,byrow=TRUE),1,c2s)      
      Seqcodons =   apply(codonSeq,2,c2s)
      
      mutationInfo = apply(rbind(GLcodons , Seqcodons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})     
      names(mutationInfo) = mutationPos     
      
      mutationInfoWhole = apply(rbind(GLcodons , SeqcodonsWhole),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})           
      names(mutationInfoWhole) = mutationPos

      mutationInfo <- mutationInfo[!is.na(mutationInfo)]
      mutationInfoWhole <- mutationInfoWhole[!is.na(mutationInfoWhole)]
      
      if(any(!is.na(mutationInfo))){       
  
        #Filter based on Stop (at the codon level)
        if(seqWithStops==1){
          nucleotidesAtStopCodons = names(mutationInfoWhole[mutationInfoWhole!="Stop"])
          mutationInfo = mutationInfo[nucleotidesAtStopCodons]
          mutationInfoWhole = mutationInfo[nucleotidesAtStopCodons]
        }else{
          countStops = sum(mutationInfoWhole=="Stop")
          if(seqWithStops==2 & countStops==0) mutationInfo = NA
          if(seqWithStops==3 & countStops>0) mutationInfo = NA
        }         
        
        if(any(!is.na(mutationInfo))){
          #Filter mutations based on multipleMutation
          if(multipleMutation==1 & !is.na(mutationInfo)){
            mutationCodons = getCodonNumb(as.numeric(names(mutationInfoWhole)))
            tableMutationCodons <- table(mutationCodons)
            codonsWithMultipleMutations <- as.numeric(names(tableMutationCodons[tableMutationCodons>1]))
            if(any(codonsWithMultipleMutations)){
              #remove the nucleotide mutations in the codons with multiple mutations
              mutationInfo <- mutationInfo[!(mutationCodons %in% codonsWithMultipleMutations)]
              #replace those codons with Ns in the input sequence
              paramSeq[unlist(lapply(codonsWithMultipleMutations, getCodonNucs))] = "N"
              matInput[inputMatrixIndex,1] <<- c2s(paramSeq)
            }
          }

          #Filter mutations based on the model
          if(any(mutationInfo)==T | is.na(any(mutationInfo))){        
            
            if(model==1 & !is.na(mutationInfo)){
              mutationInfo <- mutationInfo[mutationInfo=="S"]
            }  
            if(any(mutationInfo)==T | is.na(any(mutationInfo))) return(mutationInfo)
            else return(NA)
          }else{
            return(NA)
          }
        }else{
          return(NA)
        }
        
        
      }else{
        return(NA)
      }
    
    
    }else{
      return (NA)
    }    
  }  

   analyzeMutationsFixed <- function( inputArray, model = 0 , multipleMutation=0, seqWithStops=0){

    paramGL = s2c(inputArray[2])
    paramSeq = s2c(inputArray[1])            
    inputSeq <- inputArray[1]
    #if( any(paramSeq=="N") ){
    #  gapPos_Seq =  which(paramSeq=="N")
    #  gapPos_Seq_ToReplace = gapPos_Seq[paramGL[gapPos_Seq] != "N"]
    #  paramSeq[gapPos_Seq_ToReplace] =  paramGL[gapPos_Seq_ToReplace]
    #}        
    mutations_val = paramGL != paramSeq   
    
    if(any(mutations_val)){
      mutationPos = which(mutations_val)#{1:length(mutations_val)}[mutations_val]  
      length_mutations =length(mutationPos)
      mutationInfo = rep(NA,length_mutations)
                          
      pos<- mutationPos
      pos_array<-array(sapply(pos,getCodonPos))
      codonGL =  paramGL[pos_array]
      codonSeqWhole =  paramSeq[pos_array]
      codonSeq = sapply(pos,function(x){
                                seqP = paramGL[getCodonPos(x)]
                                muCodonPos = {x-1}%%3+1 
                                seqP[muCodonPos] = paramSeq[x]
                                return(seqP)
                              })
      GLcodons =  apply(matrix(codonGL,length_mutations,3,byrow=TRUE),1,c2s)
      SeqcodonsWhole =  apply(matrix(codonSeqWhole,length_mutations,3,byrow=TRUE),1,c2s)      
      Seqcodons =   apply(codonSeq,2,c2s)
      
      mutationInfo = apply(rbind(GLcodons , Seqcodons),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})     
      names(mutationInfo) = mutationPos     
      
      mutationInfoWhole = apply(rbind(GLcodons , SeqcodonsWhole),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})           
      names(mutationInfoWhole) = mutationPos

      mutationInfo <- mutationInfo[!is.na(mutationInfo)]
      mutationInfoWhole <- mutationInfoWhole[!is.na(mutationInfoWhole)]
      
      if(any(!is.na(mutationInfo))){       
  
        #Filter based on Stop (at the codon level)
        if(seqWithStops==1){
          nucleotidesAtStopCodons = names(mutationInfoWhole[mutationInfoWhole!="Stop"])
          mutationInfo = mutationInfo[nucleotidesAtStopCodons]
          mutationInfoWhole = mutationInfo[nucleotidesAtStopCodons]
        }else{
          countStops = sum(mutationInfoWhole=="Stop")
          if(seqWithStops==2 & countStops==0) mutationInfo = NA
          if(seqWithStops==3 & countStops>0) mutationInfo = NA
        }         
        
        if(any(!is.na(mutationInfo))){
          #Filter mutations based on multipleMutation
          if(multipleMutation==1 & !is.na(mutationInfo)){
            mutationCodons = getCodonNumb(as.numeric(names(mutationInfoWhole)))
            tableMutationCodons <- table(mutationCodons)
            codonsWithMultipleMutations <- as.numeric(names(tableMutationCodons[tableMutationCodons>1]))
            if(any(codonsWithMultipleMutations)){
              #remove the nucleotide mutations in the codons with multiple mutations
              mutationInfo <- mutationInfo[!(mutationCodons %in% codonsWithMultipleMutations)]
              #replace those codons with Ns in the input sequence
              paramSeq[unlist(lapply(codonsWithMultipleMutations, getCodonNucs))] = "N"
              #matInput[inputMatrixIndex,1] <<- c2s(paramSeq)
              inputSeq <- c2s(paramSeq)
            }
          }
          
          #Filter mutations based on the model
          if(any(mutationInfo)==T | is.na(any(mutationInfo))){        
            
            if(model==1 & !is.na(mutationInfo)){
              mutationInfo <- mutationInfo[mutationInfo=="S"]
            }  
            if(any(mutationInfo)==T | is.na(any(mutationInfo))) return(list(mutationInfo,inputSeq))
            else return(list(NA,inputSeq))
          }else{
            return(list(NA,inputSeq))
          }
        }else{
          return(list(NA,inputSeq))
        }
        
        
      }else{
        return(list(NA,inputSeq))
      }
    
    
    }else{
      return (list(NA,inputSeq))
    }    
  }  
 
  # triMutability Background Count
  buildMutabilityModel <- function( inputMatrixIndex, model=0 , multipleMutation=0, seqWithStops=0, stopMutations=0){
    
    #rowOrigMatInput = matInput[inputMatrixIndex,]    
    seqGL =  gsub("-", "", matInput[inputMatrixIndex,2])
    seqInput = gsub("-", "", matInput[inputMatrixIndex,1])    
    #matInput[inputMatrixIndex,] <<- cbind(seqInput,seqGL)
    tempInput <- cbind(seqInput,seqGL)
    seqLength = nchar(seqGL)      
    list_analyzeMutationsFixed<- analyzeMutationsFixed(tempInput, model, multipleMutation, seqWithStops)
    mutationCount <- list_analyzeMutationsFixed[[1]]
    seqInput <- list_analyzeMutationsFixed[[2]]
    BackgroundMatrix = mutabilityMatrix
    MutationMatrix = mutabilityMatrix    
    MutationCountMatrix = mutabilityMatrix    
    if(!is.na(mutationCount)){
      if((stopMutations==0 & model==0) | (stopMutations==1 & (sum(mutationCount=="Stop")<length(mutationCount))) | (model==1 & (sum(mutationCount=="S")>0)) ){ 
                  
        fivermerStartPos = 1:(seqLength-4)
        fivemerLength <- length(fivermerStartPos)
        fivemerGL <- substr(rep(seqGL,length(fivermerStartPos)),(fivermerStartPos),(fivermerStartPos+4))
        fivemerSeq <- substr(rep(seqInput,length(fivermerStartPos)),(fivermerStartPos),(fivermerStartPos+4))
    
        #Background
        for(fivemerIndex in 1:fivemerLength){
          fivemer = fivemerGL[fivemerIndex]
          if(!any(grep("N",fivemer))){
            fivemerCodonPos = fivemerCodon(fivemerIndex)
            fivemerReadingFrameCodon = substr(fivemer,fivemerCodonPos[1],fivemerCodonPos[3]) 
            fivemerReadingFrameCodonInputSeq = substr(fivemerSeq[fivemerIndex],fivemerCodonPos[1],fivemerCodonPos[3])          
            
            # All mutations model
            #if(!any(grep("N",fivemerReadingFrameCodon))){
              if(model==0){
                if(stopMutations==0){
                  if(!any(grep("N",fivemerReadingFrameCodonInputSeq)))
                    BackgroundMatrix[fivemer] <- (BackgroundMatrix[fivemer] + 1)              
                }else{
                  if( !any(grep("N",fivemerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(fivemerReadingFrameCodon)!="*" ){
                    positionWithinCodon = which(fivemerCodonPos==3)#positionsWithinCodon[(fivemerCodonPos[1]%%3)+1]
                    BackgroundMatrix[fivemer] <- (BackgroundMatrix[fivemer] + probNonStopMutations[fivemerReadingFrameCodon,positionWithinCodon])
                  }
                }
              }else{ # Only silent mutations
                if( !any(grep("N",fivemerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(fivemerReadingFrameCodon)!="*" & translateCodonToAminoAcid(fivemerReadingFrameCodonInputSeq)==translateCodonToAminoAcid(fivemerReadingFrameCodon)){
                  positionWithinCodon = which(fivemerCodonPos==3)
                  BackgroundMatrix[fivemer] <- (BackgroundMatrix[fivemer] + probSMutations[fivemerReadingFrameCodon,positionWithinCodon])
                }
              }
            #}
          }
        }
        
        #Mutations
        if(stopMutations==1) mutationCount = mutationCount[mutationCount!="Stop"]
        if(model==1) mutationCount = mutationCount[mutationCount=="S"]  
        mutationPositions = as.numeric(names(mutationCount))
        mutationCount = mutationCount[mutationPositions>2 & mutationPositions<(seqLength-1)]
        mutationPositions =  mutationPositions[mutationPositions>2 & mutationPositions<(seqLength-1)]
        countMutations = 0 
        for(mutationPosition in mutationPositions){
          fivemerIndex = mutationPosition-2
          fivemer = fivemerSeq[fivemerIndex]
          GLfivemer = fivemerGL[fivemerIndex]
          fivemerCodonPos = fivemerCodon(fivemerIndex)
          fivemerReadingFrameCodon = substr(fivemer,fivemerCodonPos[1],fivemerCodonPos[3]) 
          fivemerReadingFrameCodonGL = substr(GLfivemer,fivemerCodonPos[1],fivemerCodonPos[3])
          if(!any(grep("N",fivemer)) & !any(grep("N",GLfivemer))){
            if(model==0){
                countMutations = countMutations + 1              
                MutationMatrix[GLfivemer] <- (MutationMatrix[GLfivemer] + 1)
                MutationCountMatrix[GLfivemer] <- (MutationCountMatrix[GLfivemer] + 1)             
            }else{
              if( translateCodonToAminoAcid(fivemerReadingFrameCodonGL)!="*" ){
                  countMutations = countMutations + 1
                  positionWithinCodon = which(fivemerCodonPos==3)
                  glNuc =  substr(fivemerReadingFrameCodonGL,positionWithinCodon,positionWithinCodon)
                  inputNuc =  substr(fivemerReadingFrameCodon,positionWithinCodon,positionWithinCodon)
                  MutationMatrix[GLfivemer] <- (MutationMatrix[GLfivemer] + substitution[glNuc,inputNuc])
                  MutationCountMatrix[GLfivemer] <- (MutationCountMatrix[GLfivemer] + 1)                                    
              }                
            }                  
          }              
        }
        
        seqMutability = MutationMatrix/BackgroundMatrix
        seqMutability = seqMutability/sum(seqMutability,na.rm=TRUE)
        #cat(inputMatrixIndex,"\t",countMutations,"\n")
        return(list("seqMutability"  = seqMutability,"numbMutations" = countMutations,"seqMutabilityCount" = MutationCountMatrix, "BackgroundMatrix"=BackgroundMatrix))      
        
      }        
    }
  
  }  
  
  #Returns the codon position containing the middle nucleotide
  fivemerCodon <- function(fivemerIndex){
    codonPos = list(2:4,1:3,3:5)
    fivemerType = fivemerIndex%%3
    return(codonPos[[fivemerType+1]])
  }

  #returns probability values for one mutation in codons resulting in R, S or Stop
  probMutations <- function(typeOfMutation){    
    matMutationProb <- matrix(0,ncol=3,nrow=125,dimnames=list(words(alphabet = c(NUCLEOTIDES,"N"), length=3),c(1:3)))   
    for(codon in rownames(matMutationProb)){
        if( !any(grep("N",codon)) ){
        for(muPos in 1:3){
          matCodon = matrix(rep(s2c(codon),3),nrow=3,ncol=3,byrow=T)
          glNuc = matCodon[1,muPos]
          matCodon[,muPos] = canMutateTo(glNuc) 
          substitutionRate = substitution[glNuc,matCodon[,muPos]]
          typeOfMutations = apply(rbind(rep(codon,3),apply(matCodon,1,c2s)),2,function(x){mutationType(c2s(x[1]),c2s(x[2]))})        
          matMutationProb[codon,muPos] <- sum(substitutionRate[typeOfMutations==typeOfMutation])
        }
      }
    }
    
    return(matMutationProb) 
  }
  
  
  
  
#Mapping Trinucleotides to fivemers
mapTriToFivemer <- function(triMutability=triMutability_Literature_Human){
  rownames(triMutability) <- triMutability_Names
  Fivemer<-rep(NA,1024)
  names(Fivemer)<-words(alphabet=NUCLEOTIDES,length=5)
  Fivemer<-sapply(names(Fivemer),function(Word)return(sum( c(triMutability[substring(Word,3,5),1],triMutability[substring(Word,2,4),2],triMutability[substring(Word,1,3),3]),na.rm=TRUE)))
  Fivemer<-Fivemer/sum(Fivemer)
  return(Fivemer)
}

collapseFivemerToTri<-function(Fivemer,Weights=MutabilityWeights,position=1,NUC="A"){
  Indices<-substring(names(Fivemer),3,3)==NUC
  Factors<-substring(names(Fivemer[Indices]),(4-position),(6-position))
  tapply(which(Indices),Factors,function(i)weighted.mean(Fivemer[i],Weights[i],na.rm=TRUE))
}



CountFivemerToTri<-function(Fivemer,Weights=MutabilityWeights,position=1,NUC="A"){
  Indices<-substring(names(Fivemer),3,3)==NUC
  Factors<-substring(names(Fivemer[Indices]),(4-position),(6-position))
  tapply(which(Indices),Factors,function(i)sum(Weights[i],na.rm=TRUE))
}

#Uses the real counts of the mutated fivemers
CountFivemerToTri2<-function(Fivemer,Counts=MutabilityCounts,position=1,NUC="A"){
  Indices<-substring(names(Fivemer),3,3)==NUC
  Factors<-substring(names(Fivemer[Indices]),(4-position),(6-position))
  tapply(which(Indices),Factors,function(i)sum(Counts[i],na.rm=TRUE))
}

bootstrap<-function(x=c(33,12,21),M=10000,alpha=0.05){
N<-sum(x)
if(N){
p<-x/N
k<-length(x)-1
tmp<-rmultinom(M, size = N, prob=p)
tmp_p<-apply(tmp,2,function(y)y/N)
(apply(tmp_p,1,function(y)quantile(y,c(alpha/2/k,1-alpha/2/k))))
}
else return(matrix(0,2,length(x)))
}




bootstrap2<-function(x=c(33,12,21),n=10,M=10000,alpha=0.05){

N<-sum(x)
k<-length(x)
y<-rep(1:k,x)
tmp<-sapply(1:M,function(i)sample(y,n))
if(n>1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[,j]==i)))/n
if(n==1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[j]==i)))/n
(apply(tmp_p,1,function(z)quantile(z,c(alpha/2/(k-1),1-alpha/2/(k-1)))))
}



p_value<-function(x=c(33,12,21),M=100000,x_obs=c(2,5,3)){
n=sum(x_obs)
N<-sum(x)
k<-length(x)
y<-rep(1:k,x)
tmp<-sapply(1:M,function(i)sample(y,n))
if(n>1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[,j]==i)))
if(n==1)tmp_p<-sapply(1:M,function(j)sapply(1:k,function(i)sum(tmp[j]==i)))
tmp<-rbind(sapply(1:3,function(i)sum(tmp_p[i,]>=x_obs[i])/M),
sapply(1:3,function(i)sum(tmp_p[i,]<=x_obs[i])/M))
sapply(1:3,function(i){if(tmp[1,i]>=tmp[2,i])return(-tmp[2,i])else return(tmp[1,i])})
}

#"D:\\Sequences\\IMGT Germlines\\Human_SNPless_IGHJ.FASTA"
# Remove SNPs from IMGT germline segment alleles
generateUnambiguousRepertoire <- function(repertoireInFile,repertoireOutFile){
  repertoireIn <- read.fasta(repertoireInFile, seqtype="DNA",as.string=T,set.attributes=F,forceDNAtolower=F)
  alleleNames <- sapply(names(repertoireIn),function(x)strsplit(x,"|",fixed=TRUE)[[1]][2])
  SNPs <- tapply(repertoireIn,sapply(alleleNames,function(x)strsplit(x,"*",fixed=TRUE)[[1]][1]),function(x){
    Indices<-NULL
    for(i in 1:length(x)){
      firstSeq = s2c(x[[1]])
      iSeq = s2c(x[[i]])
      Indices<-c(Indices,which(firstSeq[1:320]!=iSeq[1:320] & firstSeq[1:320]!="." & iSeq[1:320]!="."  ))
    }
    return(sort(unique(Indices)))
  })
 repertoireOut <- repertoireIn
 repertoireOut <- lapply(names(repertoireOut), function(repertoireName){
                                        alleleName <- strsplit(repertoireName,"|",fixed=TRUE)[[1]][2]
                                        geneSegmentName <- strsplit(alleleName,"*",fixed=TRUE)[[1]][1]
                                        alleleSeq <- s2c(repertoireOut[[repertoireName]])
                                        alleleSeq[as.numeric(unlist(SNPs[geneSegmentName]))] <- "N"
                                        alleleSeq <- c2s(alleleSeq)
                                        repertoireOut[[repertoireName]] <- alleleSeq
                                      })
  names(repertoireOut) <- names(repertoireIn)
  write.fasta(repertoireOut,names(repertoireOut),file.out=repertoireOutFile)                                               
                                      
}






############
groupBayes2 = function(indexes, param_resultMat){
  
  BayesGDist_Focused_CDR = calculate_bayesG( x=param_resultMat[indexes,1], N=apply(param_resultMat[indexes,c(1,2,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[1]/(x[1]+x[2]+x[4])}))
  BayesGDist_Focused_FWR = calculate_bayesG( x=param_resultMat[indexes,3], N=apply(param_resultMat[indexes,c(3,2,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[3]/(x[3]+x[2]+x[4])}))
  #BayesGDist_Local_CDR = calculate_bayesG( x=param_resultMat[indexes,1], N=apply(param_resultMat[indexes,c(1,2)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[1]/(x[1]+x[2])}))
  #BayesGDist_Local_FWR = calculate_bayesG( x=param_resultMat[indexes,3], N=apply(param_resultMat[indexes,c(3,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[3]/(x[3]+x[4])}))
  #BayesGDist_Global_CDR = calculate_bayesG( x=param_resultMat[indexes,1], N=apply(param_resultMat[indexes,c(1,2,3,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[1]/(x[1]+x[2]+x[3]+x[4])}))
  #BayesGDist_Global_FWR = calculate_bayesG( x=param_resultMat[indexes,3], N=apply(param_resultMat[indexes,c(1,2,3,4)],1,sum,na.rm=T), p=apply(param_resultMat[indexes,5:8],1,function(x){x[3]/(x[1]+x[2]+x[3]+x[4])}))
  return ( list("BayesGDist_Focused_CDR"=BayesGDist_Focused_CDR,
                "BayesGDist_Focused_FWR"=BayesGDist_Focused_FWR) )
                #"BayesGDist_Local_CDR"=BayesGDist_Local_CDR,
                #"BayesGDist_Local_FWR" = BayesGDist_Local_FWR))
#                "BayesGDist_Global_CDR" = BayesGDist_Global_CDR,
#                "BayesGDist_Global_FWR" = BayesGDist_Global_FWR) )


}


calculate_bayesG <- function( x=array(), N=array(), p=array(), max_sigma=20, length_sigma=4001){
  G <- max(length(x),length(N),length(p))
  x=array(x,dim=G)
  N=array(N,dim=G)
  p=array(p,dim=G)

  indexOfZero = N>0 & p>0
  N = N[indexOfZero]
  x = x[indexOfZero]
  p = p[indexOfZero]  
  G <- length(x)
  
  if(G){
    
    cons<-array( dim=c(length_sigma,G) )
    if(G==1) {
    return(calculate_bayes(x=x[G],N=N[G],p=p[G],max_sigma=max_sigma,length_sigma=length_sigma))
    }
    else {
      for(g in 1:G) cons[,g] <- calculate_bayes(x=x[g],N=N[g],p=p[g],max_sigma=max_sigma,length_sigma=length_sigma)
      listMatG <- convolutionPowersOfTwoByTwos(cons,length_sigma=length_sigma)
      y<-calculate_bayesGHelper(listMatG,length_sigma=length_sigma)
      return( y/sum(y)/(2*max_sigma/(length_sigma-1)) )
    }
  }else{
    return(NA)
  }
}


calculate_bayesGHelper <- function( listMatG,length_sigma=4001 ){
  matG <- listMatG[[1]]  
  groups <- listMatG[[2]]
  i = 1  
  resConv <- matG[,i]
  denom <- 2^groups[i]
  if(length(groups)>1){
    while( i<length(groups) ){
      i = i + 1
      resConv <- weighted_conv(resConv, matG[,i], w= {{2^groups[i]}/denom} ,length_sigma=length_sigma)
      #cat({{2^groups[i]}/denom},"\n")
      denom <- denom + 2^groups[i]
    }
  }
  return(resConv)  
}

weighted_conv<-function(x,y,w=1,m=100,length_sigma=4001){
lx<-length(x)
ly<-length(y)
if({lx<m}| {{lx*w}<m}| {{ly}<m}| {{ly*w}<m}){
if(w<1){
y1<-approx(1:ly,y,seq(1,ly,length.out=m))$y
x1<-approx(1:lx,x,seq(1,lx,length.out=m/w))$y
lx<-length(x1)
ly<-length(y1)
}
else {
y1<-approx(1:ly,y,seq(1,ly,length.out=m*w))$y
x1<-approx(1:lx,x,seq(1,lx,length.out=m))$y
lx<-length(x1)
ly<-length(y1)
}
}
else{
x1<-x
y1<-approx(1:ly,y,seq(1,ly,length.out=floor(lx*w)))$y
ly<-length(y1)
}
tmp<-approx(x=1:(lx+ly-1),y=convolve(x1,rev(y1),type="open"),xout=seq(1,lx+ly-1,length.out=length_sigma))$y
tmp[tmp<=0] = 0 
return(tmp/sum(tmp))
}

########################




mutabilityMatrixONE<-rep(0,4)
names(mutabilityMatrixONE)<-NUCLEOTIDES

  # triMutability Background Count
  buildMutabilityModelONE <- function( inputMatrixIndex, model=0 , multipleMutation=0, seqWithStops=0, stopMutations=0){
    
    #rowOrigMatInput = matInput[inputMatrixIndex,]    
    seqGL =  gsub("-", "", matInput[inputMatrixIndex,2])
    seqInput = gsub("-", "", matInput[inputMatrixIndex,1])    
    matInput[inputMatrixIndex,] <<- c(seqInput,seqGL)
    seqLength = nchar(seqGL)      
    mutationCount <- analyzeMutations(inputMatrixIndex, model, multipleMutation, seqWithStops)
    BackgroundMatrix = mutabilityMatrixONE
    MutationMatrix = mutabilityMatrixONE    
    MutationCountMatrix = mutabilityMatrixONE    
    if(!is.na(mutationCount)){
      if((stopMutations==0 & model==0) | (stopMutations==1 & (sum(mutationCount=="Stop")<length(mutationCount))) | (model==1 & (sum(mutationCount=="S")>0)) ){ 
                  
#         ONEmerStartPos = 1:(seqLength)
#         ONEmerLength <- length(ONEmerStartPos)
        ONEmerGL <- s2c(seqGL)
        ONEmerSeq <- s2c(seqInput)
    
        #Background
        for(ONEmerIndex in 1:seqLength){
          ONEmer = ONEmerGL[ONEmerIndex]
          if(ONEmer!="N"){
            ONEmerCodonPos = getCodonPos(ONEmerIndex)
            ONEmerReadingFrameCodon = c2s(ONEmerGL[ONEmerCodonPos]) 
            ONEmerReadingFrameCodonInputSeq = c2s(ONEmerSeq[ONEmerCodonPos] )         
            
            # All mutations model
            #if(!any(grep("N",ONEmerReadingFrameCodon))){
              if(model==0){
                if(stopMutations==0){
                  if(!any(grep("N",ONEmerReadingFrameCodonInputSeq)))
                    BackgroundMatrix[ONEmer] <- (BackgroundMatrix[ONEmer] + 1)              
                }else{
                  if( !any(grep("N",ONEmerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(ONEmerReadingFrameCodonInputSeq)!="*"){
                    positionWithinCodon = which(ONEmerCodonPos==ONEmerIndex)#positionsWithinCodon[(ONEmerCodonPos[1]%%3)+1]
                    BackgroundMatrix[ONEmer] <- (BackgroundMatrix[ONEmer] + probNonStopMutations[ONEmerReadingFrameCodon,positionWithinCodon])
                  }
                }
              }else{ # Only silent mutations
                if( !any(grep("N",ONEmerReadingFrameCodonInputSeq)) & translateCodonToAminoAcid(ONEmerReadingFrameCodonInputSeq)!="*" & translateCodonToAminoAcid(ONEmerReadingFrameCodonInputSeq)==translateCodonToAminoAcid(ONEmerReadingFrameCodon) ){
                  positionWithinCodon = which(ONEmerCodonPos==ONEmerIndex)
                  BackgroundMatrix[ONEmer] <- (BackgroundMatrix[ONEmer] + probSMutations[ONEmerReadingFrameCodon,positionWithinCodon])
                }
              }
            }
          }
        }
        
        #Mutations
        if(stopMutations==1) mutationCount = mutationCount[mutationCount!="Stop"]
        if(model==1) mutationCount = mutationCount[mutationCount=="S"]  
        mutationPositions = as.numeric(names(mutationCount))
        mutationCount = mutationCount[mutationPositions>2 & mutationPositions<(seqLength-1)]
        mutationPositions =  mutationPositions[mutationPositions>2 & mutationPositions<(seqLength-1)]
        countMutations = 0 
        for(mutationPosition in mutationPositions){
          ONEmerIndex = mutationPosition
          ONEmer = ONEmerSeq[ONEmerIndex]
          GLONEmer = ONEmerGL[ONEmerIndex]
          ONEmerCodonPos = getCodonPos(ONEmerIndex)
          ONEmerReadingFrameCodon = c2s(ONEmerSeq[ONEmerCodonPos])  
          ONEmerReadingFrameCodonGL =c2s(ONEmerGL[ONEmerCodonPos])  
          if(!any(grep("N",ONEmer)) & !any(grep("N",GLONEmer))){
            if(model==0){
                countMutations = countMutations + 1              
                MutationMatrix[GLONEmer] <- (MutationMatrix[GLONEmer] + 1)
                MutationCountMatrix[GLONEmer] <- (MutationCountMatrix[GLONEmer] + 1)             
            }else{
              if( translateCodonToAminoAcid(ONEmerReadingFrameCodonGL)!="*" ){
                  countMutations = countMutations + 1
                  positionWithinCodon = which(ONEmerCodonPos==ONEmerIndex)
                  glNuc =  substr(ONEmerReadingFrameCodonGL,positionWithinCodon,positionWithinCodon)
                  inputNuc =  substr(ONEmerReadingFrameCodon,positionWithinCodon,positionWithinCodon)
                  MutationMatrix[GLONEmer] <- (MutationMatrix[GLONEmer] + substitution[glNuc,inputNuc])
                  MutationCountMatrix[GLONEmer] <- (MutationCountMatrix[GLONEmer] + 1)                                    
              }                
            }                  
          }              
        }
        
        seqMutability = MutationMatrix/BackgroundMatrix
        seqMutability = seqMutability/sum(seqMutability,na.rm=TRUE)
        #cat(inputMatrixIndex,"\t",countMutations,"\n")
        return(list("seqMutability"  = seqMutability,"numbMutations" = countMutations,"seqMutabilityCount" = MutationCountMatrix, "BackgroundMatrix"=BackgroundMatrix))      
#         tmp<-list("seqMutability"  = seqMutability,"numbMutations" = countMutations,"seqMutabilityCount" = MutationCountMatrix)
      }        
    }
  
################
# $Id: trim.R 989 2006-10-29 15:28:26Z ggorjan $

trim <- function(s, recode.factor=TRUE, ...)
  UseMethod("trim", s)

trim.default <- function(s, recode.factor=TRUE, ...)
  s

trim.character <- function(s, recode.factor=TRUE, ...)
{
  s <- sub(pattern="^ +", replacement="", x=s)
  s <- sub(pattern=" +$", replacement="", x=s)
  s
}

trim.factor <- function(s, recode.factor=TRUE, ...)
{
  levels(s) <- trim(levels(s))
  if(recode.factor) {
    dots <- list(x=s, ...)
    if(is.null(dots$sort)) dots$sort <- sort
    s <- do.call(what=reorder.factor, args=dots)
  }
  s
}

trim.list <- function(s, recode.factor=TRUE, ...)
  lapply(s, trim, recode.factor=recode.factor, ...)

trim.data.frame <- function(s, recode.factor=TRUE, ...)
{
  s[] <- trim.list(s, recode.factor=recode.factor, ...)
  s
}
#######################################
