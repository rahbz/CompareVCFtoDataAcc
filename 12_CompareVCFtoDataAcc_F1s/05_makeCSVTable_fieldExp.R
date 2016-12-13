args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1) Path to HDF5 file of whole imputed dataset
pathHDF5file <- args[1]
#2) Run the script from the same working directory as the ScoreAcc files
#workingDir <- "./"
#setwd(workingDir)
#4) Path to the accession file, where the accessions supposed to be in a tab-delimited file
#5) Path to merged ecotype IDs
#6) Output CSV file
outFile <- args[2]
folid <- args[3]



library("rhdf5")
#the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
#the1001chromosomes <- attributes(the1001positions)$chr_regions;
#the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
#the1001chr.ids <- attributes(the1001positions)$chrs;
the1001accessions <- h5read(pathHDF5file, name="accessions")



# Correction for changing the rows to coloumns
rows <- c("A","B","C","D","E","F","G","H")
cols <- c("1","2","3","4","5","6","7","8","9","10","11","12")
if(file.exists(outFile)){
  file.remove(outFile)
}

#-------
allScoreFiles <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
sampleRows <- character()
sampleCol <- character()
TopHitAcc <- character()
NextHitAcc <- character()
ThirdHit <- character()
TopHitScore <- numeric()
FracScore <- numeric()
TopHitAccSNPs <- numeric()
TopHitMatchedSNPs <- numeric()
NextHitMatchedSNPs <- numeric()
SNPscalled <- numeric()
AFfreq <- numeric()
LikeLihoodTopHit <- numeric()
LLRNextHit <- numeric()
TopHitsNumber <- numeric()
TopHits <- character()
ChoiceAcc <- numeric()
nF1Likeli <- numeric()
PerHetGenotype <- numeric()
for (file in allScoreFiles){
  ScoreAcc <- try(read.table(file, header = F))
  if(inherits(ScoreAcc, 'try-error')){
    next
  }

  name <- sub(".ScoreAcc.txt","",file)
  #Converting the name of file into rows and columns
  nums <- sub(folid, "", name)
  tempnum <- as.numeric(strsplit(nums, "_")[[1]][1])-700
  #tempidnum <- as.numeric(sub("50","",strsplit(name, "_")[[1]][2]))
  tempid <- chartr("12345678","ABCDEFGH",sub("50","",strsplit(nums, "_")[[1]][2]))
  sampleRows <- c(sampleRows, tempid)
  sampleCol <- c(sampleCol, tempnum)

## Changed the sorting order based on the score to the likelihood ratio
#  ranks <- order(-ScoreAcc$V4, na.last = TRUE, partial = ScoreAcc$V6)
  ranks <- order(ScoreAcc$V6, na.last = TRUE, partial = ScoreAcc$V4)

  topscore <- ScoreAcc$V4[ranks[1]]
  topacc <- as.character(ScoreAcc$V1[ranks[1]])
  snps <- ScoreAcc$V7[ranks[1]]
  maxlike <- ScoreAcc$V5[ranks[1]]

  newLike <- ScoreAcc$V5[ranks]/ScoreAcc$V5[ranks[1]]
  nextlike <- newLike[2]
  
#  snps <- as.numeric(numSNPs$V1[which(numSNPs$V2 == name)])
  nextacc <- as.character(ScoreAcc$V1[ranks[2]])
  nextscore <- ScoreAcc$V4[ranks[2]]
  frac <- nextscore/topscore

  topaccsnps <- ScoreAcc$V3[ranks[1]]
  topmatchsnps <- ScoreAcc$V2[ranks[1]]
  nexthitmatchsnps <- ScoreAcc$V2[ranks[2]]
  Names <- c(Names, name)
  
  thirdHit <- as.character(ScoreAcc$V1[ranks[3]])

#  if(file.exists(paste(name, ".refScore.txt", sep = ""))){
#    if(length(readLines(paste(name, ".refScore.txt", sep = ""))) > 0){
#      refScore <- read.table(paste(name,".refScore.txt", sep = ""), header = F)
#      refRanks <- order(refScore$V4)
#      pars <- as.character(refScore$V1[refRanks[1]])
#      parlikeli <- refScore$V4[refRanks[1]]
#      parnextlikeli <- refScore$V4[refRanks[2]] - refScore$V4[refRanks[1]]
#    }
#  }
#  Parents <- c(Parents, pars)
#  ParentLikeli <- c(ParentLikeli, parlikeli)
#  ParentNextLikeli <- c(ParentNextLikeli, parnextlikeli)
  nf1 <- min(ScoreAcc$V6[which(ScoreAcc$V1 %in% the1001accessions)])
  nF1Likeli <- c(nF1Likeli, nf1)

#  topnum <- length(which(ScoreAcc$V6 < 3.841))
  rtopnum <- length(which(ScoreAcc$V6 < 3.841))
  topnum <- length(which(ScoreAcc$V1[which(ScoreAcc$V6 < 3.841)] %in% the1001accessions))
  if (rtopnum != topnum){
    topnum = topnum + 1
  }
  TopHitsNumber <- c(TopHitsNumber, topnum)
  if(topnum > 20){
    TopHits <- c(TopHits, "TooManytoPrint")
  } else if(topnum > 3) {
    alltops <- paste(ScoreAcc$V1[ranks[4:topnum]], collapse = ":")
    TopHits <- c(TopHits, alltops)
  } else {
    TopHits <- c(TopHits, "NA")
  }

#  choiceacc <- which(ScoreAcc$V1[ranks] == name)
#  if(length(choiceacc)){
#    ChoiceAcc <- c(ChoiceAcc, choiceacc[1])
#  } else {
#    ChoiceAcc <- c(ChoiceAcc, -1)
#  }
   

  PerHetGenotype <- c(PerHetGenotype, ScoreAcc$V9[1])
  FracScore <- c(FracScore, frac)
  TopHitAcc <- c(TopHitAcc, topacc)
  NextHitAcc <- c(NextHitAcc, nextacc)
  ThirdHit <- c(ThirdHit, thirdHit)
  TopHitScore <- c(TopHitScore, topscore)
  TopHitAccSNPs <- c(TopHitAccSNPs, topaccsnps)
  TopHitMatchedSNPs <- c(TopHitMatchedSNPs, topmatchsnps)
  NextHitMatchedSNPs <- c(NextHitMatchedSNPs, nexthitmatchsnps)
  SNPscalled <- c(SNPscalled, snps)
  LikeLihoodTopHit <- c(LikeLihoodTopHit, maxlike)
  LLRNextHit <- c(LLRNextHit, nextlike)
}

#print(length(Names))
#print(length(sampleRows))
#print(length(sampleCol))
#print(length(accAssigned))
#print(length(TopHitAcc))
#print(length(NextHitAcc))
#print(length(TopHitScore))
#print(length(TopHitAccSNPs))
#print(length(TopHitMatchedSNPs))
#print(length(NextHitMatchedSNPs))
#print(length(SNPscalled))
#print(length(ChoiceAcc))
#print(length(AFfreq))
#print(length(FracScore))

fol = rep(folid, length(Names))

#DF <- data.frame(FOL = fol, FILENAME = Names, ROW = sampleRows, COL = sampleCol, TopHitAccession = TopHitAcc, NextHit = NextHitAcc, ThirdHit = ThirdHit, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsinfoAcc  = as.numeric(TopHitAccSNPs), SNPsCalled = as.numeric(SNPscalled), LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, TopHitsNumber = TopHitsNumber, TopHits = TopHits)

DF <- data.frame(FOL = fol, FILENAME = Names, TopHitAccession = TopHitAcc, NextHit = NextHitAcc, ThirdHit = ThirdHit, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsinfoAcc  = as.numeric(TopHitAccSNPs), SNPsCalled = as.numeric(SNPscalled), LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, TopHitsNumber = TopHitsNumber, TopHits = TopHits, HetPer = PerHetGenotype, nF1Likeli = nF1Likeli)
## Added choice acc from the top dataframe
#DF <- data.frame(FOL = fol, FILENAME = Names, TopHitAccession = TopHitAcc, NextHit = NextHitAcc, ThirdHit = ThirdHit, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsinfoAcc  = as.numeric(TopHitAccSNPs), SNPsCalled = as.numeric(SNPscalled), LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, TopHitsNumber = TopHitsNumber, TopHits = TopHits, ChoiceAcc = ChoiceAcc)

#----

write.csv(DF, file = outFile)


