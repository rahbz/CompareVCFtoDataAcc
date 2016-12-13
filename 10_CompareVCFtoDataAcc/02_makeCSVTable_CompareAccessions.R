args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1) Path to HDF5 file of whole imputed dataset
#pathHDF5file <- args[1]
#2) Run the script from the same working directory as the ScoreAcc files
#workingDir <- "./"
#setwd(workingDir)
#4) Path to the accession file, where the accessions supposed to be in a tab-delimited file
pathacc <- args[3]
#5) Path to merged ecotype IDs
pathechotypeIDs <- args[1]
#6) Output CSV file
outFile <- args[2]
folid <- args[4]



#library("rhdf5")
#the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
#the1001chromosomes <- attributes(the1001positions)$chr_regions;
#the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
#the1001chr.ids <- attributes(the1001positions)$chrs;
#the1001accessions <- h5read(pathHDF5file, name="accessions")


mergedAcc <- read.table(pathechotypeIDs, header = F)
AssignedAcc <- read.csv(pathacc, header = T, sep = "\t")

# Correction for changing the rows to coloumns
rows <- c("A","B","C","D","E","F","G","H")
cols <- c("1","2","3","4","5","6","7","8","9","10","11","12")
correctedAssAcc <- matrix(
  AssignedAcc$sample.description,
  nrow = 8,
  ncol = 12,
  dimnames = list(rows, cols))
if(file.exists(outFile)){
  file.remove(outFile)
}

#-------
allScoreFiles <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
sampleRows <- character()
sampleCol <- character()
accAssigned <- character()
TopHitAcc <- character()
NextHitAcc <- character()
TopHitScore <- numeric()
FracScore <- numeric()
TopHitAccSNPs <- numeric()
TopHitMatchedSNPs <- numeric()
NextHitMatchedSNPs <- numeric()
SNPscalled <- numeric()
ChoiceAcc <- numeric()
AFfreq <- numeric()
LikeLihoodTopHit <- numeric()
LLRNextHit <- numeric()
for (file in allScoreFiles){
  ScoreAcc <- read.table(file, header = F)
  name <- sub(".ScoreAcc.txt","",file)
  #Converting the name of file into rows and columns
  tempnum <- as.numeric(strsplit(name, "_")[[1]][1])-700
  #tempidnum <- as.numeric(sub("50","",strsplit(name, "_")[[1]][2]))
  tempid <- chartr("12345678","ABCDEFGH",sub("50","",strsplit(name, "_")[[1]][2]))
  sampleRows <- c(sampleRows, tempid)
  sampleCol <- c(sampleCol, tempnum)
## Changed the sorting order based on the score to the likelihood ratio
  ranks <- order(ScoreAcc$V5)

  topscore <- ScoreAcc$V4[ranks[1]]
  topacc <- ScoreAcc$V1[ranks[1]]
  snps <- ScoreAcc$V6[ranks[1]]
  maxlike <- ScoreAcc$V5[ranks[1]]

  newLike <- ScoreAcc$V5[ranks]/ScoreAcc$V5[ranks[1]]
  nextlike <- newLike[2]
  
#  snps <- as.numeric(numSNPs$V1[which(numSNPs$V2 == name)])
  nextacc <- ScoreAcc$V1[ranks[2]]
  nextscore <- ScoreAcc$V4[ranks[2]]
  frac <- nextscore/topscore

  topaccsnps <- ScoreAcc$V3[ranks[1]]
  topmatchsnps <- ScoreAcc$V2[ranks[1]]
  nexthitmatchsnps <- ScoreAcc$V2[ranks[2]]
  Names <- c(Names, name)
  accass <- AssignedAcc$sample.description[which(AssignedAcc$row == tempid & AssignedAcc$column == tempnum)]
  choicenum <- which(ScoreAcc$V1[ranks] == accass)
  
  if (length(choicenum) == 0){
    pmerged <- as.character(mergedAcc$V1[grep(accass, mergedAcc$V1)]) # Might be merged
    if (length(pmerged)){
      if(length(grep(topacc, pmerged))){
        accass <- topacc
        choicenum = "1"
      } else {
        emerged <- unlist(strsplit(pmerged, ","))
        preacc <- emerged[which(emerged %in% ScoreAcc$V1)]
        if (length(preacc)){
          choicenum = which(ScoreAcc$V1[ranks] == preacc)
          accass <- preacc
        } else {
          choicenum = "-2"
        }
      }
    } else {
      choicenum = "-2"
    }
  }

#   if(file.exists(paste(name, ".refScore.txt", sep = ""))){
#      if(length(readLines(paste(name, ".refScore.txt", sep = ""))) > 0){
#        refScore <- read.table(paste(name,".refScore.txt", sep = ""), header = F)
#        refRanks <- order(refScore$V5)
#        topacc <- refScore$V1[refRanks[1]]
#        nextacc <- refScore$V1[refRanks[2]]
#        if(length(which(refScore$V1[refRanks] == accass))){
#          choicenum <- which(refScore$V1[refRanks] == accass)
#        }
#     }
#   }

  FracScore <- c(FracScore, frac)
  accAssigned <- c(accAssigned, accass)
  TopHitAcc <- c(TopHitAcc, topacc)
  NextHitAcc <- c(NextHitAcc, nextacc)
  TopHitScore <- c(TopHitScore, topscore)
  TopHitAccSNPs <- c(TopHitAccSNPs, topaccsnps)
  TopHitMatchedSNPs <- c(TopHitMatchedSNPs, topmatchsnps)
  NextHitMatchedSNPs <- c(NextHitMatchedSNPs, nexthitmatchsnps)
  SNPscalled <- c(SNPscalled, snps)
  ChoiceAcc <- c(ChoiceAcc, choicenum)
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

DF <- data.frame(FOL = fol, FILENAME = Names, ROW = sampleRows, COL = sampleCol, AssignedAcccession = accAssigned,TopHitAccession = TopHitAcc, NextHit = NextHitAcc, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsCalled = as.numeric(SNPscalled), SNPsinfoAcc  = as.numeric(TopHitAccSNPs), LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, ChoiceofAcc = ChoiceAcc)

#DF <- data.frame(FILENAME = Names, ROW = sampleRows, COL = sampleCol, AssignedAcccession = accAssigned,TopHitAccession = TopHitAcc, NextHit = NextHitAcc, Score = as.numeric(TopHitScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsCalled = as.numeric(SNPscalled), SNPsinAcc  = as.numeric(TopHitAccSNPs), tStat  = as.numeric(ScoreSDs), HetFreq = as.numeric(AFfreq), TopHitsNumber = ScoreNOs, ChoiceofAcc = ChoiceAcc)

#DF <- data.frame(FILENAME = Names, ROW = sampleRows, COL = sampleCol, AssignedAcccession = accAssigned,TopHitAccession = TopHitAcc, Score = as.numeric(TopHitScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), SNPsCalled = as.numeric(SNPscalled), SNPsinAcc  = as.numeric(TopHitAccSNPs), tStat  = as.numeric(ScoreSDs), HetFreq = as.numeric(AFfreq),TopHitsNumber = ScoreNOs, ChoiceofAcc = ChoiceAcc)
#fix(DF)
#----

write.csv(DF, file = outFile)



