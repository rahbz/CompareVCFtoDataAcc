args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1) Path to HDF5 file of whole imputed dataset
pathHDF5file <- args[1]
kinfile <- args[2]
outFile <- args[3]
#2) Run the script from the same working directory as the ScoreAcc files
#workingDir <- "./"
#setwd(workingDir)
#3) File ID, starting name of the file.


library("rhdf5")
the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
the1001chromosomes <- attributes(the1001positions)$chr_regions;
the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
the1001chr.ids <- attributes(the1001positions)$chrs;
kinaccessions <- h5read(kinfile, name="accessions")
kinship <- h5read(kinfile, name="kinship")

list <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
GivParents <- character()
FouParents1 <- character()
FouParents2 <- character()
FouParents3 <- character()
FouParents4 <- character()
FouParents5 <- character()
FouParents6 <- character()
ScoreP1 <- numeric()
ScoreP2 <- numeric()
FracSNPs1 <- numeric()
FracSNPs2 <- numeric()
KinEstimate <- numeric()
SNPsmat1 <- numeric()
SNPsmat2 <- numeric()
SNPsCalled <- numeric()
NextDist <- numeric()
HetFreq <- character()
pVal <- numeric()
ChoiceNum <- numeric()

for (echscore in list){
  name <- sub("[.]ScoreAcc[.]txt", "", echscore)
  ScoreAcc <- read.table(echscore, header = F)
  filename <- unlist(strsplit(name, "[.]"))[1]
  givparents <- sub("F2_p1_#", "", filename)
  tophits <- ScoreAcc$V1[order(-ScoreAcc$V4)][1:6]
  tophitsacc <- which(kinaccessions %in% tophits)
  numsnps <- ScoreAcc$V5[order(-ScoreAcc$V4)][1]
  #numsnps <- numSNPs$V1[which(numSNPs$V2 == name)]
  matsnps <- ScoreAcc$V2[order(-ScoreAcc$V4)][1:2]
  
  choicenum <- which(ScoreAcc$V1[order(-ScoreAcc$V4)] == name)
  if(length(choicenum) == 0){
    choicenum = -1
  }
  
  ChoiceNum <- c(ChoiceNum, choicenum)
  ScoreP1 <- c(ScoreP1, ScoreAcc$V4[order(-ScoreAcc$V4)][1])
  ScoreP2 <- c(ScoreP2, ScoreAcc$V4[order(-ScoreAcc$V4)][2])
  FracSNPs1 <- c(FracSNPs1, matsnps[2]/matsnps[1])
  FracSNPs2 <- c(FracSNPs2, ScoreAcc$V4[order(-ScoreAcc$V4)][2]/ScoreAcc$V4[order(-ScoreAcc$V4)][1])
  KinEstimate <- c(KinEstimate, kinship[tophitsacc[1], tophitsacc[2]])
  SNPsmat1 <- c(SNPsmat1, matsnps[1])
  SNPsmat2 <- c(SNPsmat2, matsnps[2])
  FouParents1 <- c(FouParents1, tophits[1])
  FouParents2 <- c(FouParents2, tophits[2])
  FouParents3 <- c(FouParents3, tophits[3])
  FouParents4 <- c(FouParents4, tophits[4])
  FouParents5 <- c(FouParents5, tophits[5])
  FouParents6 <- c(FouParents6, tophits[6])

#   posfile <- paste(name, ".pos.txt",sep = "")
#   targetSNPs <- read.table(posfile, header = F)
#   numhet <- 0
#   totnum <- 0
#   for (i in the1001chr.ids){
#     start <- the1001chromosomes[,as.numeric(i)][1]
#     end <- the1001chromosomes[,as.numeric(i)][2]
#     pchrtargetSNPs <- targetSNPs[which(targetSNPs$V1 == as.numeric(i)),]
#     chrpositions <- the1001positions[start:end]
#     matchedINFO <- pchrtargetSNPs$V3[which(pchrtargetSNPs$V2 %in% chrpositions)]
#     allINFO <- unlist(strsplit(as.character(matchedINFO), ";"))
#     afs <- as.numeric(gsub("AF=", "",allINFO[grep("AF=", allINFO)]))
#     #afs <- as.numeric(sub("AF=", "", matchedAFs))
#     numhet <- numhet + length(which(afs == "0.5"))
#     totnum <- totnum + length(afs)
#   }
#   perhet <- 100*(numhet/totnum)
#   HetFreq <- c(HetFreq, perhet)
  
  Names <- c(Names, name)
  SNPsCalled <- c(SNPsCalled, numsnps)
  GivParents <- c(GivParents, givparents)
}


print(length(Names))
print(length(GivParents))
print(length(FouParents1))
print(length(FouParents2))
print(length(FouParents3))
print(length(FouParents4))
print(length(FouParents5))
print(length(FouParents6))
print(length(ScoreP1))
print(length(ScoreP2))
print(length(FracSNPs1))
print(length(FracSNPs2))
print(length(KinEstimate))
print(length(SNPsmat1))
print(length(SNPsmat2))
print(length(SNPsCalled))
print(length(NextDist))
print(length(HetFreq))


#DF <- data.frame(FILENAME = Names, ExpectedParent1 = FouParents1, ExpectedParent2 = FouParents2, Score = ScoreP1, FracEstimate = FracSNPs2, KinEstimate = as.numeric(KinEstimate), SNPsmatchedP1 = SNPsmat1, SNPsmatchedP2 = SNPsmat2, SNPsCalled = SNPsCalled, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6)
DF <- data.frame(FILENAME = Names, TopHit = FouParents1, NextHit = FouParents2, Score = ScoreP1, FracEstimate = FracSNPs2, ChoiceNum = ChoiceNum, SNPsmatchedP1 = SNPsmat1, SNPsmatchedP2 = SNPsmat2, SNPsCalled = SNPsCalled, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6)

#DF <- data.frame(FILENAME = Names, ExpectedParent1 = FouParents1, ExpectedParent2 = FouParents2, ScoreP1 = ScoreP1, ScoreP2 = ScoreP2, FracEstimate = as.numeric(FracSNPs2), KinEstimate = as.numeric(KinEstimate), HetFreq = as.numeric(HetFreq), SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), SNPsCalled = SNPsCalled, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6)

write.csv(DF, file = outFile)

