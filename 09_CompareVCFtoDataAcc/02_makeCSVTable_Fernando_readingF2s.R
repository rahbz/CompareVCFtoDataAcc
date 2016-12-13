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
#kinaccessions <- h5read(kinfile, name="accessions")
#kinship <- h5read(kinfile, name="kinship")

list <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
GivParents <- character()
FouParents1 <- character()
FouParents2 <- character()
FouParents3 <- character()
FouParents4 <- character()
FouParents5 <- character()
FouParents6 <- character()
FouParents7 <- character()
ScoreP1 <- numeric()
ScoreP2 <- numeric()
ScoreP3 <- numeric()
ScoreP4 <- numeric()
ScoreP5 <- numeric()
ScoreP6 <- numeric()
ScoreP7 <- numeric()

SNPsmat1 <- numeric()
SNPsmat2 <- numeric()
SNPsmat3 <- numeric()
SNPsmat4 <- numeric()
SNPsmat5 <- numeric()
SNPsmat6 <- numeric()
SNPsmat7 <- numeric()

FracSNPsobs <- numeric()
FracSNPsexp <- numeric()
KinEstimate <- character()
SNPsCalled <- numeric()
HetFreq <- character()

for (echscore in list){
  name <- sub("[.]ScoreAcc[.]txt", "", echscore)
  ScoreAcc <- read.table(echscore, header = F)
  filename <- unlist(strsplit(name, "[.]"))[1]
  givparents <- sub("F2_p1_#", "", filename)
  expP1 <- unlist(strsplit(givparents, "x"))[1]
  expP2 <- unlist(strsplit(givparents, "x"))[2]

  ranks <- order(-ScoreAcc$V4)

  tophits <- ScoreAcc$V1[ranks]
  numsnps <- ScoreAcc$V5[ranks][1]
  scores <- ScoreAcc$V4[ranks]
  matsnps <- ScoreAcc$V2[ranks]
  
#  tophitsacc <- which(kinaccessions %in% tophits)
#  KinEstimate <- c(KinEstimate, kinship[tophitsacc[1], tophitsacc[2]])
  SNPsmat1 <- c(SNPsmat1, matsnps[1])
  SNPsmat2 <- c(SNPsmat2, matsnps[2])
  SNPsmat3 <- c(SNPsmat3, matsnps[3])
  SNPsmat4 <- c(SNPsmat4, matsnps[4])
  SNPsmat5 <- c(SNPsmat5, matsnps[5])
  SNPsmat6 <- c(SNPsmat6, matsnps[6])
  SNPsmat7 <- c(SNPsmat7, matsnps[7])
  ScoreP1 <- c(ScoreP1, scores[1])
  ScoreP2 <- c(ScoreP2, scores[2])
  ScoreP3 <- c(ScoreP3, scores[3])
  ScoreP4 <- c(ScoreP4, scores[4])
  ScoreP5 <- c(ScoreP5, scores[5])
  ScoreP6 <- c(ScoreP6, scores[6])
  ScoreP7 <- c(ScoreP7, scores[7])
  FouParents1 <- c(FouParents1, tophits[1])
  FouParents2 <- c(FouParents2, tophits[2])
  FouParents3 <- c(FouParents3, tophits[3])
  FouParents4 <- c(FouParents4, tophits[4])
  FouParents5 <- c(FouParents5, tophits[5])
  FouParents6 <- c(FouParents6, tophits[6])
  FouParents7 <- c(FouParents7, tophits[7])

  SNPsCalled <- c(SNPsCalled, numsnps)

  FracSNPsobs <- c(FracSNPsobs, scores[2]/scores[1])
  if(length(which(tophits == expP1)) >= 1 & length(which(tophits == expP2)) >= 1){
    fracExp <- scores[which(tophits == expP1)]/scores[which(tophits == expP2)]
    if(fracExp > 1){
      fracExp <- scores[which(tophits == expP2)]/scores[which(tophits == expP1)]
    }
  } else {
    fracExp = "NA"
  }
  FracSNPsexp <- c(FracSNPsexp, fracExp)

   posfile <- paste(name, ".pos.txt",sep = "")
   targetSNPs <- read.table(posfile, header = F)
   numhet <- 0
   totnum <- 0
   for (chr in the1001chr.ids){
     perchrTarget <- targetSNPs[which(targetSNPs$V1 == chr),]
     chrpositions <- the1001positions[the1001chromosomes[,as.numeric(chr)][1]:the1001chromosomes[,as.numeric(chr)][2]]
     matchedSNPs <- which(perchrTarget$V2 %in% chrpositions)
     numhet <- numhet + length(which(perchrTarget$V3[matchedSNPs] == "0/1"))
     totnum <- totnum + length(which(perchrTarget$V3[matchedSNPs] == "1/1"))
   }
   perhet <- 100*(numhet/(totnum+numhet))
   HetFreq <- c(HetFreq, perhet)
  
  Names <- c(Names, name)
  GivParents <- c(GivParents, givparents)
  print("Done one file")
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
print(length(FracSNPsobs))
print(length(FracSNPsexp))
print(length(KinEstimate))
print(length(SNPsmat1))
print(length(SNPsmat2))
print(length(SNPsCalled))
print(length(HetFreq))


#DF <- data.frame(FILENAME = Names, ExpectedParents = GivParents, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, Score = ScoreP1, ObservedFracEstimate = as.numeric(FracSNPsobs), ExpectedFracEstimate = FracSNPsexp, SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), KinEstimate = as.numeric(KinEstimate), HetFreq = HetFreq, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6, SeventhHit = FouParents7, SNPsCalled = SNPsCalled, ScoreH3 = ScoreP3, SNPsmatH3 = SNPsmat3, ScoreH4 = ScoreP4, SNPsmatH4 = SNPsmat4, ScoreH5 = ScoreP5, SNPsmatH5 = SNPsmat5, ScoreH6 = ScoreP6, SNPsmatH6 = SNPsmat6, ScoreH7 = ScoreP7, SNPsmatH7 = SNPsmat7)

DF <- data.frame(FILENAME = Names, ExpectedParents = GivParents, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, Score = ScoreP1, ObservedFracEstimate = as.numeric(FracSNPsobs), ExpectedFracEstimate = FracSNPsexp, SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), HetFreq = HetFreq, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6, SeventhHit = FouParents7, SNPsCalled = SNPsCalled, ScoreH3 = ScoreP3, SNPsmatH3 = SNPsmat3, ScoreH4 = ScoreP4, SNPsmatH4 = SNPsmat4, ScoreH5 = ScoreP5, SNPsmatH5 = SNPsmat5, ScoreH6 = ScoreP6, SNPsmatH6 = SNPsmat6, ScoreH7 = ScoreP7, SNPsmatH7 = SNPsmat7)

#DF <- data.frame(FILENAME = Names, ExpectedParents = GivParents, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, Score = ScoreP1, ObservedFracEstimate = as.numeric(FracSNPsobs), ExpectedFracEstimate = FracSNPsexp, SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6, SeventhHit = FouParents7, SNPsCalled = SNPsCalled, ScoreH3 = ScoreP3, SNPsmatH3 = SNPsmat3, ScoreH4 = ScoreP4, SNPsmatH4 = SNPsmat4, ScoreH5 = ScoreP5, SNPsmatH5 = SNPsmat5, ScoreH6 = ScoreP6, SNPsmatH6 = SNPsmat6, ScoreH7 = ScoreP7, SNPsmatH7 = SNPsmat7)


write.csv(DF, file = outFile)

