args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1) Path to HDF5 file of whole imputed dataset
#pathHDF5file <- args[1]
#kinfile <- args[2]
outFile <- args[1]
#2) Run the script from the same working directory as the ScoreAcc files
#workingDir <- "./"
#setwd(workingDir)
#3) File ID, starting name of the file.
id <- args[2]

#library("rhdf5")
#the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
#the1001chromosomes <- attributes(the1001positions)$chr_regions;
#the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
#the1001chr.ids <- attributes(the1001positions)$chrs;
#kinaccessions <- h5read(kinfile, name="accessions")
#kinship <- h5read(kinfile, name="kinship")

list <- list.files("./", pattern = "[.]matLikeliAcc.txt$")
Names <- character()
plate <- character()
GivParents <- character()
FouParents1 <- character()
FouParents2 <- character()

ScoreP1 <- numeric()
ScoreP2 <- numeric()
NumWindows <- numeric()
CountP1 <- numeric()
CountP2 <- numeric()

#SNPsCalled <- numeric()
#HetFreq <- character()
AmbiWindows <- character()
AmbiWindowCount <- character()

for (echscore in list){
  name <- sub("[.]matLikeliAcc[.]txt", "", echscore)
  likeLiwind <- try(read.table(echscore, header = F, as.is = TRUE))
  if(inherits(likeLiwind, 'try-error')){
    next
  } 
  filename <- unlist(strsplit(name, "[.]"))[1]
  givparents <- sub("F2_p1_#", "", filename)
  expP1 <- unlist(strsplit(givparents, "x"))[1]
  expP2 <- unlist(strsplit(givparents, "x"))[2]
 
  homowind <- which(likeLiwind$V7 == 1) 
  clean <- sort(table(likeLiwind$V1[homowind]))
  clean <- as.data.frame(clean)
  topHits <- rownames(tail(clean, 2))
  
  ScoreP2 <- c(ScoreP2, mean(likeLiwind$V4[homowind[which(likeLiwind$V1[homowind] == topHits[2])]], na.rm = TRUE))
  ScoreP1 <- c(ScoreP1, mean(likeLiwind$V4[homowind[which(likeLiwind$V1[homowind] == topHits[1])]], na.rm = TRUE))
  FouParents1 <- c(FouParents1, as.character(topHits[1]))
  FouParents2 <- c(FouParents2, as.character(topHits[2]))
  CountP2 <- c(CountP2, as.numeric(tail(clean, 2)[2,]))
  CountP1 <- c(CountP1, as.numeric(tail(clean, 2)[1,]))

  ambWind <- which(likeLiwind$V7 < 20)
  nclean <- sort(-table(likeLiwind$V1[ambWind]))
  nclean <- as.data.frame(nclean)
  nclean <- nclean[which(nclean$nclean <= -1),]
  
  AmbiWindows <- c(AmbiWindows, paste(rownames(nclean), collapse = ":"))
  AmbiWindowCount <- c(AmbiWindowCount, paste(as.numeric(-nclean), collapse = ":"))
  
  NumWindows <- c(NumWindows, nrow(likeLiwind))
  plate <- c(plate, id)
#  SNPsCalled <- c(SNPsCalled, numsnps)
#
#
#   posfile <- paste(name, ".pos.txt",sep = "")
#   targetSNPs <- read.table(posfile, header = F)
#   numhet <- 0
#   totnum <- 0
#   for (chr in the1001chr.ids){
#     perchrTarget <- targetSNPs[which(targetSNPs$V1 == chr),]
#     chrpositions <- the1001positions[the1001chromosomes[,as.numeric(chr)][1]:the1001chromosomes[,as.numeric(chr)][2]]
#     matchedSNPs <- which(perchrTarget$V2 %in% chrpositions)
#     numhet <- numhet + length(which(perchrTarget$V3[matchedSNPs] == "0/1"))
#     totnum <- totnum + length(which(perchrTarget$V3[matchedSNPs] == "1/1"))
#   }
#   perhet <- 100*(numhet/(totnum+numhet))
#   HetFreq <- c(HetFreq, perhet)
#  
  Names <- c(Names, name)
  GivParents <- c(GivParents, givparents)
#  print("Done one file")
}


#print(length(Names))
#print(length(GivParents))
#print(length(FouParents1))
#print(length(FouParents2))
#print(length(ScoreP1))
#print(length(ScoreP2))
#print(length(FracSNPsobs))
#print(length(FracSNPsexp))
#print(length(KinEstimate))
#print(length(SNPsmat1))
#print(length(SNPsmat2))
#print(length(SNPsCalled))
#print(length(HetFreq))



DF <- data.frame(PLATE = plate, FILENAME = Names, ExpectedParents = GivParents, NumberWindows = NumWindows, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, CountP1 = CountP1, CountP2 = CountP2, ScoreP1 = ScoreP1, ScoreP2 = ScoreP2, AmbigousWindows = AmbiWindows, AmbigousWindowCount = AmbiWindowCount)

#DF <- data.frame(FILENAME = Names, ExpectedParents = GivParents, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, Score = ScoreP1, ObservedFracEstimate = as.numeric(FracSNPsobs), ExpectedFracEstimate = FracSNPsexp, SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), KinEstimate = as.numeric(KinEstimate), HetFreq = HetFreq, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6, SeventhHit = FouParents7, SNPsCalled = SNPsCalled, ScoreH3 = ScoreP3, SNPsmatH3 = SNPsmat3, ScoreH4 = ScoreP4, SNPsmatH4 = SNPsmat4, ScoreH5 = ScoreP5, SNPsmatH5 = SNPsmat5, ScoreH6 = ScoreP6, SNPsmatH6 = SNPsmat6, ScoreH7 = ScoreP7, SNPsmatH7 = SNPsmat7)

#DF <- data.frame(FILENAME = Names, ExpectedParents = GivParents, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, Score = ScoreP1, ObservedFracEstimate = as.numeric(FracSNPsobs), ExpectedFracEstimate = FracSNPsexp, SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), HetFreq = HetFreq, ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6, SeventhHit = FouParents7, SNPsCalled = SNPsCalled, ScoreH3 = ScoreP3, SNPsmatH3 = SNPsmat3, ScoreH4 = ScoreP4, SNPsmatH4 = SNPsmat4, ScoreH5 = ScoreP5, SNPsmatH5 = SNPsmat5, ScoreH6 = ScoreP6, SNPsmatH6 = SNPsmat6, ScoreH7 = ScoreP7, SNPsmatH7 = SNPsmat7)

#DF <- data.frame(FILENAME = Names, ExpectedParents = GivParents, ObservedParent1 = FouParents1, ObservedParent2 = FouParents2, Score = ScoreP1, ObservedFracEstimate = as.numeric(FracSNPsobs), ExpectedFracEstimate = FracSNPsexp, SNPsmatchedP1 = as.numeric(SNPsmat1), SNPsmatchedP2 = as.numeric(SNPsmat2), ThirdHit = FouParents3, FourthHit = FouParents4, FifthHit = FouParents5, SixthHit = FouParents6, SeventhHit = FouParents7, SNPsCalled = SNPsCalled, ScoreH3 = ScoreP3, SNPsmatH3 = SNPsmat3, ScoreH4 = ScoreP4, SNPsmatH4 = SNPsmat4, ScoreH5 = ScoreP5, SNPsmatH5 = SNPsmat5, ScoreH6 = ScoreP6, SNPsmatH6 = SNPsmat6, ScoreH7 = ScoreP7, SNPsmatH7 = SNPsmat7)


write.csv(DF, file = outFile)

