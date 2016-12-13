#Refining of the Scores in close accessions

args <- commandArgs(TRUE)
workingDir <- args[1]
setwd(workingDir)
inname <- args[2]
# Read if there is ambiguity
ScoreAcc <- read.table(paste(inname, ".ScoreAcc.txt", sep = ""), sep = "\t", header = F)
TopHitsNumber <- length(which(ScoreAcc$V4 > max(ScoreAcc$V4)-sd(ScoreAcc$V4)));
reqacc_file <- paste(inname, ".ReqAcc.txt", sep = "")
if(file.exists(reqacc_file)){
  file.remove(reqacc_file)
}
if (TopHitsNumber != 1){
  thres <- max(ScoreAcc$V4) - 2*sd(ScoreAcc$V4)
  AccessionsToCheck <- ScoreAcc$V1[which(ScoreAcc$V4 > thres)]
  write(AccessionsToCheck, file = reqacc_file, ncolumns = length(AccessionsToCheck),append = F,sep = ",")
#  posfile <- paste(inname, ".pos.txt", sep = "")
#  outfile <- paste(inname,".refScoreAcc.txt", sep = "")
#  system(paste("python ~/MyScripts/03_CompareVCFtoDataAcc/02_02_CalcHammingDist_perAcc.py -r", reqacc_file, "-p", posfile, "-o", outfile, sep = " "))
}
