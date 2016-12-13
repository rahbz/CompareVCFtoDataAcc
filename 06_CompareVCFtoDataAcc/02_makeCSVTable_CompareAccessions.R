args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1) Path to HDF5 file of whole imputed dataset
pathHDF5file <- args[1]
#2) Run the script from the same working directory as the ScoreAcc files
#workingDir <- "./"
#setwd(workingDir)
#3) File ID, starting name of the file.
fileid <- args[2]
#4) Path to the accession file, where the accessions supposed to be in a tab-delimited file
pathacc <- args[3]
#5) Path to merged ecotype IDs
pathechotypeIDs <- args[4]
#6) Output CSV file
outFile <- args[5]



library("rhdf5")
the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
the1001chromosomes <- attributes(the1001positions)$chr_regions;
the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
the1001chr.ids <- attributes(the1001positions)$chrs;

mergedAcc <- read.table(pathechotypeIDs, header = F)
# Create a file NumberSNPs_samples.txt
# using a bash script
#system("wc -l *.pos.txt|head -n -1 |tr -s ' '|sed 's/ /\t/g'|cut -f2,3|sed 's/.pos.txt$//' >NumberSNPs_samples.txt")
numSNPs <- read.csv("NumberSNPs_samples.txt", header = F, sep = "\t")
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

# ChoiceofAcc == -2  ====> It is not at all present in the database (whole imputed)
# ChoiceofAcc == -1  ====> merged

#-------
allScoreFiles <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
sampleRows <- character()
sampleCol <- character()
accAssigned <- character()
TopHitAcc <- character()
NextHitAcc <- character()
TopHitScore <- character()
TopHitAccSNPs <- character()
TopHitMatchedSNPs <- character()
NextHitMatchedSNPs <- character()
ScoreSDs <- character()
ScoreNOs <- character()
SNPscalled <- character()
ChoiceAcc <- character()
AFfreq <- character()
for (file in allScoreFiles){
  ScoreAcc <- read.table(file, header = F)
  name <- sub(".ScoreAcc.txt","",file)
  #Converting the name of file into rows and columns
  tempnum <- as.numeric(strsplit(sub(fileid, "", name), "_")[[1]][1])-700
  tempidnum <- as.numeric(sub("50","",strsplit(sub(fileid, "", name), "_")[[1]][2]))
  tempid <- chartr("12345678","ABCDEFGH",sub("50","",strsplit(sub(fileid, "", name), "_")[[1]][2]))
  sampleRows <- c(sampleRows, tempid)
  sampleCol <- c(sampleCol, tempnum)
  maxscore <- which(ScoreAcc$V4 == max(ScoreAcc$V4))
  topacc <- ScoreAcc$V1[maxscore]
  snps <- ScoreAcc$V5[maxscore]
#  snps <- as.numeric(numSNPs$V1[which(numSNPs$V2 == name)])
  nextacc <- ScoreAcc$V1[order(-ScoreAcc$V4)][2]
  topaccsnps <- as.numeric(ScoreAcc$V3[maxscore])
  topmatchsnps <- as.numeric(ScoreAcc$V2[maxscore])
  nexthitmatchsnps <- as.numeric(ScoreAcc$V2[order(-ScoreAcc$V4)][2])
  Names <- c(Names, name)
  # accass <- AssignedAcc$sample.description[which(AssignedAcc$row == tempid & AssignedAcc$column == tempnum)]
  accass <- correctedAssAcc[tempidnum, tempnum]
  choicenum <- which(ScoreAcc$V1[order(-ScoreAcc$V4)] == accass)
  if (length(choicenum)){
    if(is.na(accass)){
      accAssigned <- c(accAssigned, "NA")
    }
  } else {
    pmerged <- as.character(mergedAcc$V1[grep(accass, mergedAcc$V1)]) # Might be merged
    if (length(pmerged)){
      if(length(grep(topacc, pmerged))){
        choicenum = "1"
      } else {
        emerged <- unlist(strsplit(pmerged, ","))[unlist(strsplit(pmerged, ",")) != accass]
        preacc <- emerged[which(emerged %in% ScoreAcc$V1)]
        choicenum = which(ScoreAcc$V1[order(-ScoreAcc$V4)] == preacc)
        if(length(choicenum) == 0){
          choicenum <- which(ScoreAcc$V1[order(-ScoreAcc$V4)] == preacc)
        }
      }
      if (!is.na(accass)){
        merlist <- unlist(strsplit(pmerged, ","))[unlist(strsplit(pmerged, ",")) != accass]
        accass_corrected <- merlist[which(merlist %in% ScoreAcc$V1)]
        accass <- accass_corrected
      } else {
        accass <- "NA"
      }
    } else {
      choicenum = "-2"
    }
  }
 posfile <- paste(name, ".pos.txt",sep = "")
 targetSNPs <- read.table(posfile, header = F)
 numhet <- 0
 totnum <- 0
 for (i in the1001chr.ids){
    start <- the1001chromosomes[,as.numeric(i)][1]
    end <- the1001chromosomes[,as.numeric(i)][2]
    pchrtargetSNPs <- targetSNPs[which(targetSNPs$V1 == as.numeric(i)),]
    chrpositions <- the1001positions[start:end]
    matchedINFO <- pchrtargetSNPs$V3[which(pchrtargetSNPs$V2 %in% chrpositions)]
    allINFO <- unlist(strsplit(as.character(matchedINFO), ";"))
    afs <- as.numeric(gsub("AF=", "",allINFO[grep("AF=", allINFO)]))
    numhet <- numhet + length(which(afs == "0.5"))
    totnum <- totnum + length(afs)
  }
  perhet <- 100*(numhet/totnum)
  AFfreq <- c(AFfreq, as.numeric(perhet))


  topscore <- as.numeric(max(ScoreAcc$V4))
  ScoreAccSD = as.numeric((max(ScoreAcc$V4) - mean(ScoreAcc$V4))/sd(ScoreAcc$V4));
  ScoreAccNO = length(which(ScoreAcc$V4 > max(ScoreAcc$V4)-2*sd(ScoreAcc$V4)))
  if(file.exists(paste(name, ".refScoreAcc.txt", sep = ""))){
     refScore <- read.table(paste(name,".refScoreAcc.txt", sep = ""), header = F)
     topacc <- refScore$V1[which(refScore$V4 == max(refScore$V4))]
     nextacc <- refScore$V1[order(-refScore$V4)][2]
     if(length(which(refScore$V1[order(-refScore$V4)] == accass))){
       choicenum <- which(refScore$V1[order(-refScore$V4)] == accass)
     } else {
       choicenum <- which(ScoreAcc$V1[order(-ScoreAcc$V4)] == accass)
     }
  }
  if(length(choicenum) == 0){
    choicenum = "-2"
  }
  if (length(accass) == 0){
    accass <- "NA"
  }
  accAssigned <- c(accAssigned, accass)
  TopHitAcc <- c(TopHitAcc, topacc)
  NextHitAcc <- c(NextHitAcc, nextacc)
  TopHitScore <- c(TopHitScore, topscore)
  TopHitAccSNPs <- c(TopHitAccSNPs, topaccsnps)
  TopHitMatchedSNPs <- c(TopHitMatchedSNPs, topmatchsnps)
  NextHitMatchedSNPs <- c(NextHitMatchedSNPs, nexthitmatchsnps)
  ScoreSDs <- c(ScoreSDs, ScoreAccSD)
  ScoreNOs <- c(ScoreNOs, ScoreAccNO)
  SNPscalled <- c(SNPscalled, snps)
  ChoiceAcc <- c(ChoiceAcc, choicenum)
}

#DF <- data.frame(FILENAME = Names, ROW = sampleRows, COL = sampleCol, AssignedAcccession = accAssigned,TopHitAccession = TopHitAcc, NextHit = NextHitAcc, Score = as.numeric(TopHitScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsCalled = as.numeric(SNPscalled), SNPsinAcc  = as.numeric(TopHitAccSNPs), tStat  = as.numeric(ScoreSDs),TopHitsNumber = ScoreNOs, ChoiceofAcc = ChoiceAcc)

DF <- data.frame(FILENAME = Names, ROW = sampleRows, COL = sampleCol, AssignedAcccession = accAssigned,TopHitAccession = TopHitAcc, NextHit = NextHitAcc, Score = as.numeric(TopHitScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), NextHitMatchedSNPs = as.numeric(NextHitMatchedSNPs), SNPsCalled = as.numeric(SNPscalled), SNPsinAcc  = as.numeric(TopHitAccSNPs), tStat  = as.numeric(ScoreSDs), HetFreq = as.numeric(AFfreq), TopHitsNumber = ScoreNOs, ChoiceofAcc = ChoiceAcc)

#DF <- data.frame(FILENAME = Names, ROW = sampleRows, COL = sampleCol, AssignedAcccession = accAssigned,TopHitAccession = TopHitAcc, Score = as.numeric(TopHitScore), MatchedSNPs = as.numeric(TopHitMatchedSNPs), SNPsCalled = as.numeric(SNPscalled), SNPsinAcc  = as.numeric(TopHitAccSNPs), tStat  = as.numeric(ScoreSDs), HetFreq = as.numeric(AFfreq),TopHitsNumber = ScoreNOs, ChoiceofAcc = ChoiceAcc)
#fix(DF)
#----

write.csv(DF, file = outFile)



