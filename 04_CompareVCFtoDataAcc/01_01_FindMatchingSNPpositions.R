#Get the matched positions from a VCF file
#taking into the hdf5 file

library("rhdf5")
args <- commandArgs(TRUE)
pathHDF5file <- args[1]

accessions <- as.numeric(h5read(pathHDF5file, "accessions"));
positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
chromosomes <- attributes(positions)$chr_regions;
chromosomes[1,] <- chromosomes[1,] + 1;
chr.ids <- attributes(positions)$chrs;

pathposfile <- args[2]
targetSNPs <- read.table(pathposfile, header = F)
outfile <- args[3]
if(file.exists(outfile)){
  file.remove(outfile)
}
for (i in chr.ids){
  start <- chromosomes[,as.numeric(i)][1]
  end <- chromosomes[,as.numeric(i)][2]
  pchrtargetSNPs <- targetSNPs[which(targetSNPs$V1 == as.numeric(i)),]
  chrpositions <- positions[start:end]
  matchedposind <- which(chrpositions %in% pchrtargetSNPs$V2)
  matchedposind <- matchedposind + start - 2
  write(matchedposind, file = outfile, ncolumns = 1000,append = T,sep = ",")
}
