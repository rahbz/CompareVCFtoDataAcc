# Filter the VCF file with the required filters

# Input :
#1) VCF file
#2) HDF5 file

library("rhdf5")
library(VariantAnnotation)
args <- commandArgs(TRUE)
pathHDF5file <- args[1]
pathVCFfile <- args[2]
qthres <- 100

#the1001accessions <- as.numeric(h5read(pathHDF5file, "accessions"));
#the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
#the1001chromosomes <- attributes(the1001positions)$chr_regions;
#the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
#the1001chr.ids <- attributes(the1001positions)$chrs;

vcf <- readVcf(pathVCFfile, "ARA")
matchedSNPs <- character()

outPOSfile <- args[3]
if(file.exists(outfile)){
  file.remove(outfile)
}

filter <- which(info(vcf)$  qual(vcf) > qthres    )

CHR <- sub("Chr", "", as.character(seqnames(rowData(vcf))[filter]))
POS <- start(rowData(vcf))[filter]
GT <- info(vcf)$GT[filter]

outPOS <- data.frame()

write.table(outPOS, outPOSfile)

