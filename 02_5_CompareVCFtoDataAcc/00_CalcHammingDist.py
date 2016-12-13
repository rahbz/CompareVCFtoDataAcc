#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-p", "--vcf_pos", dest="posFile", help="Input Position File trimmed from the VCF based on quality", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-d", "--hdf5_dile", dest="hdf5File", help="Path to HDF5 SNP matrix", type="string")
(options, args) = inOptions.parse_args()


GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

posfile = open(options.posFile, 'r')
TargetSNPs = tuple(posfile.read().rstrip().split("\n"))

NumSNP = len(TargetSNPs)
ScoreList = numpy.zeros(len(GenotypeData.accessions))
for i in range(0, len(GenotypeData.chr_regions)):
  start = GenotypeData.chr_regions[i][0]
  end = GenotypeData.chr_regions[i][1]
  chrNo = i + 1
  perChrPos = GenotypeData.positions[start:end]
  perChrTargetSNPs = [float(j.split("\t")[1]) for j in TargetSNPs if int(j.split("\t")[0]) == chrNo]
  perChrmatchedSNPind = numpy.where(numpy.in1d(perChrPos, perChrTargetSNPs))[0] + start
  ScoreList = ScoreList + numpy.sum(GenotypeData.snps[perChrmatchedSNPind,:], axis = 0)
  
print "Total number of SNPs scanned", NumSNP
#print ScoreList
# Calculate the number of SNPs in all the accessions
filetotal = open(options.file_num_snps, 'r')
numberSNPs = filetotal.read().split()
filetotal.close()
# Calculate the final probability based on the score count and total count
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = numberSNPs[numberSNPs.index(str(GenotypeData.accessions[i]))+1]
  outScore = (float(ScoreList[i]) * float(ScoreList[i]))/(NumSNP * int(numsnp))
  outfile.write(GenotypeData.accessions[i])
  outfile.write("\t")
  outfile.write("%s" % int(ScoreList[i]))
  outfile.write("\t")
  outfile.write("%s" % int(numsnp))
  outfile.write("\t")
  outfile.write("%s\n" % outScore)
outfile.close()
