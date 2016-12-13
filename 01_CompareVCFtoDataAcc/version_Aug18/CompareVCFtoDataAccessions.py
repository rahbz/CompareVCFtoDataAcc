#!/usr/bin/python
import sys
import math
from subprocess import call
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype

def getAllele(GenotypeData, chr, position):
#  Returns the allele information of all the accessions, given the position
  SearchPos = [(float(chr),float(position))]
  SNPlist = GenotypeData.get_snps_from_pos(SearchPos)
  if (GenotypeData.get_chr_pos_from_index(SNPlist[0]) == SearchPos[0]):
    return SNPlist[1][0]
  else:
    return '0'

def GetSNPsnum(GenotypeData):
# Calulates the number of SNPs of particular accession from reference
  SNPsnum = numpy.zeros(len(GenotypeData.accessions))
  for i in range(0, len(GenotypeData.accessions)):
    SNPsnum[i] = sum(GenotypeData.snps[:,i])
  return SNPsnum

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-e", "--error_rate", dest="error", help="Include error rate while calculating", default = 0.001, type="float")
inOptions.add_option("-t", "--qual_threshold", dest="qual", help="Quality threshold for the SNPs", default = 100, type="float")
inOptions.add_option("-i", "--input_vcf", dest="inFile", help="Input VCF file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.set_defaults(error=0.001, qual=100)
(options, args) = inOptions.parse_args()


inputVCFfile = open(options.inFile, 'r')
GenotypeData = genotype.load_hdf5_genotype_data('/lustre/scratch/users/rahul.pisupati/all_chromosomes_binary.hdf5')

# Create a numpy array containing all the positions
ScoreList = numpy.zeros(len(GenotypeData.accessions))
NumSNP = 0
CheckStatus = 0
for vcfLine in inputVCFfile.readlines()[0:]:
  if (vcfLine[0][0] != '#'):
    if(float(vcfLine.split()[5]) > options.qual and len(vcfLine.split()[3]) == 1 and len(vcfLine.split()[4]) == 1):
      dataSNPlist = getAllele(GenotypeData, vcfLine.split()[0].replace("Chr", ""), vcfLine.split()[1])
      NumSNP += 1
      if (dataSNPlist != '0'):
        for i in range(0, len(GenotypeData.accessions)):
          ScoreList[i] += dataSNPlist[i]
      if (NumSNP % 1000 == 0):
        print "Scanned", NumSNP, "SNPs in the VCF file"

print "Total number of SNPs scanned", NumSNP
# Calculate the number of SNPs in all the accessions
# Takes a really long time
numSNPacc = GetSNPsnum(GenotypeData)
#numSNPacc = numpy.ones(len(GenotypeData.accessions))

# Calculate the final probability based on the score count and total count
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  outScore = (float(ScoreList[i]) * float(ScoreList[i]))/(NumSNP * numSNPacc[i])
  outfile.write(GenotypeData.accessions[i])
  outfile.write("\t")
  outfile.write("%s\n" % outScore)

