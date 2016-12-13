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


def getAlleleBatch(GenotypeData, SNPpos):
#  Returns the allele information of all the accessions, given the position
  SNPlist = GenotypeData.get_snps_from_pos(SNPpos)
  ScoreList = numpy.zeros(len(GenotypeData.accessions))
  for i in range(0, len(SNPlist[0])):
    if (GenotypeData.get_chr_pos_from_index(SNPlist[0][i]) == SNPpos[i]):
      for j in range(0, len(SNPlist[1][i])):
        ScoreList[j] = ScoreList[j] + int(SNPlist[1][i][j])
  return ScoreList

def AddScoreList(ScoreList, tempScore):
# Add the temporary score list to the new one
  newScoreList = numpy.zeros(len(ScoreList))
  for i in range(0, len(ScoreList)):
    newScoreList[i] = ScoreList[i] + tempScore[i]
#    print i, "\t", newScoreList[i], "\t", tempScore[i], "\t", ScoreList[i]
  return newScoreList

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-i", "--input_vcf", dest="inFile", help="Input VCF file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-n", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-t", "--qual_threshold", dest="qual", help="Quality threshold for the SNPs", default = 100, type="float")
inOptions.add_option("-e", "--error_rate", dest="error", help="Include error rate while calculating", default = 0.001, type="float")
#inOptions.add_option("-s", "--score_list", dest="score_list", help="Printing ScoreList", type="string")
inOptions.set_defaults(error=0.001, qual=100)
(options, args) = inOptions.parse_args()


inputVCFfile = open(options.inFile, 'r')
GenotypeData = genotype.load_hdf5_genotype_data('/lustre/scratch/users/rahul.pisupati/all_chromosomes_binary.hdf5')

# Create a numpy array containing all the positions
NumSNP = 0
ScoreList = numpy.zeros(len(GenotypeData.accessions))
#Modified to load the entire positions of vcf into a array
SNPpos = []
for vcfLine in inputVCFfile.readlines()[0:]:
  if (vcfLine[0][0] != '#'):
    if(float(vcfLine.split()[5]) > options.qual and len(vcfLine.split()[3]) == 1 and len(vcfLine.split()[4]) == 1):
      SNPpos.append((int(vcfLine.split()[0].replace("Chr", "")), int(vcfLine.split()[1])))
      NumSNP += 1
      if (NumSNP % 1 == 0):
        tempScoreList = numpy.zeros(len(GenotypeData.accessions))
#        print NumSNP
        tempScoreList = getAlleleBatch(GenotypeData, SNPpos)
#        print tempScoreList
        ScoreList = AddScoreList(ScoreList, tempScoreList)
        SNPpos = []
if (len(SNPpos) != 0):
  tempScoreList = numpy.zeros(len(GenotypeData.accessions))
  tempScoreList = getAlleleBatch(GenotypeData, SNPpos)
  ScoreList = AddScoreList(ScoreList, tempScoreList)
  SNPpos = []

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
#  print GenotypeData.accessions[i], "\t", ScoreList[i],"\t", outScore
  outfile.write(GenotypeData.accessions[i])
  outfile.write("\t")
  outfile.write("%s\n" % outScore)
outfile.close()
