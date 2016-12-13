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
inOptions.add_option("-m", "--matched_pos", dest="inFile", help="Input Matched Pos file from Rscript", type="string")
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs taken for analysis from VCF file", type="float")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")

(options, args) = inOptions.parse_args()


GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)

# Create a numpy array containing all the positions

NumSNPs = int(options.numsnps)
ScoreList = numpy.zeros(len(GenotypeData.accessions))

matchedSNP = open(options.inFile, 'r')
for matsnp in matchedSNP.readlines()[0:]:
  modLine = tuple(map(float, matsnp.rstrip().split(",")))
  ScoreList = ScoreList + numpy.sum(GenotypeData.snps[modLine,:], axis = 0)

matchedSNP.close()
    
print "Total number of SNPs scanned", NumSNPs

# Calculate the number of SNPs in all the accessions
TotnumberSNPs = open(options.file_num_snps, 'r').read().split("\n")

FinalScore = numpy.zeros(len(GenotypeData.accessions))
FinalScore = [(float(ScoreList[i])*float(ScoreList[i]))/(NumSNPs * float(TotnumberSNPs[i].split("\t")[1])) for i in range(0, len(GenotypeData.accessions))]
TopHitAccs = numpy.where(FinalScore > (max(FinalScore) - 2*numpy.std(FinalScore)))[0]
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = TotnumberSNPs[i].split("\t")[1]
  outfile.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), int(numsnp), FinalScore[i]))

outfile.close()

if len(TopHitAccs) > 1:
  outrefScore=open(options.refScore,'w')
# Get the positions where the ambiguous acc actually differ
  AccSNPs = GenotypeData_acc.snps[:, TopHitAccs]
  max= len(TopHitAccs)
  sumNumpy = numpy.sum(AccSNPs, axis = 1)
  RefPosInd = numpy.where((sumNumpy != 0) & (sumNumpy != max))[0]
  totalAccSNPs = numpy.sum(AccSNPs[RefPosInd, :], axis = 0)
  refScoreList = numpy.zeros(len(TopHitAccs))
  numMatched = 0
  for matsnp in open(options.inFile, 'r').readlines()[0:]:
    modLine = map(float, matsnp.rstrip().split(","))
    matSNPsind = RefPosInd[numpy.where(numpy.in1d(RefPosInd, modLine))]
    numMatched += len(matSNPsind)
    refScoreList = refScoreList + numpy.sum(AccSNPs[matSNPsind,], axis = 0)
  for i in range(0, len(TopHitAccs)):
    outScore = (refScoreList[i]*refScoreList[i])/(numMatched*totalAccSNPs[i])   
    outrefScore.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], int(refScoreList[i]), totalAccSNPs[i], outScore))
  outrefScore.close()
