#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
import pandas
from pygwas.core import genotype


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-p", "--pos_file", dest="posFile", help="Position file removing the header from VCF", type="string")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")

(options, args) = inOptions.parse_args()


GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)

# Create a numpy array containing all the positions
targetSNPs = pandas.read_table(options.posFile, header=None)
NumSNPs = len(targetSNPs)
ScoreList = numpy.zeros(len(GenotypeData.accessions))
TotMatchedSNPind = numpy.zeros(0, dtype="uint32")
NumMatSNPs = 0
for i in range(1,6):
  perchrtarSNPpos = targetSNPs[1][numpy.where(targetSNPs[0] == i)[0]]
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedSNPind = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  TotMatchedSNPind = numpy.append(TotMatchedSNPind, matchedSNPind)
  NumMatSNPs = NumMatSNPs + len(matchedSNPind)
  for j in range(0, len(matchedSNPind), 1000):
    subSNPs = GenotypeData.snps[matchedSNPind[j:j+1000],:]
    subSNPs[subSNPs < 0] = 0
    ScoreList = ScoreList + numpy.sum(subSNPs, axis = 0)

print "Total number of SNPs scanned:", NumSNPs, "matched ", NumMatSNPs

# Calculate the number of SNPs in all the accessions
TotnumberSNPs = open(options.file_num_snps, 'r').read().split("\n")

FinalScore = numpy.zeros(len(GenotypeData.accessions))
FinalScore = [(float(ScoreList[i])*float(ScoreList[i]))/(NumMatSNPs * float(TotnumberSNPs[i].split("\t")[1])) for i in range(0, len(GenotypeData.accessions))]
TopHitAccs = numpy.where(FinalScore > (max(FinalScore) - 2*numpy.std(FinalScore)))[0]
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = TotnumberSNPs[i].split("\t")[1]
  outfile.write("%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), int(numsnp), FinalScore[i], NumMatSNPs))

outfile.close()

if len(TopHitAccs) > 1:
  outrefScore=open(options.refScore,'w')
# Get the positions where the ambiguous acc actually differ
  AccSNPs = GenotypeData_acc.snps[:, TopHitAccs]
  totAccSNPs = numpy.copy(AccSNPs)
  max= len(TopHitAccs)
  AccSNPs[AccSNPs == -1] = max + 100
  sumNumpy = numpy.sum(AccSNPs, axis = 1)
#  This doesnt work because there are -1s in the SNP matrix which screw up the above
  RefPosInd = numpy.where((sumNumpy > 0) & (sumNumpy < max))[0]
  totAccSNPs[totAccSNPs < 0] = 0
  
  totalAccSNPs = numpy.sum(totAccSNPs[RefPosInd, :], axis = 0)
  refScoreList = numpy.zeros(len(TopHitAccs))
  refMatSNPs = RefPosInd[numpy.where(numpy.in1d(RefPosInd, TotMatchedSNPind))[0]]
#  print len(refMatSNPs)
  for i in range(0, len(refMatSNPs), 1000):
    subsampleSNPs = totAccSNPs[tuple(refMatSNPs[i:i+1000]),]
    subsampleSNPs[subsampleSNPs < 0] = 0
    refScoreList = refScoreList + numpy.sum(subsampleSNPs, axis = 0)

  for i in range(0, len(TopHitAccs)):
    outScore = (refScoreList[i]*refScoreList[i])/(len(refMatSNPs)*totalAccSNPs[i])
    outrefScore.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], int(refScoreList[i]), totalAccSNPs[i], outScore))
  outrefScore.close()
