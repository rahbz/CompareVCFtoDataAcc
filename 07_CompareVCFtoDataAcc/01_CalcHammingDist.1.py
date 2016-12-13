#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
import logging
import numpy
import pandas
from pygwas.core import genotype
import scipy

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-p", "--pos_file", dest="posFile", help="Position file removing the header from VCF", type="string")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)

# Create a numpy array containing all the positions
logging.info("Reading the position file")
targetSNPs = pandas.read_table(options.posFile, header=None, usecols=[0,1,2])

ScoreList = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")

TotMatchedAccInd = numpy.zeros(0, dtype="uint32")
TotMatchedTarGTs = numpy.zeros(0, dtype="uint32")
NumMatSNPs = 0
chunk_size = 1000
num_lines = len(GenotypeData.accessions)

for i in range(1,6):
  perchrTarPos = numpy.where(targetSNPs[0] == i)[0]
  perchrtarSNPpos = numpy.array(targetSNPs[1][perchrTarPos])
  logging.info("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarGTs = targetSNPs[2][perchrTarPos[matchedTarInd]]
  TarGTs = numpy.zeros(len(matchedTarGTs), dtype="int8")
  TarGTs[numpy.where(matchedTarGTs != "0/0")[0]] = 1
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  logging.debug("%s", ScoreList)
  for j in range(0, len(matchedAccInd), chunk_size):
    t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
    samSNPs = numpy.reshape(numpy.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
    tempBool = t1001SNPs == samSNPs
    tempBool = tempBool.astype(int)
    tempScore  = numpy.array(numpy.sum(tempBool, axis=0), dtype="uint16")
    ScoreList = ScoreList + tempScore
  TotMatchedAccInd = numpy.append(TotMatchedAccInd, matchedAccInd)
  TotMatchedTarGTs = numpy.append(TotMatchedTarGTs, TarGTs)
  logging.info("Done analysing %s positions", NumMatSNPs)

logging.info("Done calculating the scores for each accession")

#ScoreList = [len(numpy.where(TotMatchedTarGTs == GenotypeData_acc.snps[:, j][TotMatchedAccInd])[0]) for j in range(0, len(GenotypeData.accessions))] 
# 1.66 second for each iteration == 1.66 * 1135 seconds that is huge

#  for j in range(0, len(matchedAccInd), 1000):
#    t1001SNPs = scipy.array(numpy.copy(GenotypeData.snps[matchedAccInd[j:j+1000],:]))
#    t1001SNPs = t1001SNPs.astype("float")
#    t1001SNPs[t1001SNPs < 0] = 0.5
#    SNPs = scipy.mat(t1001SNPs * 2.0 - 1.0)
#    samSNPs = scipy.array(numpy.copy(TarGTs[j:j+1000]))
#    tSNPs = scipy.array(samSNPs * 2.0 - 1.0)
#    ScoreList = ScoreList + numpy.array(tSNPs * SNPs)[0]

TotnumberSNPs = open(options.file_num_snps, 'r').read().split("\n")

FinalScore = numpy.zeros(len(GenotypeData.accessions))
FinalScore = [float(ScoreList[i]) / float(TotnumberSNPs[i].split("\t")[1]) for i in range(0, len(GenotypeData.accessions))]
#FinalScore = [(float(ScoreList[i]) * float(ScoreList[i]))/(NumMatSNPs * float(TotnumberSNPs[i].split("\t")[1])) for i in range(0, len(GenotypeData.accessions))]
TopHitAccs = numpy.where(FinalScore > (max(FinalScore) - 2*numpy.std(FinalScore)))[0]
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = TotnumberSNPs[i].split("\t")[1]
  outfile.write("%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), int(numsnp), FinalScore[i], NumMatSNPs))

outfile.close()
  
#if len(TopHitAccs) > 1:
#  outrefScore=open(options.refScore,'w')
## Get the positions where the ambiguous acc actually differ
#  AccSNPs = GenotypeData_acc.snps[:, TopHitAccs]
#  totAccSNPs = numpy.copy(AccSNPs)
#  max= len(TopHitAccs)
#  AccSNPs[AccSNPs == -1] = max + 100
#  sumNumpy = numpy.sum(AccSNPs, axis = 1)
#  RefPosInd = numpy.where((sumNumpy > 0) & (sumNumpy < max))[0]
## This is to get the positions where the access actually differ in
#  
#  totAccSNPs[totAccSNPs < 0] = 0
#  totalAccSNPs = numpy.sum(totAccSNPs[RefPosInd, :], axis = 0)
#  refScoreList = numpy.zeros(len(TopHitAccs))
#  refMatSNPs = RefPosInd[numpy.where(numpy.in1d(RefPosInd, TotMatchedSNPind))[0]]
#
##  print len(refMatSNPs)
#  for i in range(0, len(refMatSNPs), 1000):
#    subsampleSNPs = totAccSNPs[tuple(refMatSNPs[i:i+1000]),]
#    subsampleSNPs[subsampleSNPs < 0] = 0
#    refScoreList = refScoreList + numpy.sum(subsampleSNPs, axis = 0)
#
#  for i in range(0, len(TopHitAccs)):
#    outScore = (refScoreList[i]*refScoreList[i])/(len(refMatSNPs)*totalAccSNPs[i])
#    outrefScore.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], int(refScoreList[i]), totalAccSNPs[i], outScore))
#  outrefScore.close()
