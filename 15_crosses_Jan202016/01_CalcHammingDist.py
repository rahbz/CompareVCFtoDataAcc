#!/usr/bin/python
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
import logging
import numpy
import numpy.ma
import pandas
from pygwas.core import genotype
import scipy
import math

def lnCr(n, r):
  f = math.factorial
  num = f(n) / (f(r) * f(n-r))
  return math.log(num)

def likeliTest(n, y):
  p = 0.99999999999999 # This is as close as 1, cannot be 1 because of log
  if n != y and n > 0:
    pS = float(y)/n
    a = y * scipy.log(pS/p)
    b = (n - y) * scipy.log((1-pS)/(1-p))
#    c = lnCr(n, y)
    return(a + b)
  elif n == y and n > 0:
    return(1)
  else:
    return(numpy.nan)

def getBins(g, binLen):
  chrlen = numpy.array((30427671, 19698289, 23459830, 18585056, 26975502))  # Arabidopsis chromosome length
  ChrIndex = numpy.zeros(0,dtype="int8")
  LenBins = numpy.zeros(0, dtype="int16")
  for i in range(5):
    start = g.chr_regions[i][0]
    end = g.chr_regions[i][1]
    chr_pos = g.positions[start:end]
    for j in range(0, chrlen[i], binLen):
      bin_pos = len(numpy.where((chr_pos >= j) & (chr_pos < j + binLen))[0])
      LenBins = numpy.append(LenBins, bin_pos)
      ChrIndex = numpy.append(ChrIndex, i+1)
  return(ChrIndex, LenBins)


#def likeliGenotype(gsnp, tarGT):

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-p", "--pos_file", dest="posFile", help="Position file removing the header from VCF", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-s", "--scoreFile", dest="scoreFile", help="Output of score files combining all the windows", type="string")

inOptions.add_option("-b", "--binLength", dest="binLen", help="Length of bins", type="int", default=300000)
inOptions.add_option("-t", "--LikeliThreshold", dest="LRthres", help="Threshold for the likelihood ratio to call ambiguous", type="float", default=3.841)


(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)
logging.info("Getting bins with %s bp length in the genome", options.binLen)
(ChrBins, PosBins) = getBins(GenotypeData, options.binLen)
# Create a numpy array containing all the positions
logging.info("Reading the position file")
targetSNPs = pandas.read_table(options.posFile, header=None, usecols=[0,1,2])

NumMatSNPs = 0
chunk_size = 1000

TotScoreList = numpy.zeros(num_lines, dtype="uint32")
TotNumInfoSites = numpy.zeros(num_lines, dtype="uint32")

outfile = open(options.outFile, 'w')
for i in range(len(PosBins)):
  start = numpy.sum(PosBins[0:i])
  end = start + PosBins[i]
  pos = GenotypeData.positions[start:end]
  perchrTarPos = numpy.where(targetSNPs[0] == ChrBins[i])[0]
  perchrtarSNPpos = numpy.array(targetSNPs[1][perchrTarPos])
  matchedAccInd = numpy.where(numpy.in1d(pos, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, pos))[0]
  matchedTarGTs = targetSNPs[2][perchrTarPos[matchedTarInd]]
  TarGTs = numpy.zeros(len(matchedTarGTs), dtype="int8")
  TarGTs[numpy.where(matchedTarGTs == "1/1")[0]] = 1
  TarGTs[numpy.where(matchedTarGTs == "0/1")[0]] = 2
#  TarGTs1 = numpy.ones(len(matchedTarGTs), dtype="int8")


  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  ## Initializing the score matrix and info sites
  ScoreList = numpy.zeros(num_lines, dtype="uint32")
  NumInfoSites = numpy.zeros(num_lines, dtype="uint32")
  for j in range(0, len(matchedAccInd), chunk_size):
    t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
    samSNPs = numpy.reshape(numpy.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
#    samSNPs1 = numpy.reshape(numpy.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
    ScoreList = ScoreList + numpy.sum(t1001SNPs == samSNPs, axis=0)
#    ScoreList = ScoreList + (numpy.sum(t1001SNPs == samSNPs, axis=0) + numpy.sum(t1001SNPs == samSNPs1, axis=0))/2
    if(len(TarGTs[j:j+chunk_size]) > 1):
      NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
    elif(len(TarGTs[j:j+chunk_size]) == 1): 
      NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
  if i % 10 == 0:
    logging.info("Done analysing %s positions", NumMatSNPs)
  
  TotScoreList = TotScoreList + ScoreList
  TotNumInfoSites = TotNumInfoSites + NumInfoSites
  likeliScore = [likeliTest(NumInfoSites[k], ScoreList[k]) for k in range(num_lines)]
  likeliScore = numpy.array(likeliScore)
  if len(likeliScore) > 0:
    TopHit = numpy.nanmin(likeliScore)
    likeliHoodRatio = [likeliScore[k]/TopHit for k in range(num_lines)]
    likeliHoodRatio = numpy.array(likeliHoodRatio)
    TopHitAcc = numpy.where(likeliHoodRatio == 1)[0]
    NumAmb = numpy.where(likeliHoodRatio < options.LRthres)[0]
    if len(TopHitAcc) == 1:
      t = TopHitAcc[0]
      score = float(ScoreList[t])/NumInfoSites[t]
      if len(numpy.where(likeliHoodRatio >= options.LRthres)[0]) > 0:
        nextLikeli = numpy.nanmin(likeliHoodRatio[numpy.where(likeliHoodRatio >= options.LRthres)[0]])
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
      else:
        nextLikeli = numpy.nanmin(likeliHoodRatio[numpy.where(likeliHoodRatio >= 1)[0]])
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
    elif len(TopHitAcc) > 1:
      if len(numpy.where(likeliHoodRatio >= options.LRthres)[0]) > 0:
        nextLikeli = numpy.nanmin(likeliHoodRatio[numpy.where(likeliHoodRatio >= options.LRthres)[0]])
        for k in range(len(TopHitAcc)):
          t = TopHitAcc[k]
          score = float(ScoreList[t])/NumInfoSites[t]
          outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))
      else:
        nextLikeli = numpy.nanmin(likeliHoodRatio[numpy.where(likeliHoodRatio >= 1)[0]])
        for k in range(len(TopHitAcc)):
          t = TopHitAcc[k]
          score = float(ScoreList[t])/NumInfoSites[t]
          outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[t], ScoreList[t], NumInfoSites[t], score, likeliScore[t], nextLikeli, len(NumAmb), i+1))

logging.info("Done calculating the scores for each accession")
outfile.close()


LikeLiHoods = [likeliTest(TotNumInfoSites[i], TotScoreList[i]) for i in range(num_lines)]
LikeLiHoods = numpy.array(LikeLiHoods).astype("float")

TopHit = numpy.nanmin(LikeLiHoods)
LikeLiHoodRatio = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
LikeLiHoodRatio = numpy.array(LikeLiHoodRatio).astype("float")

outfile = open(options.scoreFile, 'w')
for i in range(num_lines):
  score = float(TotScoreList[i])/TotNumInfoSites[i]
  outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(TotScoreList[i]), TotNumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatio[i], NumMatSNPs))
outfile.close()
