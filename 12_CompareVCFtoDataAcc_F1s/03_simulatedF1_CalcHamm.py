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
import random

def nCr(n, r):
  f = math.factorial
  return f(n) / (f(r) * f(n-r))

def likeliTest(n, y):
  p = 0.99999  ### Since computing the right likelihood is troubling
  pS = float(y)/n
  a = y * scipy.log(pS/p)
  b = (n - y) * scipy.log((1-pS)/(1-p))
#  c = scipy.log(nCr(n, y))
  return(a+b)

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-a", "--parent1", dest="p1", help="Accession of the parent 1", type="string")
inOptions.add_option("-b", "--parent2", dest="p2", help="Accession of the parent 2", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-n", "--numSNPs", dest="numSNPs", help="Number of SNPs to be sampled from the simulated F1s", type="int")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)

# Create a numpy array containing all the positions
logging.info("Getting the segregating sites between two parents")

p1ind = numpy.where(GenotypeData.accessions == options.p1)[0]
p2ind = numpy.where(GenotypeData.accessions == options.p2)[0]
p1 = GenotypeData_acc.snps[:,p1ind]
p2 = GenotypeData_acc.snps[:,p2ind]

f1 = p1 + p2
f1[numpy.where((p1 == -1) & (p2 == -1))[0]] = -1
f1[f1 == 2] = 4
f1[f1 == 1] = 2
f1[f1 == 4] = 1
ChrPos = numpy.ones(len(GenotypeData.positions), dtype="int8")
for i in range(5):
  start = GenotypeData.chr_regions[i][0]
  end = GenotypeData.chr_regions[i][1]
  ChrPos[start:end] = i + 1


a = ChrPos[numpy.where(f1 >= 0)[0]]
b = GenotypeData.positions[numpy.where(f1 >= 0)[0]]
c = f1[numpy.where(f1 >= 0)[0]]

allSNPs = pandas.DataFrame({0: a, 1: b, 2: c})
reqSNPs = numpy.sort(random.sample(range(allSNPs.shape[0]), options.numSNPs))

targetSNPs = allSNPs.loc[reqSNPs, :]

ScoreList = numpy.zeros(num_lines, dtype="uint32")
NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")

NumMatSNPs = 0
chunk_size = 1000

for i in range(1,6):
  perchrTarPos = numpy.where(targetSNPs[0] == i)[0]
  perchrtarSNPpos = numpy.array(targetSNPs[1][perchrTarPos])
  logging.info("Loaded %s chromosome positions", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  TarGTs = numpy.array(targetSNPs[2][perchrTarPos[matchedTarInd]], dtype="int8")
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  for j in range(0, len(matchedAccInd), chunk_size):
    t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
    samSNPs = numpy.reshape(numpy.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
    ScoreList = ScoreList + numpy.sum(t1001SNPs == samSNPs, axis=0)
    if(len(TarGTs[j:j+chunk_size]) > 1):
      NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
    elif(len(TarGTs[j:j+chunk_size]) == 1): 
      NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
  logging.info("Done analysing %s positions", NumMatSNPs)

logging.info("Done calculating the scores for each accession")

LikeLiHoods = [likeliTest(NumInfoSites[i], ScoreList[i]) for i in range(num_lines)]
LikeLiHoods = numpy.array(LikeLiHoods).astype("float")

TopHit = numpy.amin(LikeLiHoods)
LikeLiHoodRatio = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
LikeLiHoodRatio = numpy.array(LikeLiHoodRatio).astype("float")
TopHitAcc = numpy.where(LikeLiHoodRatio < 3.14)[0]
logging.info("Number of ambiguous accessions: %s", len(TopHitAcc))

outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  score = float(ScoreList[i])/NumInfoSites[i]
  outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatio[i], NumMatSNPs))
outfile.close()
