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
import vcfnp

def nCr(n, r):
  f = math.factorial
  return f(n) / (f(r) * f(n-r))

def likeliTest(n, y):
  p = 0.9999  ### Since computing the right likelihood is troubling
  if n > 0:
    pS = float(y)/n
    a = y * scipy.log(pS/p)
    b = (n - y) * scipy.log((1-pS)/(1-p))
#   c = scipy.log(nCr(n, y))
    return(a+b)
  else:
    return numpy.nan


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-p", "--pos_file", dest="posFile", help="Position file removing the header from VCF", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-v", "--the1001VCF", dest="the1001VCF", help="Path to 1001 genomes VCF", type="string")
inOptions.add_option("-a", "--accession", dest="accToCheck", help="Accession to check to get the position information", type="string")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
VCF1001 = vcfnp.variants(options.the1001VCF, cache=True).view(numpy.recarray)
accID = numpy.where(GenotypeData.accessions == options.accToCheck)[0]

num_lines = len(GenotypeData.accessions)

# Create a numpy array containing all the positions
allTargetSNPs = pandas.read_table(options.posFile, header=None, usecols=[0,1,2])
logging.info("Done reading position file")

the1001chrs = numpy.zeros(GenotypeData.num_snps, dtype="int8")
for i in range(len(GenotypeData.chrs)):
  the1001chrs[GenotypeData.chr_regions[i][0]:GenotypeData.chr_regions[i][1]] = i + 1
AccSNPs = GenotypeData_acc.snps[:,accID[0]]
nonC_ref = numpy.where((AccSNPs == 0) & (VCF1001.REF != "C"))[0]
nonC_alt = numpy.where((AccSNPs == 1) & (VCF1001.ALT != "C"))[0]
nonC_ind = numpy.sort(numpy.append(nonC_ref, nonC_alt))
nonC_chr = the1001chrs[nonC_ind]
nonC_pos = GenotypeData.positions[nonC_ind]

logging.info("Filtering on the1001SNP matrix")
### Club the columns and try to filter those sites out before taking them into consideration
npfilterSNP = numpy.core.records.fromarrays(numpy.array(numpy.column_stack((nonC_chr, nonC_pos)).T, dtype="int32"))
nptargetSNP = numpy.core.records.fromarrays(numpy.column_stack((numpy.array(allTargetSNPs[0], dtype="int32"), numpy.array(allTargetSNPs[1], dtype="int32"))).T)
filterInd = numpy.where(numpy.in1d(nptargetSNP, npfilterSNP))[0]

targetSNPs = allTargetSNPs.loc[filterInd]
logging.info("Filtered the required positions from the position file")

ScoreList = numpy.zeros(num_lines, dtype="uint32")
NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")

TotMatchedAccInd = numpy.zeros(0, dtype="uint32")
TotMatchedTarGTs = numpy.zeros(0, dtype="uint32")

NumMatSNPs = 0
chunk_size = 1000

for i in range(1,6):
  perchrTarPos = numpy.where(numpy.array(targetSNPs[0]) == i)[0]
  perchrtarSNPpos = numpy.array(targetSNPs[1])[perchrTarPos]
  logging.info("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarGTs = numpy.array(targetSNPs[2])[perchrTarPos[matchedTarInd]]
  TarGTs = numpy.zeros(len(matchedTarGTs), dtype="int8")
  TarGTs[numpy.where(matchedTarGTs == "1/1")[0]] = 1
#  TarGTs[numpy.where(matchedTarGTs == "0/1")[0]] = 2
  TarGTs1 = numpy.ones(len(matchedTarGTs), dtype="int8")
  TarGTs1[numpy.where(matchedTarGTs == "0/0")[0]] = 0

  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  for j in range(0, len(matchedAccInd), chunk_size):
    t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
    samSNPs = numpy.reshape(numpy.repeat(TarGTs[j:j+chunk_size], num_lines), (len(TarGTs[j:j+chunk_size]),num_lines))
#    ScoreList = ScoreList + numpy.sum(t1001SNPs == samSNPs, axis=0)
    samSNPs1 = numpy.reshape(numpy.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
    ScoreList = ScoreList + (numpy.sum(t1001SNPs == samSNPs, axis=0) + numpy.sum(t1001SNPs == samSNPs1, axis=0))/float(2)
    
    if(len(TarGTs[j:j+chunk_size]) > 1):
      NumInfoSites = NumInfoSites + len(TarGTs[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
    elif(len(TarGTs[j:j+chunk_size]) == 1): 
      NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
  TotMatchedAccInd = numpy.append(TotMatchedAccInd, matchedAccInd)
  TotMatchedTarGTs = numpy.append(TotMatchedTarGTs, TarGTs)
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

