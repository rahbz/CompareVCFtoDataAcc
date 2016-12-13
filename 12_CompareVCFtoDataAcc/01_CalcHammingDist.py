#!/usr/bin/python
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
# module load vcfnp
import logging
import numpy
import numpy.ma
from pygwas.core import genotype
import vcfnp

import scipy
import math

def nCr(n, r):
  f = math.factorial
  return f(n) / (f(r) * f(n-r))

def likeliTest(n, y):
  p = 0.99999  ### Since computing the right likelihood is troubling
  if n > 0:
    pS = float(y)/n
    a = y * scipy.log(pS/p)
    b = (n - y) * scipy.log((1-pS)/(1-p))
#   c = scipy.log(nCr(n, y))
    return(a+b)
  elif n == y:
    return 1
  else:
    return numpy.nan


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-i", "--vcf_file", dest="vcfFile", help="VCF file for the sample", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)

logging.info("Reading the VCF file")
vcf = vcfnp.variants(options.vcfFile, cache=True).view(numpy.recarray)
vcfD = vcfnp.calldata_2d(options.vcfFile, cache=True).view(numpy.recarray)

## Doubtful .... whether there should be a threshold based on just mean of std
#DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) + numpy.std(vcf.DP[numpy.where(vcf.DP > 0)[0]])
DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) * 4
print "Threshold for depth is set at: ", DPthres
DPmean = DPthres/4

snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres))[0]
snpCHR = numpy.array(numpy.char.replace(vcf.CHROM[snpsREQ], "Chr", "")).astype("int8")
snpPOS = numpy.array(vcf.POS[snpsREQ])
snpGT = vcfD.GT[snpsREQ, 0]   ## since one sample 
snpPL = vcfD.PL[snpsREQ, 0]   
snpDP = vcf.DP[snpsREQ]

snpWEI = numpy.copy(snpPL)
snpWEI = snpWEI.astype(float)
snpWEI = snpWEI/(-10)
snpWEI = numpy.exp(snpWEI)

ScoreList = numpy.zeros(num_lines, dtype="float")
NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")

NumMatSNPs = 0
chunk_size = 1000

for i in numpy.array(GenotypeData.chrs, dtype=int):
  perchrTarPos = numpy.where(snpCHR == i)[0]
  perchrtarSNPpos = snpPOS[perchrTarPos]
  logging.info("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarWei = snpWEI[perchrTarPos[matchedTarInd],]

  TarGTs0 = numpy.zeros(len(matchedTarInd), dtype="int8")
  TarGTs1 = numpy.ones(len(matchedTarInd), dtype="int8") + 1
  TarGTs2 = numpy.ones(len(matchedTarInd), dtype="int8")
  
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  for j in range(0, len(matchedAccInd), chunk_size):
    t1001SNPs = GenotypeData.snps[matchedAccInd[j:j+chunk_size],:]
    samSNPs0 = numpy.reshape(numpy.repeat(TarGTs0[j:j+chunk_size], num_lines), (len(TarGTs0[j:j+chunk_size]),num_lines))
    samSNPs1 = numpy.reshape(numpy.repeat(TarGTs1[j:j+chunk_size], num_lines), (len(TarGTs1[j:j+chunk_size]),num_lines))
    samSNPs2 = numpy.reshape(numpy.repeat(TarGTs2[j:j+chunk_size], num_lines), (len(TarGTs2[j:j+chunk_size]),num_lines))
    tempScore0 = numpy.sum(numpy.multiply(numpy.array(t1001SNPs == samSNPs0, dtype=int).T, matchedTarWei[j:j+chunk_size,0]).T, axis=0)
    tempScore1 = numpy.sum(numpy.multiply(numpy.array(t1001SNPs == samSNPs1, dtype=int).T, matchedTarWei[j:j+chunk_size,1]).T, axis=0)
    tempScore2 = numpy.sum(numpy.multiply(numpy.array(t1001SNPs == samSNPs2, dtype=int).T, matchedTarWei[j:j+chunk_size,2]).T, axis=0)
    ScoreList = ScoreList + tempScore0 + tempScore1 + tempScore2
    if(len(TarGTs0[j:j+chunk_size]) > 1):
      NumInfoSites = NumInfoSites + len(TarGTs0[j:j+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int), axis = 0) # Number of informative sites
    elif(len(TarGTs0[j:j+chunk_size]) == 1): 
      NumInfoSites = NumInfoSites + 1 - numpy.ma.masked_less(t1001SNPs, 0).mask.astype(int)
  logging.info("Done analysing %s positions", NumMatSNPs)

logging.info("Done calculating the scores for each accession")

LikeLiHoods = [likeliTest(NumInfoSites[i], int(ScoreList[i])) for i in range(num_lines)]
LikeLiHoods = numpy.array(LikeLiHoods).astype("float")

TopHit = numpy.amin(LikeLiHoods)
LikeLiHoodRatio = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
LikeLiHoodRatio = numpy.array(LikeLiHoodRatio).astype("float")
TopHitAcc = numpy.where(LikeLiHoodRatio < 3.841)[0]
logging.info("Number of ambiguous accessions: %s", len(TopHitAcc))

outfile = open(options.outFile, 'w')
for i in range(num_lines):
  score = float(ScoreList[i])/NumInfoSites[i]
  outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], score, LikeLiHoods[i], LikeLiHoodRatio[i], NumMatSNPs, len(snpPOS), DPmean))
outfile.close()



