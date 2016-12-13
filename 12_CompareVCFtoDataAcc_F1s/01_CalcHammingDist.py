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
import pandas
from pygwas.core import genotype
import vcfnp

import scipy
import math

def nCr(n, r):
  f = math.factorial
  return f(n) / (f(r) * f(n-r))

def likeliTest(n, y):
  p = 0.999999  ### Since computing the right likelihood is troubling
  pS = float(y)/n
  a = y * scipy.log(pS/p)
  b = (n - y) * scipy.log((1-pS)/(1-p))
#  c = scipy.log(nCr(n, y))
  return(a+b)


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-i", "--vcf_file", dest="vcfFile", help="VCF file for the sample", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")

(options, args) = inOptions.parse_args()

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
num_lines = len(GenotypeData.accessions)

logging.info("Reading the VCF file")
vcf = vcfnp.variants(options.vcfFile, cache=False).view(numpy.recarray)
vcfD = vcfnp.calldata_2d(options.vcfFile, cache=False).view(numpy.recarray)

## Doubtful .... whether there should be a threshold based on just mean of std
#DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) + numpy.std(vcf.DP[numpy.where(vcf.DP > 0)[0]])
DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) * 4
logging.info("Threshold for depth is set at: %s", DPthres)

snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres))[0]
snpCHR = numpy.array(numpy.core.defchararray.replace(vcf.CHROM[snpsREQ], "Chr", ""), dtype = "int8")
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

TotalMatTarInd = numpy.zeros(0, dtype="int32")
TotalMatAccInd = numpy.zeros(0, dtype="int32")

NumMatSNPs = 0
chunk_size = 1000

for i in range(1,6):
  perchrTarPos = numpy.where(snpCHR == i)[0]
  perchrtarSNPpos = snpPOS[perchrTarPos]
  logging.info("Loaded %s chromosome positions from the position file", i)
  start = GenotypeData.chr_regions[i-1][0]
  end = GenotypeData.chr_regions[i-1][1]
  chrpositions = GenotypeData.positions[start:end]
  matchedAccInd = numpy.where(numpy.in1d(chrpositions, perchrtarSNPpos))[0] + start
  matchedTarInd = numpy.where(numpy.in1d(perchrtarSNPpos, chrpositions))[0]
  matchedTarWei = snpWEI[perchrTarPos[matchedTarInd],]
  TarGTs0 = numpy.zeros(len(matchedTarInd), dtype="int8")   ### HomRef is the first
  TarGTs1 = numpy.ones(len(matchedTarInd), dtype="int8") + 1  ## Alt is second
  TarGTs2 = numpy.ones(len(matchedTarInd), dtype="int8")     ## HomoAlt is the third in biallelic site
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  TotalMatTarInd = numpy.append(TotalMatTarInd, perchrTarPos[matchedTarInd])
  TotalMatAccInd = numpy.append(TotalMatAccInd, matchedAccInd)
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


#homref = len(numpy.where(snpGT[TotalMatTarInd] == "0/0")[0])
#alt = len(numpy.where((snpGT[TotalMatTarInd] == "0/1") & (snpDP[TotalMatTarInd] > 1))[0])
#homalt = len(numpy.where(snpGT[TotalMatTarInd] == "1/1")[0])
#hetPer = float(alt)/(homref+homalt)
#topSNP = len(numpy.where(snpDP[TotalMatTarInd] > 1)[0])
#alt = numpy.sum(snpPL[TotalMatTarInd, 1])
#hetPer = alt/float(topSNP)

reqPL = snpPL[TotalMatTarInd,][:,(0,2)]
scoreHom = numpy.sum(reqPL.min(axis=1))
hetPer = scoreHom/float(NumMatSNPs) 

#LikeLiHoods = [likeliTest(NumInfoSites[i], ScoreList[i]) for i in range(num_lines)]
#LikeLiHoods = numpy.array(LikeLiHoods).astype("float")

#TopHit = numpy.amin(LikeLiHoods)
#LikeLiHoodRatio = [LikeLiHoods[i]/TopHit for i in range(num_lines)]
#LikeLiHoodRatio = numpy.array(LikeLiHoodRatio).astype("float")
#TopHitAcc = numpy.where(LikeLiHoodRatio < 3.841)[0]
#logging.info("Number of ambiguous accessions: %s", len(TopHitAcc))

FinalScores = [float(ScoreList[i])/NumInfoSites[i] for i in range(num_lines)]
FinalScores = numpy.array(FinalScores, dtype=float)
#outfile = open(options.outFile, 'w')
#for i in range(num_lines):
#  outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], FinalScores[i], LikeLiHoods[i], LikeLiHoodRatio[i], NumMatSNPs, len(snpPOS), hetPer))
#outfile.close()

Accessions = numpy.copy(GenotypeData.accessions)
TopHitAccs = numpy.argsort(-FinalScores)[0:10]
#outfile = open(options.refScore, 'w')
for i in range(len(TopHitAccs)):
  for j in range(i+1, len(TopHitAccs)):
    p1 = GenotypeData_acc.snps[:,TopHitAccs[i]]
    p2 = GenotypeData_acc.snps[:,TopHitAccs[j]]
    gtp1 = p1[TotalMatAccInd]
    gtp2 = p2[TotalMatAccInd]
    matchedTarWEI = snpWEI[TotalMatTarInd,]
    homalt = numpy.where((gtp1 == 1) & (gtp2 == 1))[0]
    homref = numpy.where((gtp1 == 0) & (gtp2 == 0))[0]
    het = numpy.where((gtp1 != -1) & (gtp2 != -1) & (gtp1 != gtp2))[0]
    score = numpy.sum(matchedTarWEI[homalt, 2]) + numpy.sum(matchedTarWEI[homref, 0]) + numpy.sum(matchedTarWEI[het, 1])
    ScoreList = numpy.append(ScoreList, score)
    TotInfo = len(homalt) + len(homref) + len(het)
    NumInfoSites = numpy.append(NumInfoSites, TotInfo)
    acc = GenotypeData.accessions[TopHitAccs[i]] + "x" + GenotypeData.accessions[TopHitAccs[j]]
    Accessions = numpy.append(Accessions, acc)
#    outfile.write("%sx%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], GenotypeData.accessions[TopHitAccs[j]], int(score), TotInfo, float(score)/TotInfo))
#outfile.close()

LikeLiHoods = [likeliTest(NumInfoSites[i], ScoreList[i]) for i in range(len(Accessions))]
LikeLiHoods = numpy.array(LikeLiHoods).astype("float")
TopHit = numpy.amin(LikeLiHoods)
LikeLiHoodRatio = [LikeLiHoods[i]/TopHit for i in range(len(Accessions))]
LikeLiHoodRatio = numpy.array(LikeLiHoodRatio).astype("float")
TopHitAcc = numpy.where(LikeLiHoodRatio < 3.841)[0]
logging.info("Number of ambiguous accessions: %s", len(TopHitAcc))
outfile = open(options.outFile, 'w')
for i in range(len(Accessions)):
  outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (Accessions[i], int(ScoreList[i]), NumInfoSites[i], ScoreList[i]/float(NumInfoSites[i]), LikeLiHoods[i], LikeLiHoodRatio[i], NumMatSNPs, len(snpPOS), hetPer))
outfile.close()



