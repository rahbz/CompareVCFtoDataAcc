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
#GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
#num_lines = len(GenotypeData.accessions)


logging.info("Reading the VCF file")
vcf = vcfnp.variants(options.vcfFile, cache=True).view(numpy.recarray)
vcfD = vcfnp.calldata_2d(options.vcfFile, cache=True).view(numpy.recarray)


## Doubtful .... whether there should be a threshold based on just mean of std
#DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) + numpy.std(vcf.DP[numpy.where(vcf.DP > 0)[0]])
DPthres = numpy.mean(vcf.DP[numpy.where(vcf.DP > 0)[0]]) * 4
print "Threshold for depth is set at: ", DPthres

#snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0))[0]
snpsREQ = numpy.where((vcfD.is_called[:,0]) & (vcf.QUAL > 30) & (vcf.DP > 0) & (vcf.DP < DPthres))[0]
snpCHR = numpy.array(numpy.chararray.replace(vcf.CHROM[snpsREQ], "Chr", "")).astype("int8")
snpPOS = numpy.array(vcf.POS[snpsREQ])
snpGT = vcfD.GT[snpsREQ, 0]   ## since one sample 
snpPL = vcfD.PL[snpsREQ, 0]
snpDP = vcf.DP[snpsREQ]

snpWEI = numpy.copy(snpPL)
snpWEI = snpWEI.astype(float)
snpWEI = snpWEI/(-10)
snpWEI = numpy.exp(snpWEI)

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
  NumMatSNPs = NumMatSNPs + len(matchedAccInd)
  TotalMatTarInd = numpy.append(TotalMatTarInd, perchrTarPos[matchedTarInd])
  TotalMatAccInd = numpy.append(TotalMatAccInd, matchedAccInd)


#homref = numpy.where(snpGT[TotalMatTarInd] == "0/0")[0]
#alt = numpy.where(snpGT[TotalMatTarInd] == "0/1")[0]
#homalt =  numpy.where(snpGT[TotalMatTarInd] == "1/1")[0]
#lhomoref = numpy.sum(snpPL[TotalMatTarInd,][homref, 0])
#lalt = numpy.sum(snpPL[TotalMatTarInd,][alt, 1])
#lhomoalt = numpy.sum(snpPL[TotalMatTarInd,][homalt, 2])

#reqPL = snpPL[TotalMatTarInd,][:,(0,2)]
#scoreHom = numpy.sum(numpy.multiply(reqPL.min(axis=1), snpDP[TotalMatTarInd]))
#print scoreHom, "\t", NumMatSNPs, "\t", len(snpsREQ), float(scoreHom)/numpy.sum(snpDP[TotalMatTarInd])
#print Phomref, "\t", Palt, "\t", Phomalt, "\t", Palt/(Phomref+Phomalt), "\t", len(snpsREQ)

#nAlt = len(numpy.where(snpGT[TotalMatTarInd] == "0/1")[0])
#HetPer = nAlt/float(NumMatSNPs)
#print nAlt, "\t", HetPer

pHet = numpy.prod(snpWEI[TotalMatTarInd, 1], dtype="float64")
pHomo = numpy.prod(numpy.sum(snpWEI[TotalMatTarInd,][:,(0,2)], axis = 1), dtype=float)
print pHet, "\t", pHomo#, "\t", float(pHet)/(pHet+pHomo)

