#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas
import logging
import numpy
from pygwas.core import genotype
import scipy

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--rareAlleleFreq", dest="allelFreq", help="Allele frequency to consider as rare allele", default=0.05, type="float")

#inOptions.add_option("-s", "--error_rate", dest="error", help="Maximum score which is considered to be for top hit accession", default=0.98, type="float")

(options, args) = inOptions.parse_args()


logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)
NumAcc = len(GenotypeData.accessions)

snps = GenotypeData.get_snps_iterator(is_chunked=True, chunk_size=1000)
chunk_i = 0
NumRareAllele = numpy.zeros(NumAcc)
InfoPOS = numpy.zeros(NumAcc)
logging.info("Starting the calculation")
for snp in snps:
  
  chunk_i = chunk_i + 1
  snps_array = numpy.array(snp)
  info_array = numpy.copy(snps_array)
  snps_array[snps_array == -1] = 0
  info_array[info_array == 0] = 1
  info_array[info_array == -1] = 0
  AlleFreq = numpy.array(numpy.sum(snps_array, axis= 1), dtype=float)/numpy.sum(info_array, axis= 1)
  rareInd = numpy.where(AlleFreq < options.allelFreq)[0]
  
  NumRareAllele = NumRareAllele + numpy.sum(snps_array[rareInd, :], axis=0)
  InfoPOS = InfoPOS + numpy.sum(info_array, axis=0)
  logging.info("Done the chunk %s", chunk_i)

RarePer = numpy.array(NumRareAllele, dtype=float)/InfoPOS
logging.info("Writing the Rareallel % frequency into a file")

numpy.savetxt(options.outFile, RarePer)

