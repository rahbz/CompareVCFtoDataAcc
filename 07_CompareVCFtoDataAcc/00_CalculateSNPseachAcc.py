#!/usr/bin/python
import sys
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype

hdf5file=sys.argv[1]

GenotypeData = genotype.load_hdf5_genotype_data(hdf5file)

# Calculate the number of SNPs in all the accessions
# Takes a really long time

# Calculate the final probability based on the score count and total count
#outfile = open("totalSNPsNum_1001genomes.txt", 'w')
for i in range(0, len(GenotypeData.accessions)):
#  snps = GenotypeData.snps[:, i]
#  snps[snps < 0] = 0
#  outScore = numpy.sum(snps)
  outScore = len(numpy.where(GenotypeData.snps[:, i] >= 0 )[0])
  print GenotypeData.accessions[i], "\t", outScore
#  print "Written count for", i+1, "accessions", "Accession:", GenotypeData.accessions[i], "Count:", outScore

#outfile.close()
