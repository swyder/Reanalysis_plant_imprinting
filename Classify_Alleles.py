#!/usr/bin/python

# prints reference (ref_) and non-reference (nef_) reads counts per gene (counting only at most SNP per read)
#
# Usage: python Count_Variants.py example.bam example_Pos_Of_Interest > Counts_Alleles_SRRxxxxxxx
#
# 2015-09-01 Stefan Wyder


from sys import argv
import pysam
from collections import defaultdict
#from collections import OrderedDict


#######################################################
# Set quality Cut-offs
#######################################################

MinBaseQual = 20
MinMappingQual = 20




BamFile = argv[1]
SnpFile = argv[2]

ReadCounted = defaultdict(list)
countMat = 0
countPat = 0
SNPpositions = []
SNPinfo = defaultdict(list)
maternalReads = {}
paternalReads = {}


# Read in SNP file containing positions of interest
snpPos = open(SnpFile, "rb")
for line in snpPos.readlines():
    line = line.strip()
    Chrom, Pos, Pos2, MatAllele, PatAllele, OverlappingGene = line.split("\t")
    SNPpositions.append([Chrom, int(Pos)])
    SNPinfo[(Chrom, int(Pos))] = line
 

# Read and parse BAM file
samfile = pysam.Samfile(BamFile, "rb")
for QueryChrom, QueryPos in SNPpositions:
    OverlappingGene = ''
    for pileupcolumn in samfile.pileup(QueryChrom, QueryPos, QueryPos+1, truncate = True, stepper = 'all'):
        if pileupcolumn.pos == QueryPos:	# Positions are 0-based in pysam in contrast to samtools!
            Chrom, Pos, Pos2, MatAllele, PatAllele, OverlappingGene
            Chrom, Pos, Pos2, MatAllele, PatAllele, OverlappingGene = SNPinfo.get((QueryChrom, QueryPos), "N").split("\t")
            countMat = 0
            countPat = 0
            basecount = {}
            for pileupread in pileupcolumn.pileups:
                Alignment = pileupread.alignment
                ReadName = Alignment.query_name
                if pileupread.query_position == None:
                    continue
                BaseQualAsciiCode = ord(Alignment.qual[pileupread.query_position]) - 33   # Ascii code
                MappingQual = Alignment.mapping_quality
                if not pileupread.is_del and not pileupread.is_refskip and not Alignment.is_secondary and not \
                                ReadName in ReadCounted[OverlappingGene] and BaseQualAsciiCode >= MinBaseQual and MappingQual >= MinMappingQual:
                    Allele = Alignment.query_sequence[pileupread.query_position]
                    if Allele == PatAllele:
                        countPat += 1
                        ReadCounted[OverlappingGene].append(ReadName)
                    elif Allele == MatAllele:
                        countMat += 1
                        ReadCounted[OverlappingGene].append(ReadName)
                    for base in ['A', 'C', 'G', 'T']:
                        if Allele == base:
                            basecount[base] = basecount.get(base, 0) + 1
            maternalReads[OverlappingGene] = maternalReads.get(OverlappingGene, 0) + countMat
            paternalReads[OverlappingGene] = paternalReads.get(OverlappingGene, 0) + countPat
samfile.close()


# Loop through all genes and print maternal/paternal counts
print "\t".join(["Gene", "Ref_Reads", "Nef_Reads"])
for gene in maternalReads:
    if maternalReads[gene] > 0 or paternalReads[gene] > 0:
        print "\t".join([gene, str(maternalReads[gene]), str(paternalReads[gene])])
