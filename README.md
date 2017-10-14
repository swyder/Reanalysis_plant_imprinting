Scripts and data accompanying Wyder et al 2017:    
Consistent Reanalysis of Genome-wide Imprinting Studies in Plants Using Generalized Linear Models Increases Concordance across Datasets
bioRxiv https://doi.org/10.1101/180745  

***
## Data folders
 
&nbsp;   | &nbsp;
-------- | ---
RAW_COUNTS | maternal/paternal raw counts of examined datasets
Informative_SNP_Positions | helper files for Classify_Alleles.py 
 
***

## Classify_Alleles.py
This script classifies mRNA-seq reads by strain. It takes a sorted BAM file, and a file containing positions of interest (SNPs where parental lines are homozygous for a different variant).  
Requires pysam (tested with pysam v0.8.4). Installation of pysam is easiest with bioconda especially on Mac OS X (see https://github.com/pysam-developers/pysam). 

### Example:
```
python Classify_Alleles.py example_sorted.bam example_Pos_Of_Interest > Counts_Alleles_SRRxxxxxxx
```

### Input: example_Pos_Of_Interest
A text file with six columns (tab delimited) describing position of interest annotated with overlapping GeneID. For reproducibility, this file should be sorted by chromosome and position (`sort -k1,1 -k2,2nr FILE`)  
Chr1	37387	37388	G	T	AT1G01060
...
1. chromosome  
2. start  
3. end  
3. reference allele  
4. non-reference allele  
5. GeneID overlapping this position  

Annotation with GeneID can be done with `bedtools intersect`. Files ready for analysis are available in folder `Informative_SNP_Positions` (Arabidopsis Ler-Col and Maize B73-Mo17). 

### Output:
Gene	MaternalReads	PaternalReads
AT1G06190	34	55
...
A text file with three columns
1. GeneID
2. reference read counts 
3. non-reference read counts
  
*NOTE:* PCR duplicates are not excluded but can be removed before using `samtools rmdup`.

***
## run_edgeR_LerCol_Pignatta.R
This script runs GLM analysis based on edgeR to identify statistically significantly imprinted genes. Requires RNA-seq of reciprocal F1 cross samples
(at least 1 per reciprocal cross, better 2-3 samples per reciprocal cross).  
  
Assumes that allelic count tables are located in the same directory with file names "Counts_Alleles_SRRxxxxxxx" where xxxxxxx is the SRR sample ID.  
  
Example for ColxLer and LerxCol samples from Pignatta et al. (2014). Not fully generalized, some things are hard-coded.

### Example:
```
Rscript run_edgeR_LerCol_Pignatta.R
```

*NOTE:* Make sure the working directory only contains count files from reciprocal crosses of 2 strains (e.g ColxLer or LerxCol). This script uses all files whose name start with 'Counts_Alleles_'


***



## Authors

* **Stefan Wyder** 

