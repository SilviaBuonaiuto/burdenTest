# burdenTest

#### 0) Pre-processing 
Starting from the vcf file there are some pre-processing steps : 1) normalization of the variants , 2) decomposition of biallelic block substitutions into its constituent SNPs , 3) merge vcf by individuals , 4)annotation of variants and 4) creating a list of individuals
```
1) and 2) use vt (https://genome.sph.umich.edu/wiki/Vt) normalize and vt decompose
3) use bcftools merge (http://samtools.github.io/bcftools/bcftools.html)
4) annotate with VEP (https://www.ensembl.org/info/docs/tools/vep/index.html)
5) create a list of individuals from the vcf using  bcftools query -l
```
#### 1) Creating a list of SNPs 
Selecting possibly damaging and rare (MAF < 0.01) variants (qualifying variants) from the vcf file and create a tab separate list of chromosomes and positions as
```
chr1  122746558
chr3  237485948
chr14 9473625647
chr19 753454658
```
#### 2) Retrieve per-individual allele counts for all the variants in the list with DP > 10
Using SNPs and individual lists it is possible to retrieve allele counts for all the qualifying variants 
```
use vcftools (https://vcftools.github.io/index.html)
vcftools --gzvcf allsamples.$chr.vcf.gz --out $sampleid.$chr --positions positionsList.txt --min-meanDP 10 --counts --indv $sampleid  
```
#### 3) Reformat allele counts files to consider only alternate allele counts
use the script [altCounts.py](https://github.com/SilviaBuonaiuto/burdenTest/blob/main/script/altCounts.py)
```
python3 altCounts.py -i $sampleid.$chr.frq.count -id $id > $id.$chr.counts

-i path to input file (.count file from previous step)
-id sample id
```
The script produces a file (.counts) with variants positions and alternate allele, list of individuals and their counts
```
key ID  ALTcount
chr1:122746558:/A id1 1
chr1:122746558:/A id2 2
chr1:122746558:/A id3 1
chr14:9473625647:/G id2 1
chr14:9473625647:/G id4 1
chr19:753454658:/T  id1 2
chr19:753454658:/T  id3 2
chr19:753454658:/T  id4 1
```
#### 4) Merge all alternate allele count files by chromosomes and samples
Create two files one containig all allele counts for all the case samples and another one for all control samples (casesCountfile.tsv and controlCountfile.tsv) 

#### 5) Create a list of all variants and genes
The list has to be a tab separated file (genesfile.tsv)
```
key SYMBOL
chr1:122746558:/A gene1
chr3:237485948:/C gene2
chr14:9473625647:/G`gene3
chr19:753454658:/T  gene4
```
#### 6) Run burden test and generate QQplot
The script [burdenChisq.R](https://github.com/SilviaBuonaiuto/burdenTest/tree/main/script) organizes count files, run the actual burden test and generate a QQplot. Use the command line:
```
Rscript burdenChisq.R casesCountfile.tsv controlsCounfile.tsv genesfile.tsv numberofCases numberofControls

casesCountfile.tsv = file containing variants and their counts for each individual (cases)
controlCountfile.tsv = file containing variants and their counts for each individual (controls)
genesfile.tsv = file containing the list of variants and the name of the gene in wich the variants are
numberofCases = number of the cases
numberofControls = number of controls
```
