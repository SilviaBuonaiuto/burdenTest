# burdenTest

#### Required tools
vt (https://genome.sph.umich.edu/wiki/Vt)
bcftools (http://samtools.github.io/bcftools/bcftools.html)
VEP (https://www.ensembl.org/info/docs/tools/vep/index.html)
vcftools (https://vcftools.github.io/index.html)

#### 0) Pre-processing. 
1) normalization of the variants in vcf. 
2) decomposition of biallelic block substitutions into its constituent SNPs.   
3) merge vcf by individuals.  
4) annotation of variants. 
5) creating a list of individuals. 
```
1)vt normalize -n \
-r reference_genome.fa \
$sample.$chr.vcf.gz > $sample.$chr.normalized.vcf

2)vt decompose_blocksub $sample.$chr.normalized.vcf.gz > $sample.$chr.deco.vcf

3)bcftools merge -O z \
-o allsamples.$chr.vcf.gz \
sample1.$chr.deco.vcf.gz sample2.$chr.deco.vcf.gz ... sampleN.$chr.deco.vcf.gz

4)VEP

5)bcftools query -l \
allsamples.$chr.vcf.gz > samples.txt
```
#### 1) Creating a list of SNPs 
From the annotated VCF select qualifying variants (possibly damaging and rare (MAF < 0.01)) and  write it as a bed file (chromosome and position into tab separate file)

#### 2) Retrieve per-individual allele counts for all the variants with DP > 10 in the list
```
vcftools --gzvcf allsamples.$chr.vcf.gz /
--out $sampleid.$chr /
--positions qualifying_variants.bed /
--min-meanDP 10 /
--counts /
--indv $sampleid  
```
#### 3) Reformat allele counts files to consider only alternate allele counts
run script [altCounts.py](https://github.com/SilviaBuonaiuto/burdenTest/blob/main/script/altCounts.py)
```
python3 altCounts.py -i $sampleid.$chr.frq.count /
-id $id > $sampleid.$chr.counts
```
#### 4) Merge counts of all samples and chromosomes in one file (https://github.com/SilviaBuonaiuto/burdenTest/blob/main/test/qualifying_variants.idv.counts)

#### 5) Create file containing variants and counts for control samples
As the previous step, create a file with variants and counts for control samples (the ones not carrying qualifying variants)(https://github.com/SilviaBuonaiuto/burdenTest/blob/main/test/control_samples.counts) 

#### 6) Create a list of all variants and genes
Create file with alleles and genes annotations(https://github.com/SilviaBuonaiuto/burdenTest/blob/main/test/annotated_genes.tsv)

#### 7) Run burden test and generate QQplot
The script [burdenChisq.R](https://github.com/SilviaBuonaiuto/burdenTest/tree/main/script) organizes count files, run the actual burden test and generate a QQplot. Use the command line:
```
Rscript burdenChisq.R qualifying_variants.idv.counts control_samples.counts annotated_genes.tsv numberofCases numberofControls

numberofCases = number of the cases samples
numberofControls = number of controls samples
```
