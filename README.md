# Macaca-Ninteract-evolution
This document record pipelines used when dealing with macaca Ninteract genes.

# Source of Data
Data used in the project come from 4 different sources:
  1. Wild caught samples in China (N=79, mapped to version 8)
  2. Laboratory strains from India (N=89, mapped to version 10)
  3. The stump-tailed macaque (N=10, mapped to version 8)
  4. A subspecies of longtail macaque, both of whose mitochondrial DNA is introgressed from a highly diverged ancestor (N=7, mapped to version 10)

Sequence file are downloaded from:
  1. Chinese macaca(PRJNA345528): https://www.ncbi.nlm.nih.gov/bioproject/PRJNA345528
  2. Indian macaca(PRJNA662298): https://www.ncbi.nlm.nih.gov/bioproject/PRJNA662298
    - Not all individuals are included in this project. Randomly sampled individuals are listed in 01.sampled_Indian_macaca.txt

All vcf files are downloaded from:
  1. Chinese: https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100484/population.vcf.gz
  2. Indian: local under /home/zhu46/projects/rrg-ben/zhu46/03.macaca_cleaned/01.original_vcf/Indian_maca_selected_vcf
  3. M. arctoides: local under /home/zhu46/projects/rrg-ben/zhu46/03.macaca_cleaned/01.original_vcf/arctoides_vcf
  4. M. f. aureus: local under /home/zhu46/projects/rrg-ben/zhu46/03.macaca_cleaned/01.original_vcf/aureus_vcf
  
References genome and annotations come from UCSC:
  1. https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.gz
  2. https://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/bigZips/genes/rheMac8.ncbiRefSeq.gtf.gz
  3. https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz
  4. https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/genes/rheMac10.ncbiRefSeq.gtf.gz
  

# Description of work flow and detailed pipeline
## Build phylogeny
In this part, we build phylogeny based on both nuclear genome and mitochondrial genome. The mitochondrial genome are not available directly so we need to recover it from the nuclear genome.

### Get mitochondrial DNA
Here we use Novoplasty to recover organelle genome. This progame use a seed squence as a start and try to connect k-mers to the seed. 

The seed used here is NC_005943.1, which can be downloaded from NCBI Database
#### Direct from Novoplasty
Assemble with Novoplasty requires a lot of memory. According to my experience, RAM should be 2~3 times of the size of raw compressed fasta file.
```
#!/bin/sh
#SBATCH --job-name=SRR1927139_MMUL.IN-28474_mito
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --time=4:00:00
#SBATCH --mem=65gb
#SBATCH --output=popgenWindows.%J.out
#SBATCH --error=popgenWindows.%J.err
#SBATCH --account=def-ben

perl /home/zhu46/softwares/NOVOPlasty-master/NOVOPlasty4.3.1.pl \
  -c /home/zhu46/scratch/maca_800_config/SRR1927139.txt # The location of configure file
```
The config file of mito with some explanation. Examples were placed under 
```
Project:
-----------------------
Project name          = MMUL.IN-28474
Type                  = mito
Genome Range          = 12000-22000
K-mer                 = 24
Max memory            = 60
Extended log          = 0
Save assembled reads  = no
Seed Input            = /home/zhu46/softwares/mulatta_mito.fasta 
Extend seed directly  = no
Variance detection    = 

Dataset 1:
-----------------------
Read Length           = 101
Insert size           = 202
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = /home/zhu46/scratch/source/fastq/SRR1927139_1.fastq.gz
Reverse reads         = /home/zhu46/scratch/source/fastq/SRR1927139_2.fastq.gz
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   = 
HP exclude list       = 
PCR-free              = 

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Output path           = /home/zhu46/scratch/maca_800_out/mito24

```

## Calculate different indicators
## Fit the model

# Results

# Others
## Difference between two versions of reference:


# links
1. Chinese macaca : https://academic.oup.com/gigascience/article/7/9/giy106/5079661#121467851
2. Indian macaca: https://www.science.org/doi/10.1126/science.abc6617
3. M. arctoides: 
4. M. f. aureus: https://academic.oup.com/gbe/article/13/1/evaa209/5921179
