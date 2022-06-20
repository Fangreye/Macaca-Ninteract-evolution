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
The config file of mito with some explanation and it can't be feeded to the program directly. Examples were placed under ./novoplasty_example in this repository
```
Project:
-----------------------
Project name          = MMUL.IN-28474 # The prefix of outputfile
Type                  = mito
Genome Range          = 12000-22000 # The range of length of possible assembly
K-mer                 = 24 # Smaller K-mer would increase the chance of retriving the whole organelle genome while increasing the time consumption
Max memory            = 60
Extended log          = 0
Save assembled reads  = no
Seed Input            = /home/zhu46/softwares/mulatta_mito.fasta 
Extend seed directly  = no # If yes, novoplasty will use the whole sequence of seed as initiation rather than part of it. Can provide inaccurate organelle sequence.
Variance detection    = 

Dataset 1:
-----------------------
Read Length           = 101
Insert size           = 202
Platform              = illumina
Single/Paired         = PE
Combined reads        = # the locations. Used when raw reads are not paired
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

If Novoplasty successfully generate the entire organelle genome, the genome will be recorded in a fasta file with prefix 'mito{your_kmer}Circularized_assembly'. There will also be a fasta file with prefix 'mito{your_kmer}Contigs' which contains all possible contigs retrieved from raw reads. It is suggested to check the link map and concatenate them mannually.

#### Combined with GetOrganelle
If Novoplasty didn't give an acceptable output, it's still possible to retrieve the genome. GetOrganelle is a software use the strategy of mapping-then-extend, which means it can generate accurate contigs. If it didn't provide a complete organelle genome, we can feed the longest contig in the output as the seed of Novoplasty. It is also suggested to check the link map produced by SPAde to assemble a longer contig before proceeding.

```
module load StdEnv/2020  gcc/9.3.0 blast+
module load bowtie2
module load spades/3.15.3
module load scipy-stack

python3 /home/zhu46/softwares/GetOrganelle-1.7.4.1/get_organelle_from_reads.py \
-1 /home/zhu46/projects/rrg-ben/2021_Indian_rhesus/aureus_raw_data/SRR1564766_M_fasc_Maurit/SRR1564766_trim.R1.fixed.fq.gz \
-2 /home/zhu46/projects/rrg-ben/2021_Indian_rhesus/aureus_raw_data/SRR1564766_M_fasc_Maurit/SRR1564766_trim.R2.fixed.fq.gz \
-s /home/zhu46/softwares/mulatta_mito.fasta # the seed sequence \
-o /home/zhu46/scratch/Get_Organelle_test/SRR1564766 \
--continue -t 3 \
-F animal_mt -w 32 # word size when extending the contigs \
-k 25,55,75,95 # set of kmers used when detecting the connection. This and the following arguments are optimum to the mulatta dataset \
--max-n-words 5E9 \
-R 7 --reduce-reads-for-coverage 1000 \
--max-reads 9E9 --max-extending-len 9E9 &
```

### Build phylogenetic tree
The pipeline is similar to create a phylogenetic tree with nuclear genome. But mitochondrial genome is circularized, which means that the start position of a fasta record can be any site of that genome. Temporarily there is no good solution to deal with this problem. The solution here is to find the most common sequence across all samples and use it as the universal start across all samples.

## Calculate different indicators
### Length of runs of homozygosity
Long runs of homozygosity reveals that this region are under selection during a long period. 

Example code
```
#!/bin/sh
#SBATCH --job-name=brown_chr1_bcf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=1:00:00
#SBATCH --mem=20gb
#SBATCH --output=popgenWindows.%J.out
#SBATCH --error=popgenWindows.%J.err
#SBATCH --account=def-ben

module load bcftools

bcftools roh -G30 --AF-dflt 0.4 -s '34753,35495,35498,28500,38621,35723,38627,39317' # sample names in the same group inside the vcf file \
/home/zhu46/projects/rrg-ben/2021_Indian_rhesus/output.filtered.snps.removed.AB.pass.1.vep.vcf.gz | grep 'RG' > /home/zhu46/scratch/maca_800_roh/roh_chr1_brown.txt
```

### Pairwise Fst
In this step, we calculate the Fst values based on sliding windows within each chromosome between two subgroups. 

High Fst value implies that subgroups are highly differentiated. In our case, it means that the interaction between nuclear gene and mitochondrial gene is under selection in breeding and the fitness will be different in different groups. 
#### Phase and parse
Vcf files need to be phased before continue. These steps can run fast so it's better to execute them in interactive terminal.

The following pipeline comes from software 'beagle':

https://faculty.washington.edu/browning/beagle/beagle.html
```
java -Xmx28g -jar ~/softwares/beagle.29May21.d6d.jar \
gt=/home/zhu46/scratch/source/Chinese_sep_chr/population.1.vcf.gz \
out=/home/zhu46/scratch/22.windows_30k_popgen/chinese_phased_vcf/chinese_population.1.phased.vcf.gz.vcf.gz && \
```

The following pipeline comes from software 'VCF_processing':

https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing
```
python /home/zhu46/scratch/genomics_general/VCF_processing/parseVCF.py \
-i chinese_phased_vcf/chinese_population.1.phased.vcf.gz.vcf.gz | gzip \
> ./chinese_parsed_vcf/chinese_population.1.phased.parsed.vcf.gz
```
#### Calculating Fst values
The following pipeline comes from software 'genomics_general':

https://github.com/simonhmartin/genomics_general
```
#!/bin/sh
#SBATCH --job-name=popgen_chr10_{{}}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=3:00:00
#SBATCH --mem=30gb
#SBATCH --output=popgenWindows.%J.out
#SBATCH --error=popgenWindows.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 scipy-stack/2020b python/3.8.2

for i in {1..20}; do 
python /home/zhu46/scratch/genomics_general/popgenWindows.py \
-g /home/zhu46/scratch/22.windows_30k_popgen/aureus_parsed_vcf/chr${i}.filtered.cleaned.phased.parsed.vcf.gz \
-o /home/zhu46/scratch/22.aure_fst_pi_roh_redo/fst/30k/chr${i}_aureus_fascicularis.csv \
-m 1 -f phased -T 1 \
-p aureus -p fascicularis # define\
--popsFile /home/zhu46/scratch/22.aure_fst_pi_roh_redo/aureus_species.txt \
-w 30000 -s 30000 # -w define the size of window and -s define the start of next window \
--writeFailedWindows -m 10 --windType coordinate & done
```

This pipeline needs a file containing sample names and their group. The format should be 'Sample_name Group_name'.
See file 02.aureus_species.txt as an example.

### Group Nucleotide diversity {pi}
Low nucleotide diversity means that this region is under selection in recent time. 

Pi values are already contained in the output of popgenWindows.py. They are located in the 6th and 7th column.
### dNdS 
dNdS value is also an indicator of selection but it can reveal that which type of selection operates on that gene.  

#### Data preparation
We first need to build a snp-only vcf.
```
bcftools view --types snps /home/zhu46/scratch/source/Chinese_sep_chr/population.1.vcf.gz > /home/zhu46/scratch/source/Chinese_sep_chr/population.1.snp.vcf.gz
```
And an reference sequence file only containing coding regions
```
bedtools getfasta \
-fi /home/zhu46/projects/rrg-ben/2021_Indian_rhesus/2021_rhemac_v10/rheMac10.fa \
-bed /home/zhu46/scratch/111.dNdS_analysis/01.build_bed/rheMac10.dedup.bed > rheMac10_1_ref.fasta
```
Then create single sample consensus sequence.
```
for i in  SRR1024051_trim_sorted.bam  SRR2981114_trim_sorted.bam;
do cat /home/zhu46/scratch/macaca/111.dNdS_analysis/00.seq_data/rheMac10_1_ref.fasta | \
    bcftools consensus -s ${i} -M N -o /home/zhu46/scratch/macaca/111.dNdS_analysis/04.aureus/02.source_fasta/${i}.fasta \
    /home/zhu46/scratch/macaca/111.dNdS_analysis/04.aureus/01.snp_only_vcf/aureus.snp.vcf.gz & 
    [ $( jobs | wc -l ) -ge $( nproc ) ] && wait
done
```
There could be other methods but the length of coding regions should be the same across all samples.

After creating all samples
```
```

# Data analysis


# Results

# Others
## Difference between two versions of reference:


# links
1. Chinese macaca : https://academic.oup.com/gigascience/article/7/9/giy106/5079661#121467851
2. Indian macaca: https://www.science.org/doi/10.1126/science.abc6617
3. M. arctoides: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA622565/
4. M. f. aureus: https://academic.oup.com/gbe/article/13/1/evaa209/5921179
5. Novoplasty: https://github.com/ndierckx/NOVOPlasty
6. GetOrganelle: https://github.com/Kinggerm/GetOrganelle

