This pipeline implements the mapping and statistics workflow of the Lundbeck Foundation GeoGenetics Centre's bioinformatics pipeline. Mapping assumes paired end sequencing, while the statistics pipeline assumes single end (or a mixture of paired and single end) reads.

# Installation 

Install from source and use the provided conda environment 
```bash 
git clone https://github.com/abigailramsoe/lundbeck-pipeline.git
cd lundbeck-pipeline
conda env create -f environment.yml
conda activate lundbeck
```

If not working on the dandy cluster, you will need to ensure that the following are installed 
* Java >1.8
* R 3.6.1
* [ANGSD](https://github.com/ANGSD/angsd)
* [Decluster/Superduper](https://github.com/ANGSD/decluster)
* [bam2prof](https://github.com/grenaud/bam2prof) 


# Usage 

## Mapping 

First, make a config file.  There is an exampe in config.yml

```yaml
ref: /maps/projects/lundbeck/data/hg38/genome.fa  # Path to fasta reference genome 
seed: 32  # Seed (-l) for bwa
a1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # Adapter 1 sequence to remove
a2: "GGAAGAGCGTCGTGTAGGGAAAGAGTGT" # Adapter 2 sequence to remove 
units: units.tsv # Units file, more on this later
outfol: /maps/projects/lundbeck/scratch/sslib_devel/mapped/ # Folder to write mapped files to 
 ```

You need a units file that tells snakemake where your fastqs are. This is a tab seperated file with a header (R1, R2, sm, lb, id, idlane) and columns:
1. R1: full path to R1 fastq
2. R2: full path to R2 fastq
3. sm: SM tag for read groups
4. lb: LB tag for read groups
5. id: ID tag for read groups
6. just noticed this is fucked up 
units file 
how to run 
how to slurm 

how to stats 
complexity stats
stats explaanaaiotn

