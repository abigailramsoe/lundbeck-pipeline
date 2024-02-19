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

### Config file 

First, make a config file.  There is an exampe in `config.yml`

```yaml
ref: /maps/projects/lundbeck/data/hg38/genome.fa  # Path to fasta reference genome 
seed: 32  # Seed (-l) for bwa
a1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # Adapter 1 sequence to remove
a2: "GGAAGAGCGTCGTGTAGGGAAAGAGTGT" # Adapter 2 sequence to remove 
units: units.tsv # Units file, more on this later
outfol: /maps/projects/lundbeck/scratch/sslib_devel/mapped/ # Folder to write mapped files to 
 ```

### Units file 

You need a units file that tells snakemake where your fastqs are. This is a tab seperated file with a header (R1, R2, sm, lb) and columns:
1. R1: full path to R1 fastq
2. R2: full path to R2 fastq
3. sm: SM tag for read groups
4. lb: LB tag for read groups

If you are an internal user and your data follows the Lundbeck nomenclature, you can use the script `helpers/make_mapping_units.sh $fastq_folder`, where `$fastq_folder` is the full path to a folder of fastqs you want to map. This sets the SM tag to the first 6 characters of the first field (hyphen seperated) and the LB tag to the third field. 

### How to run 

Running is simple: 
```bash
snakemake -s map.smk --configfile config.yml -p -j12 --rerun-incomplete 
```
Remember to change the `-j` parameter to the amount of threads you want to use. 

If you want to run on slurm:

You need a `cluster_config.yml`, one is provided as an example.

Run this job on the headnode (in a screen). By default logs are created in `$(pwd)/log`, change this in the command if you want this to be different

```bash
snakemake -snakefile map.smk --configfile config.yml --use-conda -j 90 --cluster-config cluster_config.yml --cluster "sbatch --export=ALL -t {cluster.time} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --mem {cluster.memory} --partition {cluster.partition} --job-name {rulename}.{jobid} --output=$(pwd)/log/slurm-%j.out" --conda-frontend mamba --latency-wait 60  --rerun-incomplete --rerun-triggers mtime --conda-frontend mamba
```



how to stats 
complexity stats
stats explaanaaiotn

