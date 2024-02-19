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

## Statistics 

This is the pipeline for the statistics part. It can also estimate library complexity and do sequencing projections (`lib_complexity` option in config file), but only if the data was also mapped by this pipeline. 

### Config file 
You need a config file, there is an example in `config-stat.yml`

```yaml
ref: /maps/projects/lundbeck/data/hg38/genome.fa # Path to fasta reference genome 
units: sslib_devel_units # Units file, more on this later
outfol: /projects/lundbeck/scratch/sslib_devel/results/ # Folder to store the results 
angsd_dir: /projects/lundbeck/apps/angsd/ # Path to ANGSD directory, https://github.com/ANGSD/angsd
bam2prof_exec: /projects/lundbeck/apps/bam2prof/src/bam2prof # Path to bam2prof executable https://github.com/grenaud/bam2prof
superduper_exec: /projects/lundbeck/apps/decluster/decluster # Path to decluster executable https://github.com/ANGSD/decluster
lib_complexity: yes # Do you want to estimate library complexity? Only available if mapped by this pipeline 
```

### Units file 

The units file is a tab seperated file with the header (sampleId, path) and columns:
1. sampleId: unique sample ID for each bam
2. path: full path for each file

There is an example in `sslib_devel_units`, in the future there should be an automated way of getting this from the mapping pipeline 

### How to run 

```bash
snakemake -s stats.smk --configfile config-stat.yml -p -j12 --rerun-incomplete 
```

Or on slurm, the same way as the mapping pieline. There is a cluster-config for the statistics part in `cluster_config-stats.jml`

### Results 

A csv with results will be generated in `$outfol/results.csv`, where `$outfol` is the output folder you specified in the config file. 

Here is an explanation of each field:

```
1. ID from units file
2. Full path to bam, from units file
3. Number of single-end mapping reads without cluster duplicates (Decluster)
4. Number of single-end mapping cluster duplicates (Decluster)
5. Number of single-end mapping reads without any duplicates (Decluster)
6. Average read length of single-end reads (Decluster)
7. Average depth of coverage for autosomes (samtools depth)
8. Average depth of coverage of chrM (samtools depth)
9. Average depth of coverage of chrX (samtools depth)
10. Average depth of coverage of chrY (samtools depth)
11. Number of sequences (Skoglund)
12. Number of chrY + chrX sequences (Skoglund)
13. Number of chrY sequences (Skoglund)
14. Ry coefficient (Skoglund)
15. Standard error for Ry (Skoglund)
16. 95% confidence intervals (Skoglund)
17. Sex determination (Skoglund)
18. C->T proportion at first basepair of 5' end (bam2prof)
19. C->T proportion at second basepair of 5' end (bam2prof)
20. C->T proportion at first basepair of 3' end (bam2prof)
21. C->T proportion at second basepair of 3' end (bam2prof)
22. chrM haplogroup assignment (haplogrep v2)
23. chrM haplogroup probability (haplogrep v2)
24. Contamix MAP authentic value for approximate settings (contamix)
25. Contamix lower bound for approximate settings (contamix)
26. Contamix upper bound for approximate settings (contamix)
27. Contamix MAP authentic value for precise settings (contamix)
28. Contamix lower bound for precise settings (contamix)
29. Contamix upper bound for precise settings (contamix)
30. ANGSD chrX contamination estimate, method 1 (ANGSD)
31. ANGSD chrX contamination estimate, method 2 (ANGSD)
32. Number of chrX SNP sites (ANGSD)
33. Number of sites flanking chrX SNP sites (ANGSD)
```

If library projections are included:

```
34. Total number of reads sequenced (AdapterRemoval)
35. Total number of reads (Reads_Total) minus number of discarded mate 1 reads (AdapterRemoval)
36. Proportion of remaining reads after trimming (Reads_AfterTrim/Reads_Total)
37. Proportion of mapping reads (MappingReads_NoClusterDups/Reads_Total)
38. Clonality (1-(MappingReads_NoDups/MappingReads_NoClusterDups)
39. Cluster duplicates (MappingReads_ClusterDups/Reads_Total)
40. Endogenous (MappingReads_NoClusterDups/Reads_AfterTrim)
41. Endogenous with no duplicates (MappingReads_NoDups/Reads_AfterTrim)
42. Efficiency (MappingReads_NoDups/Reads_Total)
43. ReadsToReach_01x (Inf if unreachable or if clonality will be over 75%)
44. ClonalityAt_01x
45. ReadsToReach_05x
46. ClonalityAt_05x
47. ReadsToReach_07x
48. ClonalityAt_07x
49. ReadsToReach_1x
50. ClonalityAt_1x
51. ReadsToReach_2x
52. ClonalityAt_2x
53. ReadsToReach_4x
54. ClonalityAt_4x
55. ReadsToReach_8x
56. ClonalityAt_8x
57. ReadsToReach_12x
58. ClonalityAt_12x
```
