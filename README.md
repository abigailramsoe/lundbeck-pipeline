![image](https://github.com/abigailramsoe/lundbeck-pipeline/assets/28560412/f96f3a70-7c3d-4398-97fd-e136a129378e)

This pipeline implements the mapping and statistics workflow of the Lundbeck Foundation GeoGenetics Center's bioinformatics pipeline. Mapping assumes paired end sequencing, while the statistics pipeline assumes single end (or a mixture of paired and single end) reads.

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
* [bam2prof](https://github.com/grenaud/bam2prof)
* [Decluster/Superduper](https://github.com/ANGSD/decluster) (only if you want library complexity projections)


# Usage

## Mapping

### Config file

First, make a config file.  There is an exampe in `configs/config.yml`

```yaml
ref: /maps/projects/lundbeck/data/hg38/genome.fa  # Path to fasta reference genome
seed: 32  # Seed (-l) for bwa
a1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # Adapter 1 sequence to remove
a2: "GGAAGAGCGTCGTGTAGGGAAAGAGTGT" # Adapter 2 sequence to remove
units: units/units.tsv # Units file, more on this later
outfol: /maps/projects/lundbeck/scratch/sslib_devel/mapped/ # Folder to write mapped files to
 ```

Remember to change a1 and a2 to suit your needs.
Here are some examples. If you are an internal users and have question about what sequence to use, ask the sequencing center (seqcenter@sund.ku.dk) or the bioinformatics team (cggbinf@list.ku.dk)

|| A1 (P7) | A2 (P5) |
| --- | --- | --- |
|Double-stranded libraries and Santa Cruz Reaction|AGATCGGAAGAGCACACGTCTGAACTCCAGTCA|AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT|
Single-stranded (ss2.0) libraries and SCR with P5short|AGATCGGAAGAGCACACGTCTGAACTCCAGTCA|GGAAGAGCGTCGTGTAGGGAAAGAGTGT|



### Units file

You need a units file that tells snakemake where your fastqs are. This is a tab seperated file with a header (R1, R2, sm, lb) and columns:
1. R1: full path to R1 fastq
2. R2: full path to R2 fastq
3. sm: SM tag for read groups
4. lb: LB tag for read groups

If you are an internal user and your data follows the Lundbeck nomenclature, you can use the script `helpers/make_mapping_units.sh $fastq_folder > units/units.tsv`, where `$fastq_folder` is the full path to a folder of fastqs you want to map. This sets the SM tag to the first 6 characters of the first field (hyphen seperated) and the LB tag to the third field.

### How to run

Running is simple:
```bash
snakemake -s map.smk --configfile configs/config.yml -p -j12 --rerun-incomplete
```
Remember to change the `-j` parameter to the amount of threads you want to use.

If you want to run on slurm:

You need a `cluster_config.yml`, one is provided as an example.

Run this job on the headnode (in a screen). By default logs are created in `$(pwd)/log`, change this in the command if you want this to be different

```bash
snakemake -s map.smk --configfile configs/config.yml --use-conda -j 90 --cluster-config configs/cluster_config.yml --cluster "sbatch --export=ALL -t {cluster.time} --ntasks-per-node {cluster.ntasks_per_node} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus_per_task} --mem {cluster.memory} --partition {cluster.partition} --job-name {rulename}.{jobid} --output=$(pwd)/log/slurm-%j.out" --conda-frontend mamba --latency-wait 60  --rerun-incomplete --rerun-triggers mtime --conda-frontend mamba
```

## Statistics

This is the pipeline for the statistics part. It can also estimate library complexity and do sequencing projections (`lib_complexity` option in config file), but only if the data was also mapped by this pipeline.

### Config file
You need a config file, there is an example in `config-stats.yml`

```yaml
ref: /maps/projects/lundbeck/data/hg38/genome.fa # Path to fasta reference genome
units: units/stats_units.tsv # Units file, more on this later
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

There is an example in `units/stats_units.tsv`, in the future there should be an automated way of getting this from the mapping pipeline

### How to run

```bash
snakemake -s stats.smk --configfile configs/config-stat.yml -p -j12 --rerun-incomplete
```

Or on slurm, the same way as the mapping pipeline. There is a cluster-config for the statistics part in `configs/cluster_config-stats.jml`

### Results

A csv with results will be generated in `$outfol/results.csv`, where `$outfol` is the output folder you specified in the config file.

Here is an explanation of each field:

```
1. ID from units file
2. Full path to bam, from units file
3. Number of mapped reads (samtools flagstat)
4. Number of duplicates (samtools flagstat)
5. Number of mapped read 1 reads (samtools flagstat)
6. Number of mapped read 2 reads (samtools flagstat)
7. Number of properly paired reads (samtools flagstat)
8. Number of singletons (samtools flagstat)
9. Average depth of coverage for autosomes (samtools depth)
10. Average depth of coverage of chrM (samtools depth)
11. Average depth of coverage of chrX (samtools depth)
12. Average depth of coverage of chrY (samtools depth)
13. Number of sequences (Skoglund)
14. Number of chrY + chrX sequences (Skoglund)
15. Number of chrY sequences (Skoglund)
16. Ry coefficient (Skoglund)
17. Standard error for Ry (Skoglund)
18. 95% confidence intervals (Skoglund)
19. Sex determination (Skoglund)
20. C->T proportion at first basepair of 5' end (bam2prof)
21. C->T proportion at second basepair of 5' end (bam2prof)
22. C->T proportion at first basepair of 3' end (bam2prof)
23. C->T proportion at second basepair of 3' end (bam2prof)
24. chrM haplogroup assignment (haplogrep v2)
25. chrM haplogroup probability (haplogrep v2)
26. Contamix MAP authentic value for approximate settings (contamix)
27. Contamix lower bound for approximate settings (contamix)
28. Contamix upper bound for approximate settings (contamix)
29. Contamix MAP authentic value for precise settings (contamix)
30. Contamix lower bound for precise settings (contamix)
31. Contamix upper bound for precise settings (contamix)
32. ANGSD chrX contamination estimate, method 1 (ANGSD)
33. ANGSD chrX contamination estimate, method 2 (ANGSD)
34. Number of chrX SNP sites (ANGSD)
35. Number of sites flanking chrX SNP sites (ANGSD)
```

If library projections are included:

```
36. Total number of reads sequenced (AdapterRemoval)
37. Total number of reads (Reads_Total) minus number of discarded mate 1 reads (AdapterRemoval)
38. Number of single-end mapping reads without cluster duplicates (Decluster)
39. Number of single-end mapping cluster duplicates (Decluster)
40. Number of single-end mapping reads without any duplicates (Decluster)
41. Average read length of single-end reads (Decluster)
42. Proportion of remaining reads after trimming (Reads_AfterTrim/Reads_Total)
43. Proportion of mapping reads (MappingReads_NoClusterDups/Reads_Total)
44. Clonality (1-(MappingReads_NoDups/MappingReads_NoClusterDups)
45. Cluster duplicates (MappingReads_ClusterDups/Reads_Total)
46. Endogenous (MappingReads_NoClusterDups/Reads_AfterTrim)
47. Endogenous with no duplicates (MappingReads_NoDups/Reads_AfterTrim)
48. Efficiency (MappingReads_NoDups/Reads_Total)
49. ReadsToReach_01x (Inf if unreachable or if clonality will be over 75%)
50. ClonalityAt_01x
51. ReadsToReach_05x
52. ClonalityAt_05x
53. ReadsToReach_07x
54. ClonalityAt_07x
55. ReadsToReach_1x
56. ClonalityAt_1x
57. ReadsToReach_2x
58. ClonalityAt_2x
59. ReadsToReach_4x
60. ClonalityAt_4x
61. ReadsToReach_8x
62. ClonalityAt_8x
63. ReadsToReach_12x
64. ClonalityAt_12x
```
