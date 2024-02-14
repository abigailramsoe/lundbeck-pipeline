you need java >1.8
r 3.6.1
and some random other things 
angsd
superduper
bam2prof

conda env 
Install from source to use the development version
By cloning in a dedicated conda environment

git clone git@github.com:genomewalker/aMGSIM.git
cd aMGSIM
conda env create -f environment.yml
conda activate aMGSIM
pip install -e .

only works on PE
stats only for humans 

So you have some fasts you want to map? You've come to the right place. First, make a config file.  There is an exampe in config.yml

```yaml
ref: /maps/projects/lundbeck/data/hg38/genome.fa  # Path to fasta reference genome 
seed: 32  # Seed (-l) for bwa
a1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # Adapter 1 sequence to remove
a2: "GGAAGAGCGTCGTGTAGGGAAAGAGTGT" # Adapter 2 sequence to remove 
units: units.tsv # Units file, more on this later
outfol: /maps/projects/lundbeck/scratch/sslib_devel/mapped/ # Folder to write mapped files to 
 ```

units file 
how to run 
how to slurm 

how to stats 
complexity stats
stats explaanaaiotn

