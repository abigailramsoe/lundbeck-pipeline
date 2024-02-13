#!/usr/bin/bash

bam=$1
cov=$2
diff=0.$3
OFOL=$4
ref=$5

module unload samtools
module load samtools/1.12
mkdir -p $OFOL

if [ "$ref" != "old" ] && [ "$ref" != "new" ]
then
  echo "Please specify old or new ref"
  exit
fi
fasta=/projects/lundbeck/data/hg38/genome.fa
region=chrM
if [ "$ref" == "old" ]
then
  fasta=/projects/lundbeck/data/hs37d5/hs37d5.fa
  region=MT
fi


#get the mtDNA mapped reads
#samtools index ${1}
base=$OFOL/`basename $bam|cut -f1 -d.`
samtools view -b $bam $region -F 1 > ${base}_mt.bam
echo "########### mtDNA bam OK ###########"

x=$(samtools view ${base}_mt.bam -c )
if [ $x == 0 ];
then
  base=$OFOL/`basename ${base}_mt`

echo 0 mt reads > ${base}.summary.txt
exit
fi
#get the variable sites
base=$OFOL/`basename ${base}_mt`


samtools mpileup -d 2000 -m3 -C50 -q 30 -EQ 20 -uf $fasta ${base}.bam | bcftools call -m --ploidy 1 > ${base}.vcf
echo "########### vcf OK ###########"

#get the consensus sequence with Victor's script (possible with ANGSD)
/projects/lundbeck/apps/contamix_helpers/CnsMaj3_1.pl -i ${base}.vcf -o ${base}.fa -l 16569 -cov $cov -diff $diff -idiff 0.5 -h ${base} -callindels no > ${base}.cns
echo "########### consensus file OK ###########"


#add the consensus to the 311 mtDNA sequences and align all with mafft
cat ${base}.fa /projects/lundbeck/apps/contamix_helpers/311humans.fasta > ${base}_312.fasta
mafft ${base}_312.fasta > ${base}_312.aligned.fasta
echo "########### alignment with 311 mtDNAs OK ###########"

#index the consensus
bwa index -a bwtsw ${base}.fa
module unload samtools
module load samtools
samtools faidx ${base}.fa
rm ${base}.dict
module load jdk
java -jar /projects/lundbeck/apps/contamix_helpers/picard.jar CreateSequenceDictionary R=${base}.fa O=${base}.dict
echo "########### consensus indexed OK ###########"


#get fastq files from the bam
#### gone wrong
#java -jar /home/fvr124//bin/picard.jar SamToFastq INPUT=${base}.bam FASTQ=${base}.fq
samtools view ${base}.bam |cut -f1,10,11 | awk '{print "@"$1"\n"$2"\n+\n"$3}' > ${base}.fq

echo "########### fastq from bam OK ###########"


#map the above fastq file to the consensus
cons=${base}.fa
bwa aln -l 1000 -t 10 $cons ${base}.fq > ${base}_ra.sai
bwa samse -r '@RG\tID:${base}\tLB:${base}_L1\tPL:ILLUMINA\tSM:${base}' $cons ${base}_ra.sai ${base}.fq | samtools view -bh -q 30 | samtools sort -O BAM -o ${base}_ra.sort.bam

#samtools view -h ${base}_ra.sort.bam | awk '{if($0~/X1:i:0/||$0~/^@/){print $0}}' | samtools view -bSh - > ${base}_ra.bam
java -jar /projects/lundbeck/apps/contamix_helpers/picard.jar MarkDuplicates I=${base}_ra.sort.bam O=${base}_ra.sort.rmdup.bam TMP_DIR=temp METRICS_FILE=${base}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
samtools index ${base}_ra.sort.rmdup.bam
rm ${base}_ra.sai

java -jar /projects/lundbeck/apps/contamix_helpers/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $cons -I ${base}_ra.sort.rmdup.bam -o ${base}.intervals
java -jar /projects/lundbeck/apps/contamix_helpers/GenomeAnalysisTK.jar -T IndelRealigner -R $cons -I ${base}_ra.sort.rmdup.bam -targetIntervals ${base}.intervals -o ${base}_ra.sort.rmdup.realn.bam
samtools calmd -Erb ${base}_ra.sort.rmdup.realn.bam $cons > ${base}_ra.final.bam

samtools index ${base}_ra.final.bam



x=$(samtools view ${base}_ra.final.bam -c )
if [ $x -lt 100 ];
then

echo too few mt reads in final > ${base}.summary.txt
exit
fi
echo "########### map fastq to the consensus OK ###########"

module load gcc
module load R/3.6.1
echo ${base}_ra.final.bam
echo ${base}_312.aligned.fasta
#do the actual comtamMix part
echo /projects/lundbeck/apps/contamix_helpers/estimate.R --samFn ${base}_ra.final.bam --malnFn ${base}_312.aligned.fasta --figure ${base}_contam.pdf --trimBases 7 > ${base}.summary.txt
