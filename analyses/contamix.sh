set -e

# Initialize variables
bam_file=""
fasta_file=""
output_base=""
coverage=""
diff=""

# Function to display usage
usage() {
  echo "Usage: $0 -b <bam_file> -f <fasta_file> -o <output_base> -c <minimum coverage for consensus> -d <proportion of agree for consensus>"
  echo "Recommended: -c 1 -d 0.5 (approximate results for low coverage)"
  echo "Recommended: -c 5 -d 0.7 (precise results for high coverage)"
  exit 1
}

# Parse command-line options
while getopts 'b:f:o:c:d:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        f) fasta_file="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        c) coverage="${OPTARG}" ;;
        d) diff="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $fasta_file || -z $output_base || -z $coverage || -z $diff ]]; then
    echo "Error: All parameters are required."
    usage
fi
dir=$(dirname $0)/../

# Validate files and directories
bash $dir/helpers/check_file.sh "$bam_file"
bash $dir/helpers/check_file.sh "$fasta_file"

bash $dir/helpers/check_directory.sh $(dirname "$output_base")


if ! [[ "$coverage" =~ ^[0-9]+$ ]] ; then
   echo "Error: coverage must be an integer (current value $coverage)"
   exit 1
fi
if [[ $(bc <<< "$diff <= 0 || $diff >= 1") == 1 ]]
then
  echo "Error: diff must be a float between 0 and 1 (current value $diff)"
  exit 1
fi

region=$(bash $dir/helpers/find_mt.sh "$fasta_file")

if [ $(hostname|cut -c1-5) == "dandy" ]; then
  echo "Detected dandy system, loading jdk and the good R version (3.6.1)"
  module load gcc jdk R/3.6.1
else
  echo "Detected you are not on the UCPH Geogenetics Dandy Server"
  echo "Make sure you have Java 1.8 and R 3.6.1"
  echo "Sleeping for 10 secs so you take this seriously"
  sleep 10
fi



bash $dir/helpers/index.sh $bam_file

base=$output_base.mt
samtools view -b $bam_file $region -F 1 -T $fasta_file  > ${base}.bam

if [ $(samtools view ${base}.bam -c) == 0 ];
then
  echo "No MT reads in $bam_file" > ${base}.summary.txt
  echo "No MT reads in $bam_file"
  exit 0
fi

bcftools mpileup -d 2000 -q 30 -Q 20 -f $fasta_file ${base}.bam | bcftools call -m --ploidy 1 > ${base}.vcf
echo "########### vcf OK ###########"

#get the consensus sequence with Victor's script (possible with ANGSD)
$dir/helpers/contamix_helpers/CnsMaj3_1.pl -i ${base}.vcf -o ${base}.fa -l 16569 -cov $coverage -diff $diff -idiff 0.5 -h ${base} -callindels no > ${base}.cns
echo "########### consensus file OK ###########"


#add the consensus to the 311 mtDNA sequences and align all with mafft
cat ${base}.fa $dir/helpers/contamix_helpers/311humans.fasta > ${base}_312.fasta
mafft ${base}_312.fasta > ${base}_312.aligned.fasta
echo "########### alignment with 311 mtDNAs OK ###########"

#index the consensus
bwa index -a bwtsw ${base}.fa
samtools faidx ${base}.fa
if [ -f ${base}.dict ]; then rm ${base}.dict; fi
java -jar $dir/helpers/contamix_helpers/picard.jar CreateSequenceDictionary R=${base}.fa O=${base}.dict
echo "########### consensus indexed OK ###########"

#get fastq files from the bam
samtools view ${base}.bam |cut -f1,10,11 | awk '{print "@"$1"\n"$2"\n+\n"$3}' > ${base}.fq
echo "########### fastq from bam OK ###########"

#map the above fastq file to the consensus
cons=${base}.fa
bwa aln -l 1000 -t 10 $cons ${base}.fq > ${base}_ra.sai
bwa samse -r '@RG\tID:${base}\tLB:${base}_L1\tPL:ILLUMINA\tSM:${base}' $cons ${base}_ra.sai ${base}.fq | samtools view -bh -q 30 | samtools sort -O BAM -o ${base}_ra.sort.bam

java -jar $dir/helpers/contamix_helpers/picard.jar MarkDuplicates I=${base}_ra.sort.bam O=${base}_ra.sort.rmdup.bam TMP_DIR=temp METRICS_FILE=${base}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
samtools index ${base}_ra.sort.rmdup.bam
rm ${base}_ra.sai
java -jar $dir/helpers/contamix_helpers/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $cons -I ${base}_ra.sort.rmdup.bam -o ${base}.intervals
java -jar $dir/helpers/contamix_helpers/GenomeAnalysisTK.jar -T IndelRealigner -R $cons -I ${base}_ra.sort.rmdup.bam -targetIntervals ${base}.intervals -o ${base}_ra.sort.rmdup.realn.bam
samtools calmd -Erb ${base}_ra.sort.rmdup.realn.bam $cons > ${base}_ra.final.bam
samtools index ${base}_ra.final.bam



if [ $(samtools view ${base}_ra.final.bam -c) -lt 10 ];
then
  echo "No MT reads in final bam file ${base}_ra.final.bam" >> ${base}.summary.txt
  echo "No MT reads in final bam file ${base}_ra.final.bam"
  exit 0
fi
echo "########### map fastq to the consensus OK ###########"

#do the actual comtamMix part
Rscript $dir/helpers/contamix_helpers/estimate.R --samFn ${base}_ra.final.bam --malnFn ${base}_312.aligned.fasta --figure ${base}_contam.pdf --trimBases 7 > ${base}.summary.txt
