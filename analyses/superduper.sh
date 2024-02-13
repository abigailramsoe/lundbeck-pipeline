set -e

# Initialize variables
bam_file=""
output_base=""
superduper=""
threads=""
fasta_file=""
pixeldist=""
qscore=""
# Function to display usage
usage() {
  echo "Usage: $0 -b <bam_file> -f <fasta_file> -e <superduper executable> -o <output_base> -t <threads (1)> -q <quality score (30)> -p <pixel distance (12000)>"
    exit 1
}

# Parse command-line options
while getopts 'b:e:o:t:q:p:f:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        e) superduper="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        t) threads="${OPTARG}" ;;
        q) qscore="${OPTARG}" ;;
        p) pixeldist="${OPTARG}" ;;
        f) fasta_file="${OPTARG}" ;;
        *) usage ;;
    esac
done


# Check if all parameters are provided
if [[ -z $bam_file || -z $output_base || -z $superduper || -z $fasta_file ]]; then
  echo "Error: Bam, superduper executable, output and fasta file are required (-b,-e,-o,-f)"
  usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_file.sh "$superduper"
bash helpers/check_file.sh "$fasta_file"
threads=$(bash helpers/check_threads.sh $threads)

bash helpers/check_directory.sh $(dirname "$output_base")

if [ -z "$pixeldist" ]; then pixeldist=12000; fi
if [ -z "$qscore" ]; then qscore=30; fi

$superduper -w -b ${bam_file} -T $fasta_file -q $qscore -o $output_base -p $pixeldist -@$threads -m
