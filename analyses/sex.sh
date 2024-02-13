set -e

# Initialize variables
bam_file=""
fasta_file=""
output_file=""
qscore=""

# Function to display usage
usage() {
    echo "Usage: $0 -b <bam_file> -f <fasta_file> -o <output_base> -q <quality score (30)>"
    exit 1
}

# Parse command-line options
while getopts 'b:f:r:o:q:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        f) fasta_file="${OPTARG}" ;;
        o) output_file="${OPTARG}" ;;
        q) qscore="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $fasta_file || -z $output_file ]]; then
    echo "Error: Bam, fasta and output are required (-b,-f,-o)"
    usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_file.sh "$fasta_file"
bash helpers/check_directory.sh $(dirname "$output_file")

samtools view -q 30 $bam_file -T $fasta_file | helpers/skoglund.ry.py > $output_file
