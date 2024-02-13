set -e

# Initialize variables
bam_file=""
fasta_file=""
region=""
output_file=""

# Function to display usage
usage() {
    echo "Usage: $0 -b <bam_file> -f <fasta_file> -o <output_base>"
    exit 1
}

# Parse command-line options
while getopts 'b:f:r:o:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        f) fasta_file="${OPTARG}" ;;
        o) output_file="${OPTARG}" ;;
        r) region="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $fasta_file || -z $output_file ]]; then
    echo "Error: Bam, fasta and output are required (-b,-f,-po)"
    usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_file.sh "$fasta_file"
bash helpers/check_directory.sh $(dirname "$output_file")

region=$(bash helpers/find_mt.sh "$fasta_file")

vcf=$output_file.vcf.gz
bash helpers/index.sh $bam_file
bcftools mpileup $bam_file --ignore-RG -Ou --fasta-ref $fasta_file -r $region | bcftools call -m - -Oz -o $vcf
bcftools index $vcf
java -jar helpers/haplogrep-2.1.25.jar --format vcf --in $vcf --out $output_file
rm $vcf 
