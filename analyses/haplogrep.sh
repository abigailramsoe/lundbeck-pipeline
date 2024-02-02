set -e

# Initialize variables
bam_file=""
fasta_file=""
region=""
output_base=""

# Function to display usage
usage() {
    echo "Usage: $0 -b <bam_file> -f <fasta_file> -o <output_base>"
    exit 1
}


check_region() {
    if [[ $1 != "chrM" ]] && [[ $1 != "MT" ]]; then
        echo "Could not find what the MT reference is called automatically"
        echo "Please supply the region name with -r"
        exit 1
    fi
}


# Parse command-line options
while getopts 'b:f:r:o:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        f) fasta_file="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        r) region="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $fasta_file || -z $output_base ]]; then
    echo "Error: All parameters are required."
    usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_file.sh "$fasta_file"
bash helpers/check_directory.sh $(dirname "$output_base")

if [ -z "$region" ]
then
  region=$(grep ">" $fasta_file|grep M|head -n1|cut -c2-)
  check_region $region
else
  echo "Using user defined region $region"
fi


vcf=$output_base.vcf.gz
bash helpers/index.sh $bam_file
bcftools mpileup $bam_file --ignore-RG -Ou --fasta-ref $fasta_file -r $region | bcftools call -m - -Oz -o $vcf
bcftools index $vcf
java -jar helpers/haplogrep-2.1.25.jar --format vcf --in $vcf --out $output_base.haplo
