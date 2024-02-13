set -e

# Initialize variables
bam_file=""
output_base=""
bam2prof=""
length=""
# Function to display usage
usage() {
  echo "Usage: $0 -b <bam_file> -e <bam2prof executable> -o <output_base> -l <length>"
    exit 1
}

# Parse command-line options
while getopts 'b:e:o:l' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        e) bam2prof="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        l) length="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $output_base || -z $bam2prof ]]; then
  echo "Error: Bam, bam2prof executable and output are required (-b,-e,-o)"
  usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_file.sh "$bam2prof"

bash helpers/check_directory.sh $(dirname "$output_base")

if  [ -z "$length" ]
then
  length=30
fi

$bam2prof -length $length ${bam_file} 2>/dev/null > $output_base
echo "Done, look for $output_base"
