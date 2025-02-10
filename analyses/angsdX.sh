set -e

# Initialize variables
bam_file=""
fasta_file=""
output_base=""
angsd_dir=""

# Function to display usage
usage() {
  echo "Usage: $0 -b <bam_file> -f <fasta_file> -o <output_base> -e <path to angsd directory>"
    exit 1
}

# Parse command-line options
while getopts 'b:f:o:e:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        e) angsd_dir="${OPTARG}" ;;
        f) fasta_file="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $fasta_file || -z $output_base || -z $angsd_dir ]]; then
    echo "Error: All parameters are required."
    usage
fi

dir=$(dirname $0)/../

# Validate files and directories
bash $dir/helpers/check_file.sh "$bam_file"
bash $dir/helpers/check_file.sh "$fasta_file"
bash $dir/helpers/check_file.sh "$angsd_dir/angsd"

bash $dir/helpers/check_directory.sh $(dirname "$output_base")

x_region=$(samtools view $bam_file -H|grep ^@SQ|grep X|head -n1|cut -f2|cut -f2 -d":")

hapmap=""
if [ $x_region == "chrX" ]; then hapmap=$dir/helpers/ChrX.hg38.hapmap.gz; fi
if [ $x_region == "X" ]; then hapmap=$dir/helpers/ChrX.hs37d5.hapmap.gz; fi

if [ ! -f "$hapmap" ]
then
  echo "Could not find HAPMAP given X chromosome $x_region"
  exit 1
fi

${angsd_dir}/angsd -i ${bam_file} -r $x_region:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 20 -out ${output_base}
${angsd_dir}/misc/contamination -a ${output_base}.icnts.gz -h $hapmap &> ${output_base}






