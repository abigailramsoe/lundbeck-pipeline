set -e

# Initialize variables
bam_file=""
output_base=""
threads=""

# Function to display usage
usage() {
    echo "Usage: $0 -b <bam_file> -o <output_base> -t <threads>"
    exit 1
}

# Parse command-line options
while getopts 'b:f:o:t:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        t) region="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $output_base ]]; then
    echo "Error: All parameters are required."
    usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_directory.sh $(dirname "$output_base")
threads=$(bash helpers/check_threads.sh $threads)

bash helpers/index.sh $bam_file

bed=$(bash helpers/make_bed.sh $bam_file)
mt_region=$(samtools view $bam_file -H|grep ^@SQ|grep M|head -n1|cut -f2|cut -f2 -d":")
x_region=$(samtools view $bam_file -H|grep ^@SQ|grep X|head -n1|cut -f2|cut -f2 -d":")
y_region=$(samtools view $bam_file -H|grep ^@SQ|grep Y|head -n1|cut -f2|cut -f2 -d":")

auto=$(samtools depth $bam_file -a -Q 30 -q 20 -b $bed -@$threads|datamash mean 3)
mt=$(samtools depth $bam_file -a -Q 30 -q 20 -r $mt_region -@$threads|datamash mean 3)
x=$(samtools depth $bam_file -a -Q 30 -q 20 -r $x_region -@$threads|datamash mean 3)
y=$(samtools depth $bam_file -a -Q 30 -q 20 -r $y_region -@$threads|datamash mean 3)

echo -e "Input\tAutosomalCov\tMTCov\tXCov\tYCov"
echo -e "$bam_file\t$auto\t$mt\t$x\t$y"
