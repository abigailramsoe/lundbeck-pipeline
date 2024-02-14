set -e

# Initialize variables
bam_file=""
totreads_file=""
table=""
dup=""
output_base=""
depth=""

# Function to display usage
usage() {
  echo "Usage: $0 -b <bam_file> -t <totreads_file> -a <superduper table file> -u <superduper dupstat file> -o <output_base> -d <depth>"
    exit 1
}

# Parse command-line options
while getopts 'b:t:a:u:o:d:' flag; do
    case "${flag}" in
        b) bam_file="${OPTARG}" ;;
        t) totreads_file="${OPTARG}" ;;
        a) table="${OPTARG}" ;;
        u) dup="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        d) depth="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $bam_file || -z $totreads_file || -z $table || -z $dup || -z $output_base || -z $depth ]]; then
    echo "Error: All parameters are required."
    usage
fi

# Validate files and directories
bash helpers/check_file.sh "$bam_file"
bash helpers/check_file.sh "$totreads_file"
bash helpers/check_file.sh "$table" # might need to make this less strict, sometimes there is no table 
bash helpers/check_file.sh "$dup"

bash helpers/check_directory.sh $(dirname $output_base)

if [[ $(bc <<< "$depth <= 0") == 1 ]]
then
  echo "Error: depth must be a float greater than 0 (current value $depth)"
  exit 1
fi

pure=$(cat $dup |grep PRC| cut -f2)
#ncd is rpf - cld
rpf=$(cat $dup|grep RPF|cut -f2)
cld=$(cat $dup|grep CLD|cut -f2)
ncd=$(($rpf-$cld))
read_len=$(cat $dup|grep CMA|cut -f2)

bambase=$(basename $bam_file|cut -f1 -d.)
totreads=$(cat $totreads_file | grep $bambase | cut -f2)


out=$output_base"_DEPTH-${depth}.complexity.out"

if [ ! -s $table ]; then
 echo -e "$totreads $ncd NaN NaN" > $out
else
 Rscript analyses/complex.R $table $depth $read_len $totreads $ncd > $out
fi
