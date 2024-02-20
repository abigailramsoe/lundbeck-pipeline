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
  echo "Usage: $0 -t <totreads_file> -a <superduper table file> -u <superduper dupstat file> -o <output_base> -d <depth>"
    exit 1
}

# Parse command-line options
while getopts 't:a:u:o:d:' flag; do
    case "${flag}" in
        t) totreads_file="${OPTARG}" ;;
        a) table="${OPTARG}" ;;
        u) dup="${OPTARG}" ;;
        o) output_base="${OPTARG}" ;;
        d) depth="${OPTARG}" ;;
        *) usage ;;
    esac
done

# Check if all parameters are provided
if [[ -z $totreads_file || -z $table || -z $dup || -z $output_base || -z $depth ]]; then
    echo "Error: All parameters are required."
    usage
fi
dir=$(dirname $0)/../

# Validate files and directories
bash $dir/helpers/check_file.sh "$totreads_file"
bash $dir/helpers/check_file.sh "$table" # might need to make this less strict, sometimes there is no table
bash $dir/helpers/check_file.sh "$dup"

bash $dir/helpers/check_directory.sh $(dirname $output_base)

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

id=$(basename $output_base)
totreads=$(cat $totreads_file|awk -v i=$id '{ if ($1==i) print($2) }')


out=$output_base"_DEPTH-${depth}.complexity.out"
if [ ! -s $table ]; then
 echo -e "$totreads $ncd NaN NaN" > $out
else
 Rscript analyses/complex.R $table $depth $read_len $totreads $ncd > $out
fi
