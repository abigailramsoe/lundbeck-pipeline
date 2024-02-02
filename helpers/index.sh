set -e
bam=$1
threads=$2

suffix=$(basename $bam|rev|cut -f1 -d.|rev)
idx_suffix=""

if [ $suffix == "cram" ]; then idx_suffix="crai"; fi
if [ $suffix == "bam" ]; then idx_suffix="bai"; fi

if [ -z "$idx_suffix" ]
then
  echo "Suffix is not cram or bam, exiting"
  exit
fi

threads=$(bash helpers/check_threads.sh $threads)


idx=$bam.$idx_suffix
if [ ! -f "$idx" ]; then
  samtools index $bam -@$threads
fi
