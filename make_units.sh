
fqfol=$1

echo -e "fastq\tsm\tlb\tid"
for fastq in $(ls $fqfol/*.fastq.gz)
do
  sm=$(basename $fastq|cut -f1 -d-|cut -c1-6)
  lb=$(basename $fastq|cut -f3 -d-)
  id=$(basename $fastq|cut -f1 -d_)
  echo -e "$(realpath $fastq)\t$sm\t$lb\t$id"
done
 
