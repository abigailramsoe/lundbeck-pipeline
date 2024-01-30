
fqfol=$1

echo -e "R1\tR2\tsm\tlb\tid\tidlane"
for r1 in $(ls $fqfol/*_R1_001.fastq.gz)
do
  sm=$(basename $r1|cut -f1 -d-|cut -c1-6)
  lb=$(basename $r1|cut -f3 -d-)
  id=$(basename $r1|cut -f1 -d_)
  id_wlane=$(basename $r1|cut -f1-3 -d_)
  r2=$(echo $r1|sed 's|R1_001|R2_001|g')
  echo -e "$(realpath $r1)\t$(realpath $r2)\t$sm\t$lb\t$id\t$id_wlane"
done
