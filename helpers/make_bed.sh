bam=$1
bed=$(mktemp)
samtools view $bam -H|grep ^@SQ|head -n22|cut -f2,3|sed 's|SN:||g'|sed 's|LN:||g'|awk '{ print $1"\t"0"\t"$2 }' > $bed
echo $bed 
