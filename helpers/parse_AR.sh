set -e

settingsfol=$1
outfol=$2

arsettings=$outfol/ar.full.txt

first=true
for i in $(find $settingsfol/ -name "*settings")
do
  if [ $first == true ]; then
    header=$(grep "statistics" -A13 ${i}|sed 1d|cut -f1 -d:|tr " " "_"|tr "\n" ";")
    echo ";$header" > $arsettings
    first=false
  fi
  id=`basename $i .settings`
  val=`grep "statistics" -A13 ${i}|sed 1d|cut -f2 -d:|sed -r 's/\s+//g'|tr "\n" ";"`
  echo -e ${id}";"${val} >> $arsettings
  length=$outfol/$id.length.tsv
  grep "Length dist" -A 200 $i|sed 1d > $length
done
echo "Parsed and written AR output: $arsettings"
echo "Parsed and written length distribution, eg: $length"


totreads=$outfol/totreads.txt
disc=$outfol/discm1.txt
> $totreads
> $disc
for i in `cut -f1 -d_ $arsettings|sed 1d|sort|uniq`
do

    t=`grep ${i} ${arsettings}|cut -f2 -d";"|datamash sum 1`
    echo -e "${i}\t${t}" >> $totreads

    d=`grep ${i} ${arsettings}|cut -f5 -d";"|datamash sum 1`
    echo -e "${i}\t${d}" >> $disc

done
echo "Parsed and written totreads: $totreads"
echo "Parsed and written discarded m1 reads: $disc"
