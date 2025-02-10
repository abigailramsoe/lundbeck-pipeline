set -e
units=$1
results_fol=$2
complexity=$3


# Headers
echo -n ID,Path,
echo -n MappedReads,Duplicates,Read1,Read2,ProperlyPaired,Singletons,
echo -n Depth_Autosomes,Depth_MT,Depth_X,Depth_Y,
echo -n SKOGLUND_Nseqs,SKOGLUND_NchrY+chrX,SKOGLUND_NchrY,SKOGLUND_Ry,SKOGLUND_SE,SKOGLUND_95%CI,SKOGLUND_SEX,
echo -n "C->T_5'bp1","C->T_5'bp2",
echo -n "C->T_3'bp1","C->T_3'bp2",
echo -n mtHaplo_Assignment,mtHaplo_Probability,
echo -n "ContamixApprox-1Xdif05_MAPauthentic,ContamixApprox-1Xdif05_LowerBound,ContamixApprox-1Xdif05_UpperBound,"
echo -n "ContamixPrecise-5Xdif07_MAPauthentic,ContamixPrecise-5Xdif07_LowerBound,ContamixPrecise-5Xdif07_UpperBound,"
echo -n ANGSD-XContamination_Method1,ANGSD-XContamination_Method2,ANGSD-XContamination_nSNPsites,ANGSD-XContamination_withFlanking

if [ "$complexity" != "False" ]; then
  echo -n ,Reads_Total,
  echo -n Reads_AfterTrim,
  echo -n MappingReads_NoClusterDups,MappingReads_ClusterDups,MappingReads_NoDups,Readlen,
  echo -n AfterTrim,
  echo -n Mapping,
  echo -n Clonality,
  echo -n ClusterDups,
  echo -n Endogenous,
  echo -n EndogenousUnique,
  echo -n Efficiency,
  for depth in 0.1 0.5 0.7 1 2 4 8 12; do
      echo -n ReadsToReach_${depth}x,
      echo -n ClonalityAt_${depth}x,
  done
fi
echo ""

# Now loop through each bam in pool
while read line
do
  id=$(echo "$line"|cut -f1)
  path=$(echo "$line"|cut -f2)
  if [ $id == "sampleId" ]; then continue; fi
  echo -n $id,$path,

  flagstat=$results_fol/flagstat/$id.txt
  mapped=$(cat $flagstat|grep mapped|head -n1|cut -f1 -d" ")
  dups=$(cat $flagstat|grep dup|head -n1|cut -f1 -d" ")
  read1=$(cat $flagstat|grep read1|cut -f1 -d" ")
  read2=$(cat $flagstat|grep read2|cut -f1 -d" ")
  proper=$(cat $flagstat|grep proper|cut -f1 -d" ")
  single=$(cat $flagstat|grep single|cut -f1 -d" ")
  echo -n $mapped,$dups,$read1,$read2,$proper,$single,

  depth=$results_fol/depth/$id.depth
  auto=$(cat $depth|tail -n1|cut -f2)
  mt=$(cat $depth|tail -n1|cut -f3)
  x=$(cat $depth|tail -n1|cut -f4)
  y=$(cat $depth|tail -n1|cut -f5)
  echo -n $auto,$mt,$x,$y,

  skoglund_file=$results_fol/sex/$id.sex
  skoglund_Nseqs=$(tail -n1 $skoglund_file|awk '{ print($1) }')
  skoglund_NchrYpluschrX=$(tail -n1 $skoglund_file|awk '{ print($2) }')
  skoglund_NchrY=$(tail -n1 $skoglund_file|awk '{ print($3) }')
  skoglund_Ry=$(tail -n1 $skoglund_file|awk '{ print($4) }')
  skoglund_SE=$(tail -n1 $skoglund_file|awk '{ print($5) }')
  skoglund_95CI=$(tail -n1 $skoglund_file|awk '{ print($6) }')
  skoglund_Assignment=$(tail -n1 $skoglund_file|cut -f7-|cut -c2-)
  echo -n $skoglund_Nseqs,$skoglund_NchrYpluschrX,$skoglund_NchrY,$skoglund_Ry,$skoglund_SE,$skoglund_95CI,$skoglund_Assignment,

  prof=$results_fol/damage/$id.prof
  ct1=$(cat $prof | cut -f6 | sed '2q;d')
  ct2=$(cat $prof | cut -f6 | sed '3q;d')
  ct3=$(cat $prof | cut -f6 | sed '33q;d')
  ct4=$(cat $prof | cut -f6 | sed '34q;d')
  echo -n $ct1,$ct2,$ct3,$ct4,

  haplo_fil=$results_fol/haplo/$id.haplo
  haplo=$(tail -n1 $haplo_fil | cut -f3|tr -d '"')
  haplo_prob=$(tail -n1 $haplo_fil | cut -f5|tr -d '"')
  echo -n $haplo,$haplo_prob,

  #isnt it beautiful and not a hack
  difs=(0.5 0.7)
  covs=(1 5)
  for i in "${!difs[@]}"
  do
    cov=${covs[i]}
    dif=${difs[i]}
    contamix=$results_fol/contamix/${id}_${cov}x_dif${dif}.mt.summary.txt
    map=$(cat $contamix|grep MAP|cut -f3 -d" ")
    lo=$(cat $contamix|grep MAP -B 1|grep -v MAP|cut -f1 -d" ")
    hi=$(cat $contamix|grep MAP -B 1|grep -v MAP|cut -f1 -d" ")
    echo -n $map,$lo,$hi,
  done

  angsdres=$results_fol/angsdX/$id.res
  method1=$(cat $angsdres | grep new | grep Method1| cut -f4 -d" " | cut -f2 -d:)
  method2=$(cat $angsdres | grep new | grep Method2| cut -f4 -d" " | cut -f2 -d:)
  nsnp=$(cat $angsdres | grep nSNP | cut -f1 -d, | cut -f7 -d" ")
  flanking=$(cat $angsdres | grep nSNP | cut -f2 -d, | cut -f4 -d" ")
  echo -n $method1,$method2,$nsnp,$flanking

  if [ "$complexity" != "False" ]; then
    totread=`cat $results_fol/complexity/totreads.txt|awk -v i=$id '{ if ($1==i) print($2) }'`
    discm1=`cat $results_fol/complexity/discm1.txt|awk -v i=$id '{ if ($1==i) print($2) }'`
    trim=$((${totread}-${discm1}))
    echo -n ,$totread,$trim,

    dupstat=$results_fol/superduper/$id.dupstat.txt
    pure_reads=$(grep PRC ${dupstat}| cut -f2)
    if [ -z "$pure_reads" ]
    then
      pure_reads=0
      rpf=0
      cld=0
      len=0
      clusterdup_reads=0
    else
      rpf=$(cat $dupstat|grep RPF|cut -f2)
      cld=$(cat $dupstat|grep CLD|cut -f2)
      len=$(cat $dupstat|grep CMA|cut -f2)
      clusterdup_reads=$(cat $dupstat|grep CLD|cut -f2)
    fi
    noclusterdup_reads=$(($rpf-$cld))
    echo -n $noclusterdup_reads,$clusterdup_reads,$pure_reads,$len,

    aftertrim=`bc -l <<<"${trim}/${totread}"|cut -c1-6`
    map=`bc -l <<<"${noclusterdup_reads}/${totread}"|cut -c1-6`
    clon=`bc -l <<<"1-(${pure_reads}/${noclusterdup_reads})" |cut -c1-6`
    cdups=`bc -l <<<"${clusterdup_reads}/${totread}"|cut -c1-6`
    endo=`bc -l <<<"${noclusterdup_reads}/${trim}"|cut -c1-6`
    endouniq=`bc -l <<<"${pure_reads}/${trim}"|cut -c1-6`
    effic=`bc -l <<<"${pure_reads}/${totread}"|cut -c1-6`
    echo -n $aftertrim,$map,$clon,$cdups,$endo,$endouniq,$effic

    for depth in 0.1 0.5 0.7 1 2 4 8 12; do
      txt=${results_fol}/complexity/${id}_DEPTH-${depth}.complexity.out
      reads=`cut -f3 -d" " $txt`
      complex=`cut -f4  -d" " $txt`
      echo -n ,$reads,$complex
    done
  fi



  echo ""

done < $units
