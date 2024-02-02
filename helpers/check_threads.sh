
threads=$1

if [ -z "$threads" ]
then
  threads=1
else
  re='^[0-9]+$'
  if ! [[ $threads =~ $re ]] ; then
     threads=1
  fi
fi
echo $threads 
