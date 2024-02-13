fasta_file=$1

check_region() {
    if [[ $1 != "chrM" ]] && [[ $1 != "MT" ]]; then
        echo "Could not find what the MT reference is called automatically"
        echo "Please supply the region name with -r"
        exit 1
    fi
}

region=$(grep ">" $fasta_file|grep M|head -n1|cut -c2-|cut -f1 -d" ")
check_region $region
echo $region
