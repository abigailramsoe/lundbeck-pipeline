#!/usr/bin/perl

####
#José Víctor Moreno Mayar <morenomayar@gmail.com>#
####

##n means covered but tied under mindiff
##ATGC=high cov, high Q, maj
##atgc=high cov, low Q
##N=low/no cov
##n=no maj/alt low Q
##if cov=4, cov 4 IS taken into account
##ONLY take into account HIGHQ and HIGHCOV INDELs

use Getopt::Long;

&printUsage if (@ARGV < 1);

my %opts=();
GetOptions (\%opts,'i=s', 'o=s', 'l=i', 'cov=i', 'diff=f', 'h=s', 'idiff=f', 'callindels=s');

die("Incomplete set of arguments") if (!$opts{i} || !$opts{o} || !$opts{l} || !$opts{cov} || !$opts{diff} || !$opts{h} || !$opts{idiff} || !$opts{callindels});


for($i=0; $i<$opts{l}; $i++){
	$CNS[$i]="N";
}
$genomesize=$opts{l};
$mincov=$opts{cov};
$mindiff=$opts{diff};
$imindiff=$opts{idiff};
open(OUT, ">$opts{o}");
print "Min_depth=$mincov\nMin_diff=$mindiff\nIndel_min_diff=$imindiff\nInitial_size=$genomesize\n";
print "Used the following SNPs:\n";
@INDELS=();
@HETINDELS=();
open(VCF, "$opts{i}");
while(<VCF>){
	if($_ !~ /\#/){
		$INDEL=0;
		$HETNR=0;
		$CNSBase="N";
		$diff=0;
		$depth=0;
		$highQ=0;
		$line=$_;
		chomp($line);
		@site=split(/\t/, $line);
		$pos=$site[1];
		$ref=$site[3];
		$alt=$site[4];
		if($alt =~ /\w+,\w+/){
			$HETNR=1;
		}
		$info=$site[7];
		@INFO=split(/\;/, $info);
		foreach $field (@INFO){
			if($field =~ /DP4/){
				$depth=$field;
			}
		}
		
		if($depth==0){
			foreach $field (@INFO){
				if($field =~ /DP/){
					$depth=$field;
				}
			}
			if($depth =~ /DP=(\d+)/){
				$RefAlleles=$1;
				$DEPTH=$1;
			}
		}
		
		if($depth =~ /DP4=(\d+),(\d+),(\d+),(\d+)/){
			$highQ=1;
			$RefAlleles=$1+$2;
			$AltAlleles=$3+$4;
			$DEPTH=$1+$2+$3+$4;
		}
		if($INFO[0] eq "INDEL"){
			$INDEL=1;
		}
		
		
		if($HETNR==0 && $INDEL==0){
			if($DEPTH>=$mincov && $highQ==1){
				if($RefAlleles/$DEPTH > $mindiff){
					$diff=$RefAlleles/$DEPTH;
					$CNSBase=$ref;
				}else{
					if($AltAlleles/$DEPTH > $mindiff && $alt ne "\."){
						$diff=$AltAlleles/$DEPTH;
						$CNSBase=$alt;
						print "$pos\t$ref\t$alt\t$DEPTH\t$RefAlleles\t$AltAlleles\t$diff\t$CNSBase\n";
					}else{
						$CNSBase="n";
					}
				}
			}else{
				if($DEPTH>=$mincov && $highQ==0){
					$CNSBase=lc $ref;
				}
			}
			$CNS[$pos-1]=$CNSBase;
		}
		
##Multiallelic sites
		if($HETNR==1 && $INDEL==0){
			if($alt =~ /(\w+)((,\w+)+)/){
				$alt=$1
			}
			if($DEPTH>=$mincov && $highQ==1){
				if($RefAlleles/$DEPTH > $mindiff){
					$diff=$RefAlleles/$DEPTH;
					$CNSBase=$ref;
				}else{
					if($AltAlleles/$DEPTH > $mindiff && $alt ne "\."){
						$diff=$AltAlleles/$DEPTH;
						$CNSBase=$alt;
						print "$pos\t$ref\t$alt\t$DEPTH\t$RefAlleles\t$AltAlleles\t$diff\t$CNSBase\tMultiAllelic\n";
					}else{
						$CNSBase="n";
					}
				}
			}else{
				if($DEPTH>=$mincov && $highQ==0){
					$CNSBase=lc $ref;
				}
			}
			$CNS[$pos-1]=$CNSBase;
		}
		if($INDEL==1){
			if($opts{callindels} eq "yes"){
			push @INDELS, $line;
			}
		}		
	}
}
close(VCF);
@OriCNS=@CNS;

if($opts{callindels} eq "yes"){
##INDELS##
$shift=0;
$substshift=0;
print "Used the following INDELs:\n";
foreach $INDEL (@INDELS){
	$depth==0;
	$HETNR=0;
	$R=0;
	$NR=0;
	$highQ=0;
	@site=split(/\t/, $INDEL);
	$pos=$site[1];
	$ref=$site[3];
	$alt=$site[4];
	if($alt =~ /\w+,\w+/){
		$HETNR=1;
	}
	$info=$site[7];
	@INFO=split(/\;/, $info);
	foreach $field (@INFO){
		if($field =~ /DP4/){
			$depth=$field;
		}
	}
		
	if($depth==0){
		foreach $field (@INFO){
			if($field =~ /DP/){
				$depth=$field;
			}
		}
		if($depth =~ /DP=(\d+)/){
			$RefAlleles=$1;
			$DEPTH=$1;
		}
	}
		
	if($depth =~ /DP4=(\d+),(\d+),(\d+),(\d+)/){
		$highQ=1;
		$RefAlleles=$1+$2;
		$AltAlleles=$3+$4;
		$DEPTH=$1+$2+$3+$4;
	}
	
	if($HETNR==0 && $alt ne "\."){
		if($DEPTH>=$mincov && $highQ==1){
			if($RefAlleles/$DEPTH > $imindiff){
				$diff=$RefAlleles/$DEPTH;
				$R=1;
			}else{
				if($AltAlleles/$DEPTH > $imindiff){
					$diff=$AltAlleles/$DEPTH;
					$NR=1;
				}
			}
		}else{
			$R=1;
		}
		if($NR==1){
			$Rlength=length($ref);
			$Alength=length($alt);
			$shift=$Rlength-$Alength;
			$newsize=$genomesize-($shift);
			@modCNS=();#
			for($i=0;$i<($pos-1+$substshift);$i++){
				$modCNS[$i]=$CNS[$i];
			}
			#$f=$pos-1+$substshift;
			#print "Empiezo en $f\n";
			@altI=split(//, $alt);
			for($i=0; $i<$Alength; $i++){
				$modCNS[($pos-1+$substshift+$i)]=$altI[$i];
				#print "$altI[$i]";
			}
			print "\n";
			for($i=$pos-1+$substshift+$Rlength;$i<$genomesize;$i++){
				#$nuevita=$i-($Rlength-$Alength);
				#print "voy en $i de CNS que es una $CNS[$i] y se la agregue a la pos $nuevita de modCNS\n";
				$modCNS[$i-($shift)]=$OriCNS[$i-$substshift];
			}
			#$f=$pos-1+$substshift+$Rlength;
			#print "Termino en $f\n";
			print "$pos\t$ref\t$alt\t$DEPTH\t$RefAlleles\t$AltAlleles\t$diff\t$genomesize\t$newsize\t$Rlength\t$Alength\t";
			print @altI;
			#print "\n$shift\t$substshift\t";
			$substshift=$substshift-$shift;
			#print "$substshift\n";
			print "\n";
			@CNS=@modCNS;
			$genomesize=$newsize;
			
		}		
	}

##Diallelic INDELS
if($HETNR==1 && $alt ne "\."){
	if($alt =~ /(\w+)((,\w+)+)/){
		$HETNR=1;
		$alt=$1;
	}
		if($DEPTH>=$mincov && $highQ==1){
			if($RefAlleles/$DEPTH > $imindiff){
				$diff=$RefAlleles/$DEPTH;
				$R=1;
			}else{
				if($AltAlleles/$DEPTH > $imindiff){
					$diff=$AltAlleles/$DEPTH;
					$NR=1;
				}
			}
		}else{
			$R=1;
		}
		if($NR==1){
			$Rlength=length($ref);
			$Alength=length($alt);
			$shift=$Rlength-$Alength;
			$newsize=$genomesize-($shift);
			@modCNS=();#
			for($i=0;$i<($pos-1+$substshift);$i++){
				$modCNS[$i]=$CNS[$i];
			}
			#$f=$pos-1+$substshift;
			#print "Empiezo en $f\n";
			@altI=split(//, $alt);
			for($i=0; $i<$Alength; $i++){
				$modCNS[($pos-1+$substshift+$i)]=$altI[$i];
			}
			print "\n";
			for($i=$pos-1+$substshift+$Rlength;$i<$genomesize;$i++){
				#$nuevita=$i-($Rlength-$Alength);
				#print "voy en $i de CNS que es una $CNS[$i] y se la agregue a la pos $nuevita de modCNS\n";
				$modCNS[$i-($shift)]=$OriCNS[$i-$substshift];
			}
			#$f=$pos-1+$substshift+$Rlength;
			#print "Termino en $f\n";
			print "$pos\t$ref\t$alt\t$DEPTH\t$RefAlleles\t$AltAlleles\t$diff\t$genomesize\t$newsize\t$Rlength\t$Alength\t";
			print @altI;
			#print "\n$shift\t$substshift\t";
			$substshift=$substshift-$shift;
			print "\tMultiallelic\n";
			#print "$substshift\tMultiallelic\n";
			@CNS=@modCNS;
			$genomesize=$newsize;
			
		}		
	}
##

}
}

##
print "The final size is: $genomesize\n";
print OUT "\>$opts{h}\n";
#for($i=0; $i<$genomesize; $i++){
#	print "Si estoy entrando\n";
#	print OUT "$CNS[$i]";
#}
print OUT @CNS;#
print OUT "\n";

close(OUT);

sub printUsage{
    print "
USAGE: $0 <arguments>\n\nArguments:\n\n\t-i\tSTR\tvcf obtained with samtools mpileup -g and bcftools -cg\n\t-o\tSTR\tOutput fasta file\n\t-l\tINT\tLength of the reference that was used to map the reads (not multifasta)\n\t-cov\tINT\tMinimum depth of coverage required to call a cns base\n\t-diff\tFLOAT\tMinumum fraction of concordant read bases to call a cns base\n\t-idiff\tFLOAT\tMinumum fraction of concordant read bases to call a cns indel\n\t-h\tSTR\tFasta header for the output file\n\t-callindels\tSTR\tIf set to \"yes\", indels are called. Other values turn this option off.\n\nSample run: ./CnsMaj3.pl -i sample.vcf -o sample.fasta -l 16569 -cov 5 -diff 0.7 -idiff 0.7 -h sample\n\nNote: All of the arguments are required.\n\n";
exit 0;
}




