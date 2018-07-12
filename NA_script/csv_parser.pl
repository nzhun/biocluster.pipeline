
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use csv2map;


print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: CSVFile\n";
my $file=$ARGV[0];
open IN,"$file" or die "Error: $file cannot find!\n";
open OUT, ">$file.corrected.csv";
print "Please find out in $file.corrected.csv\n";
my $i=0;

my $ALT="ALT";
my $index1=0;
my $index2=0;
my $GT="proband";
my @header=();
while(my $line=<IN>){
	chomp($line);
	
	
	if($i==0){
		@header=@{csv2map::loadSet($line)};
		my $j=0;
		foreach my $key (@header){
			#print $key."\n";
			if($key eq $ALT){
				$index1=$j;
			}
			if($key eq $GT){
				$index2=$j;
			}
			$j=$j+1;
		}
		print OUT $line."\n";		
	    $i=$i+1;	
	}
	else{
		my @names=@{csv2map::loadSet($line)};
		my $temp_alt=$names[$index1];
		my $temp_gt=$names[$index2];
		my @alts=split(",",$temp_alt);
		#print $temp_alt."\t".$temp_gt."\n";
		if(@alts==1){
			print OUT $line."\n";
			next;
		}else{
			my @temps=split(/\(/,$temp_gt);
			my @infos=split(":",$temps[1]);
			my @alleles=split("/",$infos[0]);
			my $index=$alleles[0];
			if($alleles[0]<$alleles[1]){
				$index=$alleles[1];
			}
		
			my @ac_as=split(",",$infos[1]);
			$infos[1]=$ac_as[$alleles[0]].",".$ac_as[$alleles[1]];
		    $names[$index2]=$temps[0]."(".join(":",@infos);  
			#print $line."\n";
		    for(my $k=0;$k<(@names);$k++){
				my $str=$names[$k];
			#	print $k."\t".$str."\n";
				my @sets=split(",",$str);
				#print "$k  vary ".@sets."\t".@alts."\n";
				if(@sets == @alts){
					#print "has" .$k."\t".$str."\n";
					$names[$k]=$sets[$index-1];
					#print $k."\t".$names[$k]."\n";
				}
				
		    }
			print OUT "\"".join("\",\"",@names)."\""."\n";
		}	    	
	}
	
}


close IN;
close OUT;
