
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

my @header=();
while(my $line=<IN>){
	chomp($line);
	$line =~ s/\r//g;
	
	if($i==0){
		@header=@{csv2map::loadSet($line)};
		print OUT join("\t",@header)."\n"; 		
	    $i=$i+1;	
	}
	else{
		my @names=@{csv2map::loadSet($line)};
		#my $strl="\"".join("\"\t\"",@names)."\"\n" ;	
		my $strl=join("\t",@names)."\n" ;			
		print OUT $strl;
			
	}	    	
}
	



close IN;
close OUT;
