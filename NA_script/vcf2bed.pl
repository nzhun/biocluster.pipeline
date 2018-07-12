
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use INFO2map;


print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: vcf_file\n";
my $file=$ARGV[0];
open IN,"zcat $file|" or die "Error: $file cannot find!\n";
open OUT, ">$file.2.bed";
print "Please find out in $file.2.bed\n";
my @header=();
my $wi=0;
while(my $line=<IN>){
	chomp($line);
	#print $line."\n";
	if ($line =~ /^##INFO/){
		my @sets=split(/ID=|,/,$line);
		#print $sets[1]."\n";
		push(@header,$sets[1]);
	}elsif ($line =~ /^#/){next}
	else{
		if($wi==0){

			print OUT "#CHROM\tFROM\tTO\tREF\tALT\t".join("\t",@header)."\n";
		}
		$wi+=1;
		my @sets=split(/\s+/,$line);
		my $chr=$sets[0];
		my $start=$sets[1]-1;
		my $end=$sets[1]+length($sets[3])-1;
		#print $sets[1]."\t".$sets[3]."\n";
		my $ref=$sets[3];
		my @alts=$sets[4];
		if($sets[4] =~/,/){
			@alts=split(",",$sets[4]);
		}
	    for(my $index=0; $index<@alts;$index=$index+1){
			my $alt=$alts[$index];	
			my %names=%{INFO2map::loadINFO($sets[7],scalar(@alts),$index)};
			my $outstr="$chr\t$start\t$end\t$ref\t$alt";
			foreach my $item(@header){
				if(exists($names{$item})){
					$outstr=$outstr."\t".$names{$item};
				}else{
					$outstr=$outstr."\t.";
				}
			}
			print OUT $outstr."\n";
		}
	}
	
}


close IN;
close OUT;
