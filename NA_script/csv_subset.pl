
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use csv2map;

print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: CSVFile\n";
my $file=$ARGV[0];
my $fsubset=$ARGV[1];
if(! -f $file){
	die "Error: $file cannot find!\n";
}
if(!-f $fsubset){
	die "Error: $fsubset cannot find!\n";
}

my %subset=();
open SIN,"$fsubset";
while(my $sline=<SIN>){
	chomp($sline);
	my @sets=split(/\s+/,$sline);
	$subset{$sets[0]}=1;
}
close SIN;
print "Totally ".scalar(keys %subset)."in $fsubset\n";
open IN,"$file" or die "Error: $file cannot find!\n";

open OUT, ">$file.corrected.csv";
print "Please find out in $file.corrected.csv\n";
my $i=0;

my $ALT="ALT";
my $GT="proband";
my %keymap=();
my @header=();
while(my $line=<IN>){
	chomp($line);
	$line =~ s/\r//g;
	
	if($i==0){
		@header=@{csv2map::loadSet($line)};
		my $j=0;
		foreach my $key (@header){
			$keymap{$key}=$j;
			$j=$j+1;
		}
		print OUT $line.",GT,AD,Ind_DP,AB,GQ,ProbandName\n";		
	    $i=$i+1;	
	}
	else{
		my @names=@{csv2map::loadSet($line)};
		my $temp_alt=$names[$keymap{"ALT"}];
		my $temp_gt=$names[$keymap{"proband"}];
		my @alts=split(",",$temp_alt);
		#print $temp_alt."\t".$temp_gt."\n";
		if(@alts==1){
			my @temps=split(/\(/,$temp_gt);
			my $ProbandName=$temps[0];
			if(!exists($subset{$ProbandName})){next;}
			#print "$ProbandName  in\n";
			my @infos=split(":",$temps[1]);
			my @ac_as=split(",",$infos[1]);
			my @alleles=split("/",$infos[0]);
			my $AD=$ac_as[$alleles[1]];
		    my $DP=$ac_as[$alleles[1]]+$ac_as[$alleles[0]];
			#print $temp_gt."\n";
			if($DP==0){next}
			my $AB=$AD/$DP;
			my $GQ=$infos[3];
			my $GT=$infos[0];
			print OUT $line.",$GT,\"$AD\",$DP,$AB,$GQ,$ProbandName\n";
			next;
		}else{
			my @temps=split(/\(/,$temp_gt);
			my $ProbandName=$temps[0];
			if(!exists($subset{$ProbandName})){next;}
			my @infos=split(":",$temps[1]);
			my @alleles=split("/",$infos[0]);
			my $index=$alleles[0];
			if($alleles[0]<$alleles[1]){
				$index=$alleles[1];
			}
		
			my @ac_as=split(",",$infos[1]);
			$infos[1]=$ac_as[$alleles[0]].",".$ac_as[$alleles[1]];
		    $names[$keymap{"proband"}]=$temps[0]."(".join(":",@infos);  
			
			my $AD=$ac_as[$alleles[1]];
		    my $DP=$ac_as[$alleles[0]]+$ac_as[$alleles[1]];
			my $AB=$AD/$DP;
			my $GQ=$infos[3];
			my $GT=$infos[0];
			
						
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
			
			print OUT "\"".join("\",\"",@names)."\"".",$GT,\"$AD\",$DP,$AB,$GQ,$ProbandName\n";
		}	    	
	}
	
}


close IN;
close OUT;