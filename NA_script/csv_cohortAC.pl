
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use csv2map;


print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: CSVFile\n";
my $file=$ARGV[0];
sub traverse {
		my $file=$_[0];
		my $type=$_[1];
		my %sites=%{$_[2]};
		open IN,"$file" or die "Error: $file cannot find!\n";
		open OUT, ">$file.corrected.csv";
		print "Please find out in $file.corrected.csv\n";
		my $i=0;
		my $ALT="ALT";
		my $chr="CHROM";
		my $pos="POS";
		my $index1=0;
		my $index2=0;
		my $index3=0;
		my $index4=0;
		my $GT="GT";
		my @header=();		
		while(my $line=<IN>){
			chomp($line);
			$line =~ s/\r//g;
			if($i==0){
				@header=@{csv2map::loadSet($line)};
				my $j=0;
				foreach my $key (@header){
					#print $key."\n";
					if($key eq $chr){
						$index3=$j;
					}
					if($key eq $pos){
						$index4=$j;
					}
					if($key eq $ALT){
						$index1=$j;
					}
					if($key eq $GT){
						$index2=$j;
					}
					$j=$j+1;
				}
				if($type==2){
					print OUT $line.",Cohort_AC\n";
				}		
			    $i=$i+1;	
			}
			else{
				my @names=@{csv2map::loadSet($line)};
				my $alt=$names[$index1];
				my $gt=$names[$index2];
				my $chr=$names[$index3];
				my $pos=$names[$index4];
				my $add=1;
				my $key="$chr:$pos:$alt";
				if($type==1){
					if($gt =~ /[1-9]\/[1-9]/){$add=2;}	
			    	if(exists($sites{$key})){
						$sites{$key}=$sites{$key}+$add;
					}
			    	else{$sites{$key}=$add;}
				}else{
					print OUT $line.",".$sites{$key}."\n";
				}
			}	    	
		}
	close IN;
	close OUT;
    return \%sites;
}
		
my $type=1;
my %sites=();
my $addr_sites=traverse($file,$type,\%sites);
$type=2;
traverse($file,$type,$addr_sites);
