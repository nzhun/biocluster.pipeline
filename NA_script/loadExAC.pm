package loadExAC;

use strict;
use warnings;
my $load;
sub load{
	my $fexac="/home/local/ARCS/nz2274/Resources/ExAC_r03_0316_z_pli_rec_null_data.bed.gz";
	
	my %info=();
	open INDP, "zless $fexac|" or die "$fexac does not exist!\n";
	my $line=<INDP>;
	while($line=<INDP>){
		chomp($line);
		my @sets=split("\t",$line);
		my @gene_info=@sets[17..19];
		#print join(":",@sets[17..19])."\n";
		$info{$sets[4]}=\@gene_info;
	}
	close INDP;
	### misz LofZ PLI
	return \%info;
}
1;