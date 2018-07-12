package loadGeneInfo;

use strict;
use warnings;
my $load;
sub load{
	my $flung="/home/local/ARCS/nz2274/Resources/GeneExpression/Lung_rank_RNAseq Asselin-Labat-GPL13112_human.csv";
	my $fheart="/home/local/ARCS/nz2274/Resources/GeneExpression/mousebrain_Heart.csv";
	my $fdiap="/home/local/ARCS/nz2274/Resources/GeneExpression/diaphragm_rank.csv";
	my %info=();
	open INDP, $fdiap or die "$fdiap does not exist!\n";
	my $line=<INDP>;
	while($line=<INDP>){
		chomp($line);
		my @gene_info=(-1,-1,-1);
		my @sets=split(",",$line);
		if(exists($info{$sets[0]})){
			@gene_info=@{$info{$sets[0]}};
			$gene_info[0]=$sets[1];
		}else{
			$gene_info[0]=$sets[1]
		}
		$info{$sets[0]}=\@gene_info;
	}
	close INDP;
	
	
	open INH, $fheart or die "$fheart does not exist!\n";
	$line=<INH>;
	while($line=<INH>){
		chomp($line);
		my @gene_info=(-1,-1,-1);
		my @sets=split(",",$line);
		if(exists($info{$sets[0]})){
			@gene_info=@{$info{$sets[0]}};
			$gene_info[1]=$sets[12];
		}else{
			$gene_info[1]=$sets[12]
		}
		$info{$sets[0]}=\@gene_info;
	}
	close INH;
	
	open INL, $flung or die "$flung does not exist!\n";
	$line=<INL>;
	while($line=<INL>){
		chomp($line);
		my @gene_info=(-1,-1,-1);
		my @sets=split(",",$line);
		if(exists($info{$sets[1]})){
			@gene_info=@{$info{$sets[1]}};
			$gene_info[2]=$sets[3];
		}else{
			$gene_info[2]=$sets[3];
		}
		$info{$sets[1]}=\@gene_info;
	#	print $sets[1]."\t".join("\t",@{$info{$sets[1]}})."\t".join("\t",@gene_info)."\n";
		#exit;
	
	}
	close INL;
	### info: Diaph  Heart Lung
	
	return \%info;
}
1;