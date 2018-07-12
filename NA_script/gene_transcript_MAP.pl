#!/user/bin/perl

use strict;
use warnings;

print "Input format: perl gene_transcript_MAP.pl \$file \$index (start from 0)\n";
my $file="/home/local/ARCS/nz2274/Resources/Transcript/genecodeV19.gene.info.txt";
my $fgene=$ARGV[0];
my $index=$ARGV[1];
if(@ARGV<2){die "input is not correct!\n input format: perl gene_transcript_MAP.pl \$file \$index\n";}
my $fout=$fgene.".gene.txt";
#my $fout="/home/local/ARCS/nz2274/PSAP/script/MutationRate_correctMappability_addpDmis.gene.txt";

my %transcript=();
my %geneID=();
open IN,"$file" or die "$file cannot open!\n";
my $line=<IN>;
while($line=<IN>){
	chomp($line);
	my @sets=split(/\s+/,$line);
	my @trans=split(/\./,$sets[1]);
	my @geneIDs=split(/\./,$sets[0]);
	#print $trans[0]."\t".$sets[4]."\n";
	$transcript{$trans[0]}=$sets[4];
	$geneID{$geneIDs[0]}=$sets[4];
}
close IN;

print "start\n";

open OUT, "> $fout";

open GN," $fgene" or die "$fgene cannot open!\n";
my $count=0;
my @fails=();
while(my $linegene=<GN>){
	chomp($linegene);
    if($linegene =~ /^#/){print OUT "#GeneName\t".$linegene."\n";next}
	my @gene_sets=split("\t",$linegene);
    my $gene="-";
	if(exists($transcript{$gene_sets[$index]})){
    	$gene=$transcript{$gene_sets[$index]};
    }else{
		if(exists($geneID{$gene_sets[$index]})){
			$gene=$geneID{$gene_sets[$index]};
		}else{
			 push(@fails,$gene_sets[$index]); # ." cannot find the gene with given GeneID and TranscriptID!\n";
			$count=$count+1;
    	}
    }
	print OUT $gene."\t".$linegene."\n";
}
close GN;
close OUT;
print join(",",@fails)." totally $count geneID are not in $file\n";
