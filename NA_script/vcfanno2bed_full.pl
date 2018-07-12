
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use INFO2map;
sub cell2array {
	my @arr=($_[0]);
	return \@arr;
}

print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: vcf_file\n";
my $file=$ARGV[0];
my %samples=();
open IN,"zcat $file|" or die "Error: $file cannot find!\n";
open OUT, ">$file.2.txt";
print "Please find out in $file.2.txt\n";
my @header=();
my $wi=0;
my $lc=0;
push(@header,"FILTER");
while(my $line=<IN>){
	chomp($line);
	#print $line."\n";
	if ($line =~ /^##INFO/){
		my @sets=split(/ID=|,/,$line);
		#print $sets[1]."\n";
		push(@header,$sets[1]);
#		@header=sort(@header);
	}elsif ($line =~ /^#C/){
		my @sets=split(/\s+/,$line);
		for(my $i=9;$i<@sets;$i++){
			$samples{$i}=$sets[$i];
		}
	
	}elsif ($line =~ /^#/){next}
	else{
		
		if($wi==0){
			@header= (@header,"GeneName","Transcript","Exonloc","NucleotideChange","ProteinChange");

			print OUT "#CHROM\tFROM\tTO\tREF\tALT\t".join("\t",@header)."\tproband"."\tGenotype"."\n";
		}
		$wi+=1;
		my @sets=split(/\s+/,$line);
		my $chr=$sets[0];
		my $start=$sets[1];
		my $end=$sets[1]+length($sets[3])-1;
		#print $sets[1]."\t".$sets[3]."\n";
		my $ref=$sets[3];
		my @alts=$sets[4];
		if($sets[4] =~/,/){
			@alts=split(",",$sets[4]);
			
		}
		my %names=%{INFO2map::loadINFO($sets[7])};
                my $exon_key="Func.refGene";
	        if(!exists($names{$exon_key})){$exon_key="VarFunc";}
		my $exon="";
                if(exists($names{$exon_key})){
			my @exons=@{$names{$exon_key}};
			$exon=$exons[0];
		}
                if($exon !~ /exonic|splicing/){next}
		my @fv=($sets[6]);
		$names{"FILTER"}=\@fv;
	    for(my $index=0; $index<@alts;$index=$index+1){
			my $alt=$alts[$index];	
			my %sub_names=%{INFO2map::extractINFO(\%names,scalar(@alts),$index)};
			#my $outstr="$chr\t$start\t$end\t$ref\t$alt";
			my $nalt=$alt;
			my $nref=$ref;
			if(length($alt)>1 && length($ref)>1){
				if(length($alt)==length($ref)){
					$nref=substr($ref,0,1);
					$nalt=substr($alt,0,1);
				}elsif(length($alt)>length($ref)){
					$nref=substr($ref,0,1);
					$nalt=substr($alt,0,length($alt)-length($ref)+1)
				}else{
					
					$nalt=substr($alt,0,1);
					$nref=substr($ref,0,length($ref)-length($alt)+1)
					
				}
			}
			$end=$sets[1]+length($nref)-1;
			my $outstr="$chr\t$start\t$end\t$nref\t$nalt";
	#		if($sub_names{"AF"}>0.01){next;}
		#	print "outsearch\t".$outstr."\t".$index."\n";
			foreach my $item(@header){
				if(exists($sub_names{$item})){
					$outstr=$outstr."\t".$sub_names{$item};
					if(lc($item) eq "aachange"|| lc($item) eq "aachange.refgene"){
							my @AAVs=split(":",$sub_names{$item});
							if(@AAVs>4){
								$sub_names{"GeneName"}=$AAVs[0];
								$sub_names{"Transcript"}=$AAVs[1];
								$sub_names{"Exonloc"}=$AAVs[2];
								$sub_names{"NucleotideChange"}=$AAVs[3];
								$sub_names{"ProteinChange"}=$AAVs[4];
							}					
					}
				}else{
					$outstr=$outstr."\t.";
				}
			}
			for(my $id=9;$id<@sets;$id++){
				my $k=($index+1);
			    if($sets[$id] =~ /\/$k/){
			    	print OUT $outstr."\t".$samples{$id}."\t".$sets[$id]."\n";
			   }
			}
		}
	}
	
}


close IN;
close OUT;
