
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
my @outf=split(/\//,$file);
my %samples=();
open IN,"zcat $file|" or die "Error: $file cannot find!\n";
open OUT, ">$outf[@outf-1].2.bed";
print "Please find out in $file.2.bed\n";
my @header=();
my @csq=();
my $wi=0;
my $lc=0;
push(@header,"FILTER");
while(my $line=<IN>){
	chomp($line);
	#print $line."\n";
	if ($line =~ /^##INFO/){
		my @sets=split(/ID=|,/,$line);
		if($sets[1] eq "CSQ"){
			my @subsets=split(/:|\||\"/,$sets[@sets-1]);
			#print "sub\t".@subsets."\t".$sets[@sets-1]."\n"; 
			#exit;
			for(my $i=2;$i<@subsets-1;$i++){
			#	print $i."\t".$subsets[$i]."\n"; 
			#	exit;
				#push(@header,$subsets[$i]);
				push(@csq,$subsets[$i]);
			}
		}
		push(@header,$sets[1]);
#		@header=sort(@header);
	}elsif ($line =~ /^#C/){
		#my @sets=split(/\s+/,$line);
		#for(my $i=9;$i<@sets;$i++){
	#		$samples{$i}=$sets[$i];
	#	}
	
	}elsif ($line =~ /^#/){next}
	else{
		
		if($wi==0){
			@header= (@header,@csq);
		#	for(my $ii=0;$ii<@header;$ii++){
		#		print $ii."\t".$header[$ii]."\n";
		#	}
			print OUT "#CHROM\tFROM\tTO\tREF\tALT\t".join("\t",@header)."\n";
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
		#print "$chr\t$start\t$end\t$ref\t$sets[4]\n";
		my %names=%{INFO2map::loadINFO($sets[7])};
		my @fv=($sets[6]);
		$names{"FILTER"}=\@fv;
	    for(my $index=0; $index<@alts;$index=$index+1){
			my $alt=$alts[$index];	
			my %sub_names=%{INFO2map::extractINFO(\%names,scalar(@alts),$index)};
#			if($sub_names{"AC_NFE"}==0){next;}
			if($sub_names{"AF"} ne "."&&$sub_names{"AF"}>0.01){next;}
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
			$end=$start+length($nref)-1;
			my $outstr="$chr\t$start\t$end\t$nref\t$nalt";
			#my $outstr="$chr\t$start\t$end\t$ref\t$alt";
		#	print "outsearch\t".$outstr."\t".$index."\ttotal\t".@alts."\n";
			foreach my $item(@header){  ## simplify the output [(0..8,30..56,92,98..@header)]
			#	print $item."\n";
				if(exists($sub_names{$item})){
					if($sub_names{$item} eq ""){$sub_names{$item}="."}
					$outstr=$outstr."\t".$sub_names{$item};
					if($item eq "CSQ"){
							my @AAVs=split(/\|/,$sub_names{$item});
				#			print @AAVs."\t".$item."\t".$sub_names{$item}."\n";
							for(my $km=0;$km<@AAVs;$km++){
								$sub_names{$csq[$km]}=$AAVs[$km];
							}				
					}
				}else{
					$outstr=$outstr."\t.";
				}
			}
			my $global_accept=1;
			my $flag=0;
			#if($sub_names{"Consequence"} =~ /splice_donor|splice_acceptor_/){$flag=1;}
			#if($flag==0 && $sub_names{"Consequence"} =~ /UTR|intron|RNA|downstream|intergenic|non_coding|regulatory|TF|upstream/){next;}
			AF:{
				#print $sub_names{"Consequence"}."\n" ;
				foreach my $skey("ExAC_MAF","ExAC_AFR_MAF","ExAC_AMR_MAF","ExAC_EAS_MAF","ExAC_FIN_MAF","ExAC_NFE_MAF","ExAC_OTH_MAF","ExAC_SAS_MAF"){
					if($sub_names{$skey}) {
						my @af_sets=split(/:|\&/,$sub_names{$skey});
						#if(@af_sets<2){print "$outstr\n";next;}
						my $accept=1;
						for(my $jk=0;$jk<@af_sets;$jk=$jk+2){
							if($af_sets[$jk] eq $alt) {
								if($af_sets[$jk+1]>0.01){$accept=0;$global_accept=0;last AF;}
							}
							
						}
						if($accept==0){$global_accept=0;last AF;}
					}
			 }
	    }
		 if($global_accept==0){next}
		#	print "outsearch\t".$outstr."\t".$index."\ttotal\t".@alts."\n";
		  print OUT $outstr."\n";
			
		}
	}
	
}


close IN;
close OUT;
