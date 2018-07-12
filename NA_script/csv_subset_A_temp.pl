
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use csv2map;
use loadGeneInfo;
use loadExAC;


sub gene_anno {
	my ($addr_names,$addr_map,$outstr,$i,$addr_exp,$addr_exac)=@_;
	my %gene_exp=%{$addr_exp};
	my %gene_ExAC=%{$addr_exac};
	my @names=@{$addr_names};
	my %map=%{$addr_map};
	
	
	my @per_gene=(-1,-1,-1);
	#if(exists($gene_exp{$names[$map{"GeneName"}]})){
#		@per_gene=@{$gene_exp{$names[$map{"GeneName"}]}};
#	}
#	if(@per_gene <3){
#		print $names[$map{"ProbandName"}]."\n";
#		exit;
#	}
	my @cells=("Diaphragm","Heart","Lung");
	my $j=0;
	foreach my $cell(@cells){
		if($i==1){
			$outstr.="\t$cell";
			$map{$cell}=@names;
		}
		my $value=$per_gene[$j];
		
		push(@names,$value);
		$j+=1;
	}
	my @exac_agene=(-1,-1,-1);
#	if(exists($gene_ExAC{$names[$map{"GeneName"}]}) ){
#		my $addr=$gene_ExAC{$names[$map{"GeneName"}]};
#		@exac_agene=@{$addr};
#	}
#	if(@exac_agene <3){
#		print "2  ".$names[$map{"ProbandName"}]."\n";
#		exit;
#	}
	@cells=("misz","lofz","pLI");
	$j=0;
	foreach my $cell(@cells){
		if($i==1){
			$outstr.="\t$cell";
			$map{$cell}=@names;
		}
		my $value=$exac_agene[$j];
		push(@names,$value);
		$j+=1;
	}
		
	return (\@names,\%map,$outstr);
}



sub annotate {
	my ($addr_names,$addr_map,$outstr,$i)=@_;
	my @names=@{$addr_names};
	my %map=%{$addr_map};
	my @alts=split(",",$names[$map{"ALT"}]);
	my @temps=split(/\(/,$names[$map{"proband"}]);
	my @infos=split(":",$temps[1]);
	my $probandName=$temps[0];
	my $AB="";
	my $GT="";
	my $GQ="";
	my $DP="";
	my $AD="";
	if(@alts==1){
		my @ac_as=split(",",$infos[1]);
		my @alleles=split("/",$infos[0]);
		$AD=$ac_as[$alleles[1]];
	    $DP=$ac_as[$alleles[1]]+$ac_as[$alleles[0]];
		if($DP<5){return;}
		$AB=$AD/$DP;
		$GQ=$infos[3];
		$GT=$infos[0];
		
	}else{
		my @alleles=split("/",$infos[0]);
		my $index=$alleles[0];
		if($alleles[0]<$alleles[1]){
			$index=$alleles[1];
		}
	
		my @ac_as=split(",",$infos[1]);
		$infos[1]=$ac_as[$alleles[0]].",".$ac_as[$alleles[1]];
	    $names[$map{"proband"}]=$temps[0]."(".join(":",@infos).")";  
		$AD=$ac_as[$alleles[1]];
	    $DP=$ac_as[$alleles[0]]+$ac_as[$alleles[1]];
		if($DP<5){return;}
		$AB=$AD/$DP;
		$GQ=$infos[3];
		$GT=$infos[0];
	    for(my $k=0;$k<(@names);$k++){
			my $str=$names[$k];
			my @sets=split(",",$str);
			if(@sets == @alts){
				$names[$k]=$sets[$index-1];
			}
	    }	
	}
	if($i==1){
	  $map{"GT"}=@names;
	  $outstr.="\tGT";
	}
	push(@names,$GT);
	if($i==1){
		$map{"AD"}=@names;
		$outstr.="\tAD";
	}
	push(@names,$AD);
	if($i==1){
	$map{"Ind_DP"}=@names;
	$outstr.="\tInd_DP";
	}
	push(@names,$DP);
	if($i==1){
	$map{"AB"}=@names;
	$outstr.="\tAB";
	}
	push(@names,$AB);
	if($i==1){
		$map{"GQ"}=@names;
		$outstr.="\tGQ";
	}
	push(@names,$GQ);
	if($i==1){
		$map{"ProbandName"}=@names;
		$outstr.="\tProbandName";
	}	
	push(@names,$probandName);


    return (\@names, \%map,$outstr)	;		
	
}


### output coloumns


print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: CSVFile\n";
#my $fheader="/home/local/ARCS/nz2274/Pipeline/NA_script/csv_subset.header.txt";

my $file=$ARGV[0];
my $fheader=$ARGV[1];
if(! -e $file){
	die "Error: $file cannot find!\n";
}
if(!-e $fheader){#
	die "Error: $fheader cannot find!\n";
}

my %outkeys=();
open INH, "$fheader" or die "Error: $fheader cannot find!\n";
while(my $line=<INH>){
	chomp($line);
	$outkeys{$line}=1;
}
close INH;

open IN,"$file" or die "Error: $file cannot find!\n";

open OUT, ">$file.corrected_temp.csv";
print "Please find out in $file.corrected_temp.csv\n";



my $i=0;

my %map=();
my $freq="ExAC.nfe.freq";
my @header=();
my @indexs=();
my $outstr="";
my $flag=0;
my $addr_exp=loadGeneInfo::load();
my $addr_exac=loadExAC::load();

while(my $line=<IN>){
	chomp($line);
	$line =~ s/\r//g;
	
	if($i==0){
		@header=@{csv2map::loadSet($line)};
		my $j=0;
		foreach my $key (@header){
			$map{$key}=$j;
			if(exists($outkeys{$key})){
				$outstr=$outstr."\t".$key;
				push(@indexs,$j);
			}
			$j=$j+1;
		}
		$outstr=~ s/^\s+//;
			
		$i=$i+1;	
		next;
	}
	
    
	if($i>0){
		my @names=@{csv2map::loadSet($line)};
		my $alt=$names[$map{"ALT"}];
		my $gt=$names[$map{"proband"}];
		my $init_len=@names;
		if(!exists($map{"AB"})){
			$flag=1;	
		}
		if($flag==1){
			my @rt=annotate(\@names,\%map,$outstr,$i);
			#print @rt;
			if(@rt<1){next}
			my ($addr_name,$addr_map,$outstr1)=@rt;
			$outstr=$outstr1;
			@names=@{$addr_name};
			%map=%{$addr_map};
		}
		my $varfunc=$names[$map{"VarFunc"}];
    	if($varfunc && $varfunc ne "exonic" && $varfunc ne "exnoic,splicing"  && $varfunc ne "splicing,exnoic"  && $varfunc ne "splicing"){
			next;
		}
		### 1kg freq 
		my $kgfreq=$names[$map{"1KGfreq"}];
		#if ($kgfreq ne "." && $kgfreq ne "NA" && $kgfreq >0.01){next}
		
		##gene ##
		my $gene=$names[$map{"GeneName"}];
		if($gene eq "MUC" || $gene eq "HLA") {next;}
	#	print $init_len."\t".@names."\t".$map{"Ind_DP"}."\t".$names[197]."\n";
		my $dp=$names[$map{"Ind_DP"}];
		if($dp<10){next;}
		my $ad=$names[$map{"AD"}];
		if($ad<5){next}
	
		my $gq=$names[$map{"GQ"}];
		if($gq<30){next}
	    
		my $gt0=$names[$map{"GT"}];
		my $ab=$names[$map{"AB"}];
		if($gt0 eq "0/1" && $ab<0.1){
			next;
		}
		
		my $fs=$names[$map{"FS"}];
		if($fs>25){next}
		my $qd=$names[$map{"QD"}];
		if($qd<2){next}
		my $LAF=$names[$map{"AF"}];
		#if($LAF>0.1){next;}
		
		my ($addr_name,$addr_map,$outstr)=gene_anno(\@names,\%map,$outstr,$i,$addr_exp,$addr_exac);
		@names=@{$addr_name};
		%map=%{$addr_map};
		if($i==1){
			print OUT $outstr."\n";	
			
		}
		$i=$i+1;
		my $end_len=@names;
		print OUT join ("\t", @names[@indexs])."\t";
		if($init_len!=$end_len){
			print OUT join("\t",@names[$init_len..($end_len-1)]);
	    }	
		print OUT "\n";
		    	
	}
	
}


close IN;
close OUT;
