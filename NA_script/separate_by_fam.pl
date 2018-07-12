#/usr/bin/perl
use strict;
use warnings;


## load ped file
my $fvcf="/home/local/ARCS/nz2274/BreastCancer/Variants/CUMC_RGN_BC_FreezeOne.hardfiltered.adjnm.mappability.gnomad.vcf.gz";
my $fped="/home/local/ARCS/nz2274/BreastCancer/DOCs/BCR_Freeze5.ped";
my %fam=();
my %mem=();
my $i=0;
open PED, $fped or die "$fped cannot find!";
while(my $line=<PED>){
	if($i==0){$i=$i+1;next}
	my @arr=();
	my @sets=split(/\s+/,$line);
	if(exists($fam{$sets[0]})){
		@arr=@{$fam{$sets[0]}};
	}
	$mem{$sets[1]}=$sets[0];
	push(@arr,$sets[1]);
	$fam{$sets[0]}=\@arr;
	
}
close PED;
my %handlers=();
my @keys=keys %fam;
foreach my $key(@keys){
	my @sets=@{$fam{$key}};

	if(@sets<2){
		delete($fam{$key});
		#print $key."\t singleton\n";
	}else{
		my $fh;
		open $fh,"> fam_$key.vcf";
		$handlers{$key}=$fh;
		
	}
}
print "we have ".scalar(keys %fam)."\n";


my %index=();
open VCF, "zless $fvcf|" or die "$fvcf cannot find!\n";
while(my $line=<VCF>){
	if($line =~ /^##/){
 	   foreach my $sfam (keys %fam){
 		   my $handle=$handlers{$sfam};
 	   	   print {$handle} $line;
 	   }
		next;
	}
	chomp($line);
	my @sets=split(/\s+/,$line);
	if($line =~ /^#C/){
		for(my $i=9;$i<@sets;$i++){
			$index{$sets[$i]}=$i;
		}
	}else{
		my @freqExac=split("ExACfreq=",$sets[7]);
		if(@freqExac >1){
			my @exacValue=split(";",$freqExac[1]);
			if(@exacValue >0){			
				if($exacValue[0]!~/,/ && $exacValue[0]>0.01){next;}
			}else{
				print $sets[0]."\t".$sets[1]."\n";
			}
		}
		my @freqKG=split("1KGfreq=",$sets[7]);
		if(@freqKG>1){
			my @kgValue=split(";",$freqKG[1]);
			if(@kgValue>0){
				if($kgValue[0] !~/,/ && $kgValue[0]>0.01){next;}
			}
		}
		my @func=split("VarFunc=",$sets[7]);
		if(@func>1){
			my @funcvalue=split(";",$func[1]);
			if(@funcvalue>0){
				if($funcvalue[0] !~ /^exonic|^splicing/){next;}
			}
		}
		
		my @mapset=split("Mappability=",$sets[7]);
		my @mapvalue=split(";",$mapset[1]);
		if($mapvalue[0]<1){next}
		
	}
	foreach my $sfam (keys %fam){
	   my $handle=$handlers{$sfam};
	   my @arr=@{$fam{$sfam}};
	   my @cols=(0..8);
	   my $flag=0;
	   foreach my $id(@arr){
		   if(exists($index{$id})){
			   push(@cols,$index{$id});
			   if($sets[$index{$id}] !~ /0\/0/ && $sets[$index{$id}] !~ /\.\/\./){$flag=1;}
	   		}
	   }
	   if($flag!=0){
	   		my $ostr=join("\t",@sets[@cols])."\n";
	   		print {$handle}  $ostr;
   		}
	}	
	
	
}
close VCF;

foreach my $sfam(keys %fam){
	my $handle=$handlers{$sfam};
	close $handle;
}

## only large family are required to be repaired


###  get the columns of the sample family

## open this file handle
## set up a map => family => handle
##  set up a map => family => column
### set up a map => col => family


## write the position mutatted in this family

## close this handle
