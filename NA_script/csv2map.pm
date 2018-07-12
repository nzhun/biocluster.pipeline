package csv2map;

use strict;
use warnings;
my $loadSet;
sub loadSet{
	my @names=();
	my $part=0;
	my $fullkey="";
	my @sets=split(",",$_[0]);
	#print $_[0];
	for(my $j=0;$j<@sets;$j++){
		my $key=$sets[$j];
		my $start=substr $key, 0,1;
		my $end=substr $key,-1;
	    if($part==0){
			 if($start eq "\"" && $end ne "\""){
				 if($fullkey eq ""){ 
					 $fullkey=$key;
				 }else{
				 	$fullkey=$fullkey.",".$key;
				}
					$part=1;				
		     }elsif( ($start ne "\"" && $end eq "\"")||($start eq "\"" && $end ne "\"")){
				 print "\nError: ".$start.$end;
				 exit;
		     }
			 else{
				 $key =~ s/\"//g;
				 push(@names,$key);
			 }
			 
	     }else{
			 if($fullkey eq ""){ 
				 $fullkey=$key;
			 }else{
			 	$fullkey=$fullkey.",".$key;
			}
			#print $key."\t$fullkey\n";
			 if($end ne "\""){
				 next;
			 }else{
				 $fullkey =~ s/\"//g;
				 push(@names,$fullkey);
				 $part=0;
				 $fullkey="";
	     	 }
	     }	
	}
	return \@names;
}
1;