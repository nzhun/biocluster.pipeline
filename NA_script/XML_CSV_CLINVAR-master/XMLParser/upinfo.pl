#!/usr/bin/perl

#!/usr/nom/env perl 


#         my $alpha=1;
	my $in=$ARGV[0];##DBSNP FILE  "../original/dbsnp.Esp5400.union.noxy.txt
	my $in2=$ARGV[1]; ##local FILE "378Sample.vcf"
	my $out2f=$ARGV[2];##out2put file
	
	if(!(-e $in && -e $in2)){print "error!not all file exists!";exit();}
	open inrare, $in;#"zcat $in "; 
	open inorg,"$in2";
	open out2,"> $out2f";
	print $out2f."\n";
	my $org;
	my $rarein;
	while(($org=<inorg>)=~/^#/){print out2 $org;}
	do{$rare=<inrare>;}while($rarein=~/^#/);
	
        while($org){  ## local file
	    chomp($org);
	    my @orgstrs=split("\t",$org);
	    my $orgchr=$orgstrs[0];    ##original file format: chr pos id ref alt
	    my $orgpos=$orgstrs[1];
	    my $orgref="";
	    my $orgalt="";
    ##******complicated format 
	    $orgref=$orgstrs[2];
	    $orgalt=$orgstrs[3];
            my $id=$orgstrs[9];
	    my @rarestrs;
	    my $rarechr="";
	    my $rarepos="";
	    my $rareref="";
	    my $rarealt="";
	    my $rid=0;
	    my $i=0;
	    do{
		@rarestrs=split("\t",$rare);
		$rarechr=$rarestrs[0];			##rare format gmaf, chr, pos, id ,ref,alt.	
		$rarepos=$rarestrs[1];
		$rarealt=$rarestrs[3];
		my @rids=split(/=/,$rarestrs[7]);
		$rid=$rids[1];
		chomp($rid);
	#	print "rare:  ".$rid."\t com  ".$id."\n";
	     }while(($rid<$id)&&($rare=<inrare>));
    	    if($rare){
    	       #   print $rare."\n";
	          if($rid==$id){
	                    $orgstrs[0]=$rarestrs[0];
	                    $orgstrs[1]=$rarestrs[1];
	                    $orgstrs[2]=$rarestrs[3];
	                    $orgstrs[3]=$rarestrs[4];
	                    print out2 join("\t",@orgstrs)."\n";			     
			    $org=<inorg>;
			    $rare=<inrare>;
	       }
	      if($id<$rid){
		print out2 "$org\n";	
		$org=<inorg>;	
      	    }
            }else {
              
		print out2 "$org\n";	
		$org=<inorg>;	
            }
         }

	close out2;
	close inorg;
	close inrare;

