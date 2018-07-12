echo "usage: bash split_info.sh vcf_file"
vcf=$1;
awk 'BEGIN{
   FS="\t";OFS="\t";
   for(i=1;i<4;i++){
      nms[i]="";
   }

}
{

  if($1 ~/^##/){next;}

   if($1 ~/^#/){
	 $1="CHR";    
         $8="AC\tAN\t1KPFreq\tExacFreq\tGene\tVarFunc\tvarClass\tAAchange\tCADDphred\tMetaSVMprd\tSIFT\tMutationTaster\tMutationAssessor\tPolyPhen2_HAR\tPhyloP\tGERP\tSiPhy\tProband\tInherit"
          print;
          for(i=10;i<NF+1;i++){
             nms[i-9]=$i;
          }
	  next;
	} 
	n=split($10,gts,":")
	n=split(gts[1],allele,"/");
	if(gts[3]<10){next}
	alt_c=1;
	
	if($5 ~/,/){
		if(allele[1]>1){alt_c=allele[1]}
		if(allele[2]>1){alt_c=allele[2]}
		split($5,a,",");
		$5=a[alt_c];
	}
	n=split($8,infos,";");
	freq1kg=0;
	ac=0;
	an=0;
	exacfreq=0;
	aachange=".";
	gene=".";
	varfun=".";
	meta=".";
	sift="."
	mutt="."
        pp2="."
	muta="."
	cadd=".";
        varclass=".";
	for(i=1;i<n+1;i++){
	   if(infos[i] ~/^1KGfreq=/) {
	       m=split(infos[i],b,"=" ); 
		   if(m>1){
		      freq1kg=b[2];
              split(freq1kg,g,",");
              freq1kg=g[alt_c]
		   }
		}
 	   if(infos[i] ~/^AC=/) {
 	       m=split(infos[i],b,"=" ); 
 		   if(m>1){
		      ac=b[2];
              split(ac,g,",");
              ac=g[alt_c]
	       }
 		}
  	   if(infos[i] ~/^AN=/) {
  	       m=split(infos[i],b,"=" ); 
  		   if(m>1){
		      an=b[2];
                 #split(an,g,",");
                #an=g[alt_c]
	       }
  		} 
	 
	   if(infos[i] ~/^VarFunc=/) {
	       m=split(infos[i],b,"=" ); 
		   if(m>1){
		      varfun=b[2];
              #split(varfun,g,",");
              #varfun=g[alt_c]
	        }
	   }

           if(infos[i] ~/^VarClass=/) {
	       m=split(infos[i],b,"=" ); 
		   if(m>1){
		      varclass=b[2];
              split(varclass,g,",");
              varclass=g[alt_c]
	        }
	   }

 	   if(infos[i] ~/^GeneName=/) {
         	m=split(infos[i],b,"="); 
   		   if(m>1){
		     gene=b[2];
		  }
	  }
	   if(infos[i] ~/^AAChange=/) {
    	        m=split(infos[i],b,"=" ); 
  		if(m>1){
		  	aachange=b[2];
                        split(aachange,g,",");
                        aachange=g[alt_c]
                }
	  }
	  
       if(infos[i] ~/^CADDphred=|^CADDInDelphred=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           cadd=b[2];
               split(cadd,g,",");
               cadd=g[alt_c]
             }
       }
	  
       if(infos[i] ~/^MetaSVMprd=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           meta=b[2];
               split(meta,g,",");
               meta=g[alt_c]
             }
       }
	  
       if(infos[i] ~/^SIFTprd=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           sift=b[2];
               split(sift,g,",");
               sift=g[alt_c]
             }
       }

       if(infos[i] ~/^MutTprd=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           mutt=b[2];
               split(mutt,g,",");
               mutt=g[alt_c]
             }
       }
      
       if(infos[i] ~/^MutAprd=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           muta=b[2];
               split(muta,g,",");
               muta=g[alt_c]
             }
       }

       if(infos[i] ~/^PP2.Hvar.prd=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           pp2=b[2];
               split(pp2,g,",");
               pp2=g[alt_c]
             }
       }
     
       if(infos[i] ~/^PhyloP=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           phylop=b[2];
               split(phylop,g,",");
               phylop=g[alt_c]
             }
       }
    
     
       if(infos[i] ~/^GERP=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           gerp=b[2];
               split(gerp,g,",");
               gerp=g[alt_c]
             }
       }      
    
       if(infos[i] ~/^SiPhy=/) {
   	        m=split(infos[i],b,"=" ); 
 	    	if(m>1){
	           siphy=b[2];
               split(siphy,g,",");
               siphy=g[alt_c]
             }
       }
     if(infos[i] ~/^ExACfreq=/) {
    	        m=split(infos[i],b,"=" ); 
  		if(m>1){
		   exacfreq=b[2];
           split(exacfreq,g,",");
           exacfreq=g[alt_c]
		}
      }
  }
  $8=ac"\t"an"\t"freq1kg"\t"exacfreq"\t"gene"\t"varfun"\t"varclass"\t"aachange"\t"cadd"\t"meta"\t"sift"\t"mutt"\t"muta"\t"pp2"\t"phylop"\t"gerp"\t"siphy

  if($10 !~/0\/0|^\./){$8=$8"\t"nms[1] } 

  f=0;

  if($11 !~/0\/0|^\./){f=1;$8=$8"\t"nms[2] }

  if($12 !~/0\/0|^\./){if(f==0){ $8=$8"\t"nms[3];}else{$8=$8","nms[3];} }

 #for(i=10;i<NF+1;i++){
  #   if($i ~/^\.|^0\/0/){$i=0;}   
  #   else if($i ~/^0\/1|^0\|1|^1\|0/){$i=1}
  #   else if($i ~/^1\/1|^1\|1/){$i=2}
  #}
  print ;
  
 }' $vcf > $vcf.list
