echo "usage: bash split_info.sh vcf_file"
vcf=$1;
#labels=$2;

zcat $vcf|awk -v labels=$2 'BEGIN{
	FS="\t";OFS="\t";
}
{  
       if($1 ~/^##/){next;}
       marker="";
       if(labels!=""){marker="\t"labels}
	if($1 ~/^#/){
		$1="CHR";    
		$8="AC\tAN\tLAC\tLAN\t1KPFreq\tExacFreq\tGene\tVarFunc\tvarClass\tAAchange\tCADD\tMarker\tMetaSVMprd\tSIFT\tMutationTaster\tMutationAssessor\tPolyPhen2_HAR\tPhyloP\tGERP\tSiPhy"
	     #   $8=$8""marker	
         	print;
		next;
	} 
	alc=1;
	
	if($5 ~/,/){alc=split($5,a,",");}
		alt_c=1
	for(i=1;i<alc+1;i++){
          lac[i]=0;
        }
	n=split($8,infos,";");
   
	freq1kg=0;
	ac=0;
	an=0;
	exacfreq=0;
	aachange=".";
	gene=".";
	varfun=".";
    varclass=".";
    CADD="NA";
	meta="."
	sift="."
	gerp="."
	phylop="."
	siphy="."
	pp2="."
	mutt="."
	muta="."
	  
	for(i=1;i<n+1;i++){
	   if(infos[i] ~/^1KGfreq=/) {
	       m=split(infos[i],b,"=" ); 
		   if(m>1){
		      freq1kg=b[2];
              split(freq1kg,g,",");
              freq1kg=g[1]
		   }
		}
 	   if(infos[i] ~/^AC=/) {
 	       m=split(infos[i],b,"=" ); 
 		   if(m>1){
		      ac=b[2];
	       }
 		}
  	   if(infos[i] ~/^AN=/) {
  	       m=split(infos[i],b,"=" ); 
  		   if(m>1){
		      an=b[2];
	       }
  		} 
	   if(infos[i] ~/^CADDphred=|^CADDInDelphred=/) {
  	       m=split(infos[i],b,"=" ); 
  		   if(m>1){
		      CADD=b[2];
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

        if(infos[i] ~/^PP2.hvar.prd=/) {
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
		
	  
	   if(infos[i] ~/^VarFunc=/) {
	       m=split(infos[i],b,"=" ); 
		   if(m>1){
		      varfun=b[2];
	        }
	   }

           if(infos[i] ~/^VarClass=/) {
	       m=split(infos[i],b,"=" ); 
		   if(m>1){
		      varclass=b[2];
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
                }
	  }
	  
	  
	   if(infos[i] ~/^ExACfreq=/) {
    	        m=split(infos[i],b,"=" ); 
  		if(m>1){
		   exacfreq=b[2];
		}
      }
  }
 
  lan=0; 
  for(i=10;i<NF+1;i++){
    if($i !~/^\./){lan=lan+1;}
    if($i !~/0\/0/){
       split($i,gtinfo,":");
       m=split(gtinfo[1],alleles,"/");
       if(m>1){
               
             if(alleles[1]!=0){ lac[alleles[1]]=lac[alleles[1]]+1;}
             if(alleles[2]!=0){ lac[alleles[2]]=lac[alleles[2]]+1;}
            
       }
    }
    if($i ~/^\.|^0\/0/){  $i=0;}   
    else if($i ~/^0\/1|^0\|1|^1\|0/){ $i=1}
    else if($i ~/^1\/1|^1\|1/){  $i=2}
  }
 lstr=lac[1]
for(i=2;i<alc+1;i++){
  lstr=","lstr
} 
$8=ac"\t"an"\t"lstr"\t"lan"\t"freq1kg"\t"exacfreq"\t"gene"\t"varfun"\t"varclass"\t"aachange"\t"CADD""marker"\t"meta"\t"sift"\t"mutt"\t"muta"\t"pp2"\t"phylop"\t"gerp"\t"siphy

 
 print ;
  
 }' | bgzip -c  > $vcf.list.gz
