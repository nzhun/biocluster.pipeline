### filter GQ and DP, 90% >9 and GQ>30

zcat $1 |awk 'BEGIN{
FS="\t";
OFS="\t";
}
{
if($1 ~/^#/){print;next;}
dct=0;
gct=0;
GT=1;
AD=2;
DP=3;
GQ=4;
n=split($9,gts,":");
for(i=1;i<n+1;i++){
   if(gts[i]=="GT"){GT=i;}
   if(gts[i]=="AD"){AD=i;}
   if(gts[i]=="DP"){DP=i;}
   if(gts[i]=="GQ"){GQ=i;}  	   
}
tct=0;
hct=0;
ct=0;


for(i=10;i<NF+1;i++){
   m=split($i,infos,":");
   if(infos[GT] !~/^\./){
	   ct=ct+1;
	   dp=infos[DP];
	   ad=infos[AD];
	   n=split(ad,ac,",");
	   af=0;
	   if(dp>0&&n>1){
		   min=ac[1]>ac[2]?ac[2]:ac[1];
		   af=min/dp;
		   if(min>2){tct=tct+1}
	   }
	   gq=infos[GQ];
	   if(dp>8){dct=dct+1;}
	   if(gq>29){gct=gct+1;}
	   #if(af>0.2||af==0.2){tct=tct+1;}
	   if(infos[GT] ~/0[\/|\|][1-9]/||info[GT] ~/[1-9][\/|\|]0/){hct=hct+1;}
    }
}

#if(hct==0){print}
if(dct/ct<0.9||gct/ct<0.9||(hct>0&&tct/hct<0.9)){next}
	print 

}
' |bgzip -c  >$2
