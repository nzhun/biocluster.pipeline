#setwd("/home/local/ARCS/nz2274/");#"~/server/PAH/PAH_2017/Relatedness/")
setwd("~/server/")
library("scatterplot3d")
#uff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.gz.plus.HapMap"

#folder<-"~/server/BreastCancer/AJ/VCF/PCA_HQ/"
#f_ped<-"~/server/BreastCancer/AJ/Result/AJ_case_control.ped"
#suff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.HQ.gz.plus.HapMap"
cont<-read.table("Resources/1000Gneome_2013/source/integrated_call_samples_v3.20130502.ALL.panel",header=1,fill = T,strip.white = T,sep="\t",check.names=F,stringsAsFactors = F)



pop_classify<-function(nm_sims,pheno,evec){
  print(length(nm_sims))
  loss=0;
  plot(evec$V2,evec$V3,type='n',xlab="",ylab="")
  for(i in 1:dim(pheno)[1]){
      index<-grep(pheno$ID[i],nm_sims)
      if(length(index)>1){index<-which(nm_sims==pheno$ID[i])}
      if(length(index)>0){
        if(evec$V3[index]< -0.03 && evec$V2[index] <0 && evec$V2[index] > -0.01){
          pheno$pop[i]="EAS"
        #  print(paste(  pheno$pop[i],  pheno$sampleID[i]))
          points(evec$V2[index],evec$V3[index],col=1,pch=3)
          text(evec$V2[index],evec$V3[index],labels = "EAS",col=1)
        }else if(evec$V3[index]> 0 && evec$V2[index]> 0.004){
          pheno$pop[i]="EUR"
       #   print(paste(  pheno$pop[i],  pheno$sampleID[i]))
          points(evec$V2[index],evec$V3[index],col=2,pch=4)
          text(evec$V2[index],evec$V3[index],labels = "EUR",col=2)
        }else if(evec$V2[index]< -0.015 && evec$V3[index]>0){
          pheno$pop[i]="AFR"
      #    print(paste(  pheno$pop[i],  pheno$sampleID[i]))
          points(evec$V2[index],evec$V3[index],col=3,pch=5)
          text(evec$V2[index],evec$V3[index],labels = "AFR",col=3)
        }else if (evec$V4[index] < -0.02 && evec$V3[index] >-0.02 && evec$V3[index] <0 ){
          pheno$pop[i]="SAS"
       #   print(paste(  pheno$pop[i],  pheno$sampleID[i]))
          points(evec$V2[index],evec$V3[index],col=4,pch=6)
          
          text(evec$V2[index],evec$V3[index],labels = "SAS",col=4)
          
        }else if(evec$V4[index] > -0.02 && evec$V2[index] > -0.02 && evec$V2[index] < 0.004 && evec$V3[index]>-0.03){
          pheno$pop[i]="Dominican"
       #   print(paste(  pheno$pop[i],  pheno$sampleID[i]))
          points(evec$V2[index],evec$V3[index],col=5,pch=20)
          
          text(evec$V2[index],evec$V3[index],labels = "DOM",col=5)
          
        }else {
          pheno$pop[i]="Hispanic"
       #   print(paste(  pheno$pop[i],  pheno$sampleID[i]))
          points(evec$V2[index],evec$V3[index],col=6,pch=20)
          
          text(evec$V2[index],evec$V3[index],labels = "Hispanic",col=6)
        }
      }else{
       # print (paste("Error: ",pheno$ID[i]," cannot find in evec!"))
        loss=loss+1
        pheno$pop[i]="-"
      }
  }
  print(paste(loss," samples cannot find in Evec!"))
  return(pheno)
}



callPCA<-function(evec,pheno,nm_sims,pops,x,y,p1,p2){
 
    cols<-c("pink","orange","violet","red","darkblue",colors()[563],colors()[c(56,111)],
          colors()[c(435,430,425,371,417,466,498,551,555,574,373,375,376,382,367)])
  #cols<-c("")
  pchs<-c(3,1,4,1,2,c(1:6),c(1:6))
  plot(evec[,x], evec[,y], 
       xlab=paste("eigenvector",x,"\n",p1, "% of observed genetic variation", sep=""), 
       ylab=paste("eigenvector",y,"\n",p2, "% of observed genetic variation", sep=""), 
       main="PCA",
       col=cols[3],
       pch=pchs[3],cex=1,type='n'
       # ,xlim=xr,ylim=yr
       # ,  xlim=c(0.005,0.015),
       # ylim=c(0.002,0.012)
  )   
  pp<-c()
  lnm<-c()
  i=1
  d=0;
  for(c in unique(cont$super_pop)){
    j=1;
    for(p in unique(cont$pop[which(cont$super_pop==c)])){
      index<-which(pops==p)
      pp<-c(pp,index)
      #  d=d+length(index)
      points(evec[index,x], evec[index,y], 
             col=cols[i],
             pch=pchs[j],
             cex=1.2
      )
      j=j+1;
    }
    i=i+1
    lnm<-c(lnm,c)
  }

  index_pop<-which(nm_sims%in% pheno$ID[which(pheno$cohort==cfeature)])
  points(evec[index_pop,x], evec[index_pop,y],
         col=cols[i],
         pch=3,
         cex=1.3)
  
  legend("left",legend=c(lnm,cfeature),fill=cols,bty='n')
  

  
  #abline(h=-0.005)
  
  i=i+1;
  ph=feature
  index_pop<-which(nm_sims%in% pheno$ID) #[which(pheno$cohort==ph)]
  points(evec[index_pop,x], evec[index_pop,y],
         col=cols[i],
         pch=3,
         cex = 1.3)
  legend("topleft",legend=c(feature),fill=cols[i],bty='n')
 # sox17_samples <- read.table("~/server/PAH/PSAP/SOX17.check.txt",header=1,comment.char = "",check.names = F,stringsAsFactors = F)
  #sox17_IDs <- sox17_samples$ProbandName
  # 
 # i="darkgreen"
#  sox17_st<-unlist(lapply(as.character(sox17_IDs),function(x)(substr(x,nchar(x)-15,nchar(x)))))
#  index_pop<-which(nm_sims %in% sox17_st)
#  points(evec[index_pop,x], evec[index_pop,y],
 #          col=i,
  #         pch=20,
   #        cex=2.3)
  #text(evec[index_pop,x], evec[index_pop,y],labels =nm_sims[index_pop],pos = c(1,3),col="darkgreen",cex=0.75)
  #legend("left",legend=c(),fill=cols[i],bty='n')
  # 
  # i="yellow"
  # index_pop<-which(nm_sims%in% other_abcc8_query)
  # points(evec[index_pop,x], evec[index_pop,y],
  #        col=i,
  #        pch=20,
  #        cex=1.3)
  # legend("bottomleft",legend=c("ABCC8-syn","ABCC8-fun"),fill=c("darkgreen","yellow"),bty='n')
  # text(evec[index_pop,x], evec[index_pop,y],labels = c(1:length(index_pop)),pos = c(1,3),col="orange",cex=0.75)
  # 
}

myplot <- function(folder,suff,fpheno,fpdf,feature,cfeature,fout) {
    f_eval<-paste(folder,suff,".eval",sep="")
    f_evec<-paste(folder,suff,".pca.evec",sep="")
    if(missing(cfeature)){cfeature=c()}
    pheno<-read.table(fped,header=1,stringsAsFactors = F,check.names = F,comment.char = "",fill=T)
    #pheno <- read.table("~/server/PAH/PAH_2017/PAH_affected-excluding.ped",header=F,stringsAsFactors = F,check.names = F)
    #pheno<-read.csv(fpheno,stringsAsFactors = F)
    pheno$proband<-pheno$ID
    pheno$ID<-unlist(lapply(pheno$ID,FUN = function(x) gsub("COL-CHUNG_","",x ) ))
    
    a<-unlist(lapply(as.character(pheno$ID),function(x)(substr(x,nchar(x)-15,nchar(x)))))
    pheno$ID <- a
   # sox17_samples <- read.table("~/server/PAH/PSAP/SOX17.check.txt",header=1,comment.char = "",check.names = F)
   # sox17_IDs <- sox17_samples$ProbandName
    ## add the sox17 samples 
    # a<-unlist(lapply(as.character(syn_abcc8_samples),function(x)(substr(x,nchar(x)-15,nchar(x)))))
    #syn_abcc8_query<-a
    #a<-unlist(lapply(as.character(other_abcc8_sample),function(x)(substr(x,nchar(x)-15,nchar(x)))))
    #other_abcc8_query<-a
    
    eval <- read.table(f_eval)
    evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
    evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
    evec <- read.table(f_evec,stringsAsFactors = F)
    evec<-evec[which(evec$V3 < 0.019 & evec$V5 > -0.25),]
    nm_sims<-unlist(lapply(as.character(evec$V1),function(x)(unlist(strsplit(as.character(x),split = ":")))[1]))
  #  phenotypes = unique(pheno$cohort)
 
    pops<-c()
    inds<-as.vector(mapply(as.character(evec$V1),FUN = function(x) unlist(strsplit(x,":"))[1]))
    for(id in inds){
      index <- which(cont$sample == id)
      if(length(index)<1){
        pops <- c(pops,"InHouse")
      }else{
        pops<-c(pops,cont$pop[which(cont$sample==id)])
      }
    }
    ## pops
    hapmap <- grep("NA",evec$V1,ignore.case = T) #which(grep("HapMap",evec$V1,ignore.case = T))
    local <- grep("NA",evec$V1,ignore.case = T,invert = T) # which(grep("HapMap",evec$V1,ignore.case = T,invert = T))
    a<-unlist(lapply(as.character(evec$V1[hapmap]),function(x)(unlist(strsplit(as.character(x),split = ":")))[1]))
 #   pops<-unlist(lapply(a,function(x) unlist(strsplit(x,".",fixed = T))[2]))
    
    euro<-which(evec$V1>0.01 & evec$V2>0.003)
   # write.table(pheno$sampleID[which(pheno$ID %in% nm_sims[euro] & pheno$cohort == "PAH") ],file = "PAH_european.csv",sep="\t",quote = F,row.names = F,col.names = F)
    
 #   write.table(pheno[which(pheno$ID%in%nm_sims[euro]),c("sampleID","cohort") ],file = "RGN_european.csv",sep="\t",quote = F,row.names = F,col.names = F)
    ####
# 
     pdf(fpdf)
     par(mar = c(5.1, 5.1, 2.1, 2.1))
     callPCA(evec,pheno,nm_sims ,pops,2,3,round(eval[1,1]/sum(eval)*100,digits=2),round(eval[2,1]/sum(eval)*100,digits=2))
     abline(v=c(-0.01,0.004),h=c(-0.03,0),lty=2)
     callPCA(evec,pheno,nm_sims,pops ,3,4,round(eval[2,1]/sum(eval)*100,digits=2),round(eval[3,1]/sum(eval)*100,digits=2))
    # abline(h=0.011,lty=2)
     abline(h= -0.02,lty=2)
     callPCA(evec,pheno,nm_sims ,pops,4,5,round(eval[3,1]/sum(eval)*100,digits=2),round(eval[4,1]/sum(eval)*100,digits=2))
     abline(v=0.013,lty=2)
     abline(v= -0.02, lty=2)
    dev.off()
 #   nm_sims<-unlist(lapply(as.character(evec$V1),function(x)(unlist(strsplit(as.character(x),split = ":")))[1]))
    
    pheno<-pop_classify(nm_sims,pheno,evec)
    write.csv(pheno,fout,row.names = F)
    
}

#args<-commandArgs(TRUE)
folder="PAH/Relatedness/"#args[1];
#suff="RGN_Freeze6.rawvariants.recalibrated.hg19_multianno.mappability.vcf.gz.plus.HapMap"   #"RGN_FREEZE6.annotated.gz.plus.HapMap" #args[2];
suff="PAH_07102017_1KG_PCA.plus.HapMap"
f#pheno="PAH/PAH_2017/Relatedness/RGN_freeze6.csv"
fpdf=paste(folder,"PCA---alll.pdf",sep="")
feature="PAH"
cfeature="Domincan_control"
fout="PAH/Result/Data/PAH_pop.csv"
fped="PAH/Result/Data/source/VCX_Control.ped"
nm_sims<-c()
myplot(folder,suff,fped,fpdf,feature,cfeature,fout  )
