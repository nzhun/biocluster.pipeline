#uff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.gz.plus.HapMap"

#folder<-"~/server/BreastCancer/AJ/VCF/PCA_HQ/"
#f_ped<-"~/server/BreastCancer/AJ/Result/AJ_case_control.ped"
#suff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.HQ.gz.plus.HapMap"

myplot <- function(folder,suff,f_ped) {
    f_eval<-paste(folder,suff,".eval",sep="")
    f_evec<-paste(folder,suff,".pca.evec",sep="")
    ped<-read.table(f_ped,fill=T,check.names=F)
    cases<-unlist(lapply(ped[which(ped$V6==2),1],function(x) substr(x,nchar(as.character(x))-15,nchar(as.character(x)))))
    controls<-unlist(lapply(ped[which(ped$V6==1),1],function(x) substr(x,nchar(as.character(x))-15,nchar(as.character(x)))))
    
   
    
    eval <- read.table(f_eval)
    evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
    evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
    evec <- read.table(f_evec)
  
    case_inh<-which(unlist(lapply(evec$V1,function(x)unlist(strsplit(as.character(x),":"))[1])) %in% cases)
    control_inh<-which(unlist(lapply(evec$V1,function(x)unlist(strsplit(as.character(x),":"))[1])) %in% controls)
    
    cols<-c("red","dark blue","gray","darkgray","orange","cyan","violet",colors()[c(367,371,373,375,376,382,417,435,430,425,466,498,551,555,574)])
    pchs<-c(20,4,3,1,2,c(1:6),c(1:6))
    
    ## pops
    hapmap<-grep("HapMap",evec$V1,ignore.case = T) #which(grep("HapMap",evec$V1,ignore.case = T))
    local<- grep("HapMap",evec$V1,ignore.case = T,invert = T) # which(grep("HapMap",evec$V1,ignore.case = T,invert = T))
    a<-unlist(lapply(as.character(evec$V1[hapmap]),function(x)(unlist(strsplit(as.character(x),split = ":")))[1]))
    pops<-unlist(lapply(a,function(x) unlist(strsplit(x,".",fixed = T))[2]))
    
    
    pdf(paste(folder,"PCA.pdf",sep=""))
    par(mar = c(5.1, 5.1, 2.1, 2.1))
    plot(evec[,2], evec[,3], 
         xlab=paste("eigenvector1\n",evec1.pc, "% of observed genetic variation", sep=""), 
         ylab=paste("eigenvector2\n",evec2.pc, "% of observed genetic variation", sep=""), 
         main="PCA",
         col=cols[3],
         pch=pchs[3],cex=1,type='n'
         # ,xlim=xr,ylim=yr
         # ,  xlim=c(0.005,0.015),
         # ylim=c(0.002,0.012)
    ) 
    
    
  #  points(evec[hapmap,2], evec[hapmap,3], 
 #          col=cols[2],
 #          pch=pchs[2],
#           cex=1.3)
    k=3
    query_pops<-unique(pops)
    for(pop in query_pops){
      index_pop<-hapmap[which(pops==pop)]
      points(evec[index_pop,2], evec[index_pop,3], 
             col=cols[k],
             pch=pchs[k],
             cex=1.3)
      k=k+1;
    }
  #  points(evec[local,2], evec[local,3], 
  #         col=cols[1],
  #         pch=pchs[1],
 #          cex=1.2)
    points(evec[case_inh,2], evec[case_inh,3], 
           col=cols[1],
           pch=pchs[1],
           cex=1.2
    )
    
    
    
    points(evec[control_inh,2], evec[control_inh,3], 
           col=cols[2],
           pch=pchs[2],
           cex=1.2
    )
    legend("bottomleft",legend=c("Case","Control",query_pops),col=cols,pch=pchs,bty='n')
    #, cex=2.5, cex.axis=1.5,cex.lab=1.5)
    
   # dev.off()
    
    
   # pdf(paste(folder,"PCA-Inhouse.pdf",sep=""))
    #par(mar = c(5.1, 5.1, 2.1, 2.1))
    xr=c(min(evec[c(case_inh,control_inh),2])*0.7,max(evec[c(case_inh,control_inh),c(2)])*1.5)
    yr=c(min(evec[c(case_inh,control_inh),c(3)])*0.7,max(evec[c(case_inh,control_inh),c(3)])*1.5)
    plot(evec[,2], evec[,3], 
         xlab=paste("eigenvector1\n",evec1.pc, "% of observed genetic variation", sep=""), 
         ylab=paste("eigenvector2\n",evec2.pc, "% of observed genetic variation", sep=""), 
         main="PCA",
         col=cols[3],
         pch=pchs[3],cex=1,type='n'
          ,xlim=xr,ylim=yr
         # ,  xlim=c(0.005,0.015),
         # ylim=c(0.002,0.012)
    ) 
    
    
    #  points(evec[hapmap,2], evec[hapmap,3], 
    #          col=cols[2],
    #          pch=pchs[2],
    #           cex=1.3)
    k=3
    query_pops<-unique(pops)
    for(pop in query_pops){
      index_pop<-hapmap[which(pops==pop)]
      points(evec[index_pop,2], evec[index_pop,3], 
             col=cols[k],
             pch=pchs[k],
             cex=1.3)
      k=k+1;
    }
    #  points(evec[local,2], evec[local,3], 
    #         col=cols[1],
    #         pch=pchs[1],
    #          cex=1.2)
 
    
    
    
    points(evec[control_inh,2], evec[control_inh,3], 
           col=cols[2],
           pch=pchs[2],
           cex=1.2
    )
    
    points(evec[case_inh,2], evec[case_inh,3], 
           col=cols[1],
           pch=pchs[1],
           cex=1.2
    )
    legend("bottomleft",legend=c("Case","Control",query_pops),col=cols,pch=pchs,bty='n')
    #, cex=2.5, cex.axis=1.5,cex.lab=1.5)
    
    dev.off()
}

args<-commandArgs(TRUE)
#args<-c("~/server/PAH/VCF/PCA/","PAH_63_PCGC_parent00.plus.HapMap","~/server/PAH/Result/PAH_63_PCGC.ped")
folder=args[1];
suff=args[2];
ped=args[3];
myplot(folder,suff,ped)
