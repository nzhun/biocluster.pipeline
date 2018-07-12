folder="~/server/PAH/PAH_2017/Relatedness/"
suff<-"RGN_FREEZE6.annotated.gz.plus.HapMap"
f_ped<-"~/server/PAH/PAH_2017/PAH.ped"
#folder<-"~/server/BreastCancer/AJ/VCF/PCA_HQ/"
#f_ped<-"~/server/BreastCancer/AJ/Result/AJ_case_control.ped"
#suff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.HQ.gz.plus.HapMap"


f_eval<-paste(folder,suff,".eval",sep="")
f_evec<-paste(folder,suff,".pca.evec",sep="")
ped<-read.table(f_ped)

cases<-unlist(lapply(ped[which(ped$V6==2),1],function(x) substr(x,nchar(as.character(x))-15,nchar(as.character(x)))))
controls<-unlist(lapply(ped[which(ped$V6==1),1],function(x) substr(x,nchar(as.character(x))-15,nchar(as.character(x)))))
eval <- read.table(f_eval)
evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
evec <- read.table(f_evec)
par(mar = c(6.1, 6.1, 2.1, 2.1))

cols<-c("red","dark blue","gray","orange")
pchs<-c(20,4,3,1)
inhouse<-which(evec$V12=="Case")
case_inh<-which(unlist(lapply(evec$V1,function(x)unlist(strsplit(as.character(x),":"))[1])) %in% cases)
xr=c(min(evec[inhouse,2]),max(evec[inhouse,c(2)]))
yr=c(min(evec[inhouse,c(3)]),max(evec[inhouse,c(3)]))

pdf(paste(folder,"PCA.pdf",sep=""))
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

hapmap<-which(evec$V12!="Case")

points(evec[hapmap,2], evec[hapmap,3], 
       col=cols[4],
       pch=pchs[4],
       cex=1.3)
control_inh<-which(unlist(lapply(evec$V1,function(x)unlist(strsplit(as.character(x),":"))[1])) %in% controls)
inh<-inhouse[which(!inhouse %in%c(case_inh,control_inh))]
points(evec[inh,2], evec[inh,3], 
       col=cols[3],
       pch=pchs[3],
       cex=1.2)

legend("left",legend=c("Case","Control","Inhouse-NonGWAS","HAPMAP"),col=cols,pch=pchs,bty='n')
#, cex=2.5, cex.axis=1.5,cex.lab=1.5)
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



############# enlarge ###########
plot(evec[,2], evec[,3], 
     xlab=paste("eigenvector1\n",evec1.pc, "% of observed genetic variation", sep=""), 
     ylab=paste("eigenvector2\n",evec2.pc, "% of observed genetic variation", sep=""), 
     main="PCA",
     col=cols[3],
     pch=pchs[3],cex=1,type='n',
    # ,xlim=xr,ylim=yr
      xlim=c(0.005,0.015),
     ylim=c(0.002,0.012)
    ) 

hapmap<-which(evec$V12!="Case")

points(evec[hapmap,2], evec[hapmap,3], 
       col=cols[4],
       pch=pchs[4],
       cex=1.3)

legend("topleft",legend=c("Case","Control","Inhouse-NonGWAS","HAPMAP"),col=cols,pch=pchs,bty='n')
   #, cex=2.5, cex.axis=1.5,cex.lab=1.5)

inh<-inhouse[which(!inhouse %in%c(case_inh,control_inh))]
points(evec[inh,2], evec[inh,3], 
       col=cols[3],
       pch=pchs[3],
       cex=1.2)


points(evec[case_inh,2], evec[case_inh,3], 
       col=cols[1],
       pch=pchs[1],
       cex=1.2
)


control_inh<-which(unlist(lapply(evec$V1,function(x)unlist(strsplit(as.character(x),":"))[1])) %in% controls)

points(evec[control_inh,2], evec[control_inh,3], 
       col=cols[2],
       pch=pchs[2],
       cex=1.2
       )


dev.off()
# 
# 
# library("MASS")
# #dist<-read.table("~/server/BreastCancer/AJ/VCF/Test/similarity_matrix_hamming.txt")
# #fit <- isoMDS(as.dist(dist), k=2) # k is the number of dim
# #fit # view results
# 
# # plot solution 
# cols<-c("red","blue","gray")
# pchs<-c(20,3,4)
# x <- fit$points[,1]
# y <- fit$points[,2]
# 
# pdf("~/server/BreastCancer/AJ/VCF/Test/Genetic_Hamming_distance.pdf")
# 
# case_ids<-which(unlist(lapply(names(dist),function(x) substr(x,nchar(as.character(x))-15,nchar(as.character(x)))))%in%cases)
# control_ids<-which(unlist(lapply(names(dist),function(x) substr(x,nchar(as.character(x))-15,nchar(as.character(x)))))%in%controls)
# inhouse<-which(! c(1:length(x)) %in%  c(case_ids,control_ids))
# xr=c(min(x[c(case_ids,control_ids)])-0.005,max(x[c(case_ids,control_ids)]))
# yr=c(min(y[c(case_ids,control_ids)]),max(y[c(case_ids,control_ids)]))
# 
# plot(x, y,col=cols[3], xlab="Coordinate 1", ylab="Coordinate 2", 
#      main="Genetic Hamming Distance", pch=pchs[3],xlim=xr,ylim=yr)
# points(x[case_ids], y[case_ids],pch=pchs[1],col=cols[1])
# 
# points(x[control_ids], y[control_ids],pch=pchs[2],col=cols[2])
# 
# legend("topleft",legend =c("case","control","Inhouse-nonGwas"),col=cols,bty='n',pch=pchs )
# #text(x, y, labels = row.names(mydata), cex=.7)
# dev.off()
