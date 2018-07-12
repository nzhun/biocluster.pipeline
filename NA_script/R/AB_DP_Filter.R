#setwd('/Volumes/pe/alh2194/mosaic/')
library("bbmle")
library("emdbook") 
library("SuppDists")
library("ggplot2")
library("grid")
library("gridExtra")
#library(foreach)
#library(doParallel)

## estimate Power and FDR at different LR/BAF+DP levels
## For each LR cutoff, plot Power as function of 1-BAF vs. DP
## Take weighted average of Power across DP values (weighted by P(DP=x) based on DoC?), for each BAF?, for each LR cutoff?
findBAFcut <- function(N,  lcut, theta = 50) {
  for (af in c(100:1) / 200 ) {
    ad = round(N * af)
    lr = dbetabinom(ad, N, p=af, theta = theta) / dbetabinom(ad, N, p=0.5, theta = theta)
    if (lr >= lcut) {
      return(af)
    }
  }
  return(0)
}


# given candidate Nalt and DP, calculate likelihood ratios for M0, M1, Mx
# use theta=mtheta or thetahat, since it only depends on capture kit
llratio <- function(alt, dp, th){
  phat <- alt/dp
  px <- dbetabinom(alt, prob=phat, size=dp, theta=th) # P(D|Mx)
  p1 <- dbetabinom(alt, prob=0.5, size=dp, theta=th) # P(D|M1)
  p0 <- dbetabinom(alt, prob=0, size=dp, theta=th) # P(D|M0)
  x_1 <- px/p1
  x_0 <- px/p0 # irrelevant, since p0 = 0
  return(x_1)
}

betaBinomEst <- function(counts, total) {
  mtmp <- function(prob,size,theta) { -sum(dbetabinom(counts,prob,size,theta,log=TRUE)) } 
  m0 <- mle2(mtmp,start=list(prob=0.5,theta=10),data=list(size=total))
  # MLE of theta
  t = coef(m0)["theta"]
  return(t)
}

fitReadCounts <- function(a,fplot) {
  ## main function
#  colnames(a)[1:23] = c('id','gene','chr','pos','ref','alt','filt','refdp','altdp','aachange','exac_freq','exac_popfreq','vclass','pp2hdiv','pp2hvar', 'metasvm', 'metasvm_flag', 'cadd', 'cadd_flag', 'mcap', 'mcap_flag', 'adhoc_flag','src') 
 a$refdp<-a[,2]-a[,1]  
 names(a)<-c("altdp","N","refdp")
  
  a$refdp <- as.numeric(a$refdp)
  a$altdp <- as.numeric(a$altdp)
  x<-a
#  N = a$refdp + a$altdp 
  #x = cbind.data.frame(a, N)
 # names(x) = c('id','gene','chr','pos','ref','alt','filt','refdp','altdp','aachange','exac_freq','exac_popfreq','vclass','pp2hdiv','pp2hvar', 'metasvm', 'metasvm_flag', 'cadd', 'cadd_flag', 'mcap', 'mcap_flag', 'adhoc_flag','src','N') 
  
  x<-x[!is.na(x$N),]
  x <- x[x$refdp>0&x$altdp>0,] # ignore 0 ref or alt DP entries
  
  # generate depth table
  results = estimateTheta(x)  
  results = results[results$m  > 0, ]
  
  # theta and p estimations
  # take the median value of estimation from the parts with large enough sample size
  mp = median(results[results$N > 15&results$N < 100   & results$m > 20, ]$p)
  mtheta =  median(results[results$N > 15 &results$N < 100  & results$m > 20, ]$theta)  # alex: min(50,median(results[results$N > 12 & results$m > 20, ]$theta)) # cap mtheta to not be too large
  # take weighted average (weighted by # variants supporting theta estimation at each depth, m)
  thetahat = sum(results[results$N > 15 & results$N < 100 &
                 results$m > 20, ]$m*results[results$N > 15 & results$N <100 & 
                 results$m > 20, ]$theta)/sum(results[results$N > 15 &results$N <100 & results$m > 20, ]$m)  
  x$mp=mp
  x$mtheta=mtheta
  x$thetahat=thetahat
  # alex: min(50,sum(results[results$N > 12 & results$m > 20, ]$m*results[results$N > 12 & results$m > 20, ]$theta)/sum(results[results$N > 12 & results$m > 20, ]$m))
  
  # ssc estimations
  #mp = median(results$p) # relax for ssc set
  #mtheta = min(100,median(results$theta)) # cap mtheta to not be too large
  #thetahat = min(100,sum(results$m*results$theta)/sum(results$m))
  #mp=0.47
  #thetahat=50
  
  # calculate p values for candidates
  pvalues <- testBetaBinomial(x, mp, thetahat)
  #  generate candidates
  x = cbind.data.frame(x, pvalues)
  lrcut = 100
  x$lr <- mapply(llratio, x$altdp, x$N, thetahat)
  x$post <- x$lr *0.04 # posterior odds = LR * p(0/1*)/p(0/1)
  x$col = 'black' # for later plots
  x$col[x$lr>lrcut & x$altdp/x$N<=0.35]='red'
  
 
  #bafcut = 1-max(x$altdp/x$N)
  z = x[x$lr > lrcut & x$altdp/x$N<=0.35,]
  #y = x[x$altdp/x$N < bafcut & x$lr<=lrcut,]
  #z.1 = x[x$altdp/x$N < 1-max(x$altdp/x$N) & x$lr>max(x[x$altdp/x$N > 0.5,]$lr),]
 pdf(fplot) 
  plot(x$altdp/x$N, log10(x$lr), xlab="BAF", ylab="log10 LR", main="BAF vs. LR", cex=0.5, xlim=c(0,1))
  lines(z$altdp/z$N, log10(z$lr), type="p", col="red", cex=0.5)
  #lines(y$altdp/y$N, log10(y$lr), type="p", col="blue", cex=0.5)
  #lines(z.1$altdp/z.1$N, log10(z.1$lr), type="p", col="purple", cex=0.5)
  #abline(v=1-max(x$altdp/x$N), col="blue")
  #abline(h=log10(max(x[x$altdp/x$N > 0.5,]$lr)), col="red")
  abline(h=log10(lrcut), col="red")
  abline(v=0.35, col="red")
  legend("topright", c("all variants", "candidate mosaics"), pch=1, col=c("black", "red"), cex=0.75)
  
  ## Yufeng's DP vs alt fraction plot
  #x$lr <- mapply(llratio, x$alt, x$N, thetahat)
  
  # smoothScatter(a$ref + a$alt, a$alt/ (a$alt+a$ref),  xlab = "DP", ylab="alt fraction")
  mn =max(x$N) # min(max(x$N), 800)
  ci = matrix(0, nrow = mn, ncol=3)
  for (i in 10:mn) {
    lci = 0
    hci = mn 
    d = 0
    for (j in 1:i) {
      d = sum(dbetabinom(0:j, prob = mp, size = i, theta = thetahat  ) )
      if (d > 0.025) {
        lci = j - 1
        break
      } 
    }   
    for (j in i:1) {
      d = sum(dbetabinom(i:j, prob = mp, size = i, theta = thetahat  ) )
      if (d > 0.025) {
        hci = j + 1
        break
      } 
    }
    ci[i, 1] = lci
    ci[i, 2] = hci
    ci[i, 3] = i
  }
 
  x$FP<-0
  for(i in 10:dim(ci)[1]){
    index<-which(x$N==ci[i,3] & (x$altdp< ci[i,1] | x$altdp > ci[i,2]) )
    x$FP[index]=1
  }
  index<-which(x$FP==1)
  
  plot(x$refdp + x$altdp, x$altdp/ (x$altdp+x$refdp),  xlab = "DP", ylab="BAF", pch=20, ylim=range(c(0.05, 0.95)), col=x$col)
  points(x$N[index], x$altdp[index]/ (x$N[index]),col="red",pch=3)
  lines(1:mn, ci[,1]/ci[,3], col='red')
  lines(1:mn, ci[,2]/ci[,3], col='red')
  abline(h=mp, col='blue')
  title(paste('DP vs. BAF \n ',paste("Theta",format(thetahat,digits = 2,trim = T),sep=":")))
  legend("topright", c("all variants", "candidate mosaic"), col=c("black", "red"), pch=1, pt.cex=0.5, cex=0.75)
  
  plot(x$refdp + x$altdp, x$altdp/ (x$altdp+x$refdp),  xlab = "DP", ylab="BAF", pch=20,xlim=c(0,200), ylim=range(c(0.05, 0.95)), col=x$col)
  points(x$N[index], x$altdp[index]/ (x$N[index]),col="red",pch=3)
  lines(1:mn, ci[,1]/ci[,3], col='red')
  lines(1:mn, ci[,2]/ci[,3], col='red')
  abline(h=mp, col='blue')
  title(paste('DP vs. BAF \n ',paste("Theta",format(thetahat,digits = 2,trim = T),sep=":")))
  legend("topright", c("all variants", "candidate mosaic"), col=c("black", "red"), pch=1, pt.cex=0.5, cex=0.75)
  
  
  ## yufeng's overdispersion plot
  ## plot dp vs. var(nalt)
  plot(results$N, results$var,  xlab = "DP", ylab = "Var(BAF)", log="y", cex=0.5,
       main=paste("DP vs. Var(BAF) \n Binomial vs. Beta-Binomial\n",
                  paste("Theta",format(thetahat,digits = 2,trim = T))))
  ## if binomial
  lines(results$N, results$N * results$p * (1 - results$p), col='blue')
  ## if beta-binomial, with estimated p and theta
  lines(results$N, results$N * mp * (1 - mp) * (results$N + thetahat) / (thetahat + 1), col='red')
  legend("bottomright", c("Binomial", "Beta-Binomial"), col=c("blue", "red"), lty=1)
  
  
  plot(results$N, results$var,  xlab = "DP", ylab = "Var(BAC)", log="y", cex=0.5,xlim=c(0,200),
       main=paste("DP vs. Var(BAF) \n Binomial vs. Beta-Binomial\n",
                  paste("Theta",format(thetahat,digits = 2,trim = T))))
  ## if binomial
  lines(results$N, results$N * results$p * (1 - results$p), col='blue')
  ## if beta-binomial, with estimated p and theta
  lines(results$N, results$N * mp * (1 - mp) * (results$N + thetahat) / (thetahat + 1), col='red')
  legend("bottomright", c("Binomial", "Beta-Binomial"), col=c("blue", "red"), lty=1)
  
dev.off()
  #priv <- x[x$exac_freq==0,]
  
#  write.table(z, "parents.0216.candidates.txt", quote=F,row.names=F, sep="\t")
 # write.table(x, "parents.0216.denovo.txt", quote=F,row.names=F, sep="\t")
  #write.table(priv, "parents.0131.private.txt", quote=F, row.names=F, sep="\t")
  return(x)
  #return(results)
}

estimateTheta <- function(x) {
  mn = min(max(x$N), 500)
  results = matrix(0, nrow=length(levels(factor(x$N))) , ncol=6)
  results = as.data.frame(results)
  colnames(results) = c("N", "m", "p", "var", "theta","nalt")
  i = 1
  for (t in 1:mn) {
    v = x[x$N == t, ]$altdp
    if (length(v) > 3 ) { 
      theta = betaBinomEst(v, t)[1]
      results[i,1] = t
      results[i,2] = length(v)
      results[i,3] = mean(v)/t
      results[i,4] = var(v)
      results[i,5] = theta	
      results[i,6] = mean(v) # mean of nalts of variants used to support this
      
      i = i + 1
    }
  }
  return(results)
}

testBetaBinomial <- function(x, mp , th) {
  ### now given p and theta, test each variant against a null. The alternative hypothesis is that the variant is a mosaic
  pvalues = rep(0, nrow(x))
  for ( i in 1:nrow(x)) {
    p = sum(  dbetabinom(1:x[i,]$alt,   prob= mp, size = x[i,]$N,  theta = th ) )
    pvalues[i] = p
  }
  return(pvalues) 
}


#a <- read.table('ADfile.proband.indels.txt', sep='\t', header=T)
#fplot="AB_DP.pdf"
#results = fitReadCounts(a,fplot)
