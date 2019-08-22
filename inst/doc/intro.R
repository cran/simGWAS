## ------------------------------------------------------------------------
library(simGWAS)

## ------------------------------------------------------------------------
set.seed(42)
nsnps <- 100
nhaps <- 1000
lag <- 5 # genotypes are correlated between neighbouring variants
maf <- runif(nsnps+lag,0.05,0.5) # common SNPs
laghaps <- do.call("cbind", lapply(maf, function(f) rbinom(nhaps,1,f)))
haps <- laghaps[,1:nsnps]
for(j in 1:lag) 
    haps <- haps + laghaps[,(1:nsnps)+j]
haps <- round(haps/matrix(apply(haps,2,max),nhaps,nsnps,byrow=TRUE))

snps <- colnames(haps) <- paste0("s",1:nsnps)
freq <- as.data.frame(haps+1)
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)

## ------------------------------------------------------------------------
CV=sample(snps[which(colMeans(haps)>0.1)],2)
g1 <- c(1.4,1.2)

## ------------------------------------------------------------------------
FP <- make_GenoProbList(snps=snps,W=CV,freq=freq)
zexp <- expected_z_score(N0=10000, # number of controls
              N1=10000, # number of cases
              snps=snps, # column names in freq of SNPs for which Z scores should be generated
              W=CV, # causal variants, subset of snps
              gamma.W=log(g1), # odds ratios
              freq=freq, # reference haplotypes
              GenoProbList=FP) # FP above

## ------------------------------------------------------------------------
plot(1:nsnps,zexp); abline(v=which(snps %in% CV),col="red"); abline(h=0)

## ------------------------------------------------------------------------
zsim <- simulated_z_score(N0=10000, # number of controls
              N1=10000, # number of cases
              snps=snps, # column names in freq of SNPs for which Z scores should be generated
              W=CV, # causal variants, subset of snps
              gamma.W=log(g1), # log odds ratios
              freq=freq, # reference haplotypes
          nrep=3)

## ---- fig.height = 6-----------------------------------------------------
op <- par(mfcol=c(3,1))
for(i in 1:3) {
  plot(1:nsnps,zexp,xlab="SNP",ylab="Z score"); abline(v=which(snps %in% CV),col="red"); abline(h=0)
  title(main=paste("Replication",i))
  points(1:nsnps,zsim[i,],col="blue",pch=2)
} 
par(op)

## ------------------------------------------------------------------------
 vbetasim <- simulated_vbeta(N0=10000, # number of controls
              N1=10000, # number of cases
              snps=snps, # column names in freq of SNPs for which Z scores should be generated
              W=CV, # causal variants, subset of snps
              gamma.W=log(g1), # log odds ratios
              freq=freq, # reference haplotypes
          nrep=3)
  betasim <- zsim / sqrt(vbetasim)

## ------------------------------------------------------------------------
CV
log(g1)
betasim[,c(which(snps==CV[1]),which(snps==CV[2]))]

## ------------------------------------------------------------------------

plot(1:nsnps,betasim[1,],
     xlab="SNP",ylab="Beta (log OR)",
     ylim=c(min(betasim),max(betasim)))
abline(v=which(snps %in% CV),col="red"); abline(h=0)
for(i in 2:3) {
  points(1:nsnps,betasim[i,],col=i,pch=i)
} 
  

