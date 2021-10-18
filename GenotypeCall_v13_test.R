### =========================================================================
### GenotypeCall_v13.R
### Genotype Calling
### -------------------------------------------------------------------------
### by Jia: 10/8/2021
### Routines for calling and manipulating genotypes.
###
### This script is used to compare results with EXCEL examples

#BiocManager::install("VariantAnnotation")
# pacman::p_load(VariantAnnotation, stringr,tibble, dplyr, openxlsx, pbapply,readr,vcfR,hash,
#                data.table,progress,matrixStats,pastecs,DescTools,ggplot2)
pacman::p_load(VariantAnnotation, stringr, stringi, tibble, dplyr, data.table, tidyr)
rm(list = ls(all = TRUE))
###################################################################################################
## Set program parameters
###################################################################################################
p.error <- 0.05 #A priori estimate of read error
mLR <- 0
p.error.test <- 0.1 #Error value for binomial test of called genotype
###################################################################################################
## Set the output directory, create a log file, and print the parameters
###################################################################################################
setwd("")

out.dir <- ""
version <- "v13" #Same as in the file name 
fname <- paste0("GenotypeCall_test_", version)

sink(); #Remove existing connections - Do not worry if you get a warning message
sink.out <- file(paste0(out.dir, "\\", fname,".log"), open="wt")
sink(sink.out, type = "output", split=TRUE)
###################################################################################################
## Construct likelihood matrix for Bayesian method
###################################################################################################
#Construct likelihood matrix. pD_G is the probability of D vector of reads given the genotype (G)
#Matrix of probabilities for obtaining the indicated read (R) given the homozygous (ho) genotype
pR_ho <-matrix(p.error/3,4,4)
diag(pR_ho) <- (1-p.error)

#Matrix of probabilities for obtaining the indicated read (R) given the homozygous (he) genotype
#Example for AC heterozygote
#   Some A alleles were read correctly = (1-p.error)/2
#   Some A alleles are incorrectly read C alleles = (p.error/3)/2
pR_he <- matrix(p.error/3, nrow=4, ncol=6)
pR_he <- mapply(function(x,y) replace(x, y, 0.5*(1-p.error)+p.error/6), unname(split(pR_he,rep(1:ncol(pR_he)))),combn(1:4,2,simplify=F))

#Matrix of probabilities for obtaining the indicated read (R) given the genotype
pR_G_m <- cbind(pR_ho, pR_he)
#prepare conversion table for numeric GT
p_conv <- list(1,2,3,4,c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4),1,2,3,4)
names(p_conv) <- 1:14
a_conv <- c("AA","CC","GG","TT","AC","AG","AT","CG","CT","GT","A","C","G","T")
names(a_conv) <- 1:14
###################################################################################################
## Use the function call_GT to compute GT and binomial test
###################################################################################################
#=====Set various variables to test callGT ==========#
pG_dip <- c(0.32,0.68,0,0,0,0,0,0,0,0) #for test purpose
pG_dip <- rep(0.1,10) #for test purpose
#pG_dip <- c(0.0625, 0.0625, 0.0625, 0.0625, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125) #for test purpose
#pG_dip <- c(0.000,	0.250,	0.000,	0.250,	0.000,	0.000,	0.000,	0.000,	0.500,	0.000) #for test purpose
D <- c(16,3,1,1) #for test purpose
D.sum <- sum(D) #Total number of reads
#check if read depth==0
#=====compute GT ==========#
#if(type=="d") { #Change pR_G  for haploid calls
pR_G <- pR_G_m
#convert calculated likelihoods to ln
pd_G <- dbinom(D, D.sum, pR_G, log=T) #Probability of d(k)|G (ln)
pD_G <- colSums(pd_G) #P(D|Gj) = MULT [P(dk|Gj)] (ln)
pGxpD_G <- pD_G + log(pG_dip) #P(Gj)P(D|Gj) (ln)
pGxpD_G.sum <- ifelse(log(sum(exp(pGxpD_G)))==-Inf, max(pGxpD_G), log(sum(exp(pGxpD_G)))) #SUM[P(Gj)P(D|Gj)]
pG_D <- pGxpD_G-pGxpD_G.sum #P(Gi|D) (ln)

#likelihood ratio (ln transformed)
#The next statement handles the case where the top two calls have equal max LR
log_LR <- sapply(seq_along(pG_D), function(x){pG_D[x]-max(pG_D[-x])})
#log_LR <- sapply(pG_D,function(x){x-max(pG_D[pG_D!=x])})

#FPR (false positive risk)
FPR <- sum(sapply(pG_D,exp))-exp(pG_D)

#find the genotype with the highest ln = called genotype
#1:10={c("AA","CC","GG","TT","AC","AG","AT","CG","CT","GT")}
call <- ifelse(max(log_LR) <= mLR, NA,
               c(1:10)[which.max(log_LR)])
#Calculate the multinomial probability of this genotype
pCall <- ifelse(max(log_LR) <= mLR, NA, exp(pG_D[which.max(log_LR)]))

#=====compute Binomial test stats ==========#
if(call%in%5:10){
  #Calculate the number of reads for two called alleles and combined error alleles
  allele.num <- D[p_conv[[call]]] #Number of first called heterozygous alleles
  allele.err.num <- D.sum - as.numeric(allele.num[1] + allele.num[2]) #Sum of non-called alleles
  #Calculate the error percentage for this call
  err.perc <- allele.err.num/D.sum
  #Calculate the probability that the error percentage is less than p.error.test
  pErr <- binom.test(allele.err.num, D.sum, p.error.test, alternative="greater")$p.value
  #Calculate the probability that there's an equal distribution of reads (50/50) between two alleles
  pHet <- as.numeric(binom.test(min(allele.num[1], allele.num[2]), allele.num[1] + allele.num[2], 0.50, alternative="less")$p.value)
} else if (!is.na(call)) {
  #Calculate the number of reads for two most frequent alleles and combined error alleles
  allele1.num <- D[p_conv[[call]]] #Number of called homozygous alleles
  allele2.num <- sort(as.vector(D), decreasing=TRUE)[2] #Number of second most frequent alleles
  allele.err.num <- D.sum - as.numeric(allele1.num)
  #Calculate the error percentage for this call
  err.perc <- allele.err.num/D.sum
  #Calculate the probability that the error percentage is less than p.error.test
  pErr <- binom.test(allele.err.num, D.sum, p.error.test, alternative="greater")$p.value
  #Calculate the probability that there's an equal distribution of reads (50/50) between two alleles
  #pHet <- as.numeric(binom.test(min(allele1.num, allele2.num), allele1.num + allele2.num, 0.50, alternative="less")$p.value)
  pHet <- NA
}
#Find the minimum FPR
FPR <- FPR[which.max(log_LR)]
#Find the maximum lnLR and LR
log_LR <- max(log_LR)
LR <- exp(log_LR)
df <- c(call=call, pCall=pCall, err.perc=err.perc, pErr=pErr,pHet=pHet,LR=LR,log_LR=log_LR,FPR=FPR)
df[1] <- a_conv[call]

df
