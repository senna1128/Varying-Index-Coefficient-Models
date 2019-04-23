##############################
## Preparing genetic datasets
rm(list=ls())
# change working directory
 setwd("./Mladen/project/pro1/SIVC/4.5code/realdata/geneticdata")
# read data
dat_Pop1_SNP = read.csv("./data/Pop_Int_SNPdata.csv")
dat_Pop2_SNP = read.csv("./data/Pop_Prec_SNPdata.csv")
dat_Pop1_resp = read.csv("./data/Pop_Int_PHENOdata.csv")
dat_Pop2_resp = read.csv("./data/Pop_Prec_PHENOdata.csv")

## make observations independent in each population
N_Pop_1 = dim(dat_Pop1_SNP)[1]
halfN_Pop1 = round(N_Pop_1/2)
dat_Pop1 = dat_Pop1_resp[c(1:halfN_Pop1, (N_Pop_1+halfN_Pop1+1):(2*N_Pop_1)),]

N_Pop_2 = dim(dat_Pop2_SNP)[1]
halfN_Pop2 = round(N_Pop_2/2)
dat_Pop2 = dat_Pop2_resp[c(1:halfN_Pop2, (N_Pop_2+halfN_Pop2+1):(2*N_Pop_2)),]

## select common SNPs
Pop1_SNPnames = colnames(dat_Pop1_SNP[-1])
Pop2_SNPnames = colnames(dat_Pop2_SNP[-1])
Pop_SNPnames = intersect(Pop1_SNPnames, Pop2_SNPnames)

## construct new dataset based on common SNPs
dat_Pop1 = cbind(dat_Pop1, dat_Pop1_SNP[,Pop_SNPnames])
dat_Pop2 = cbind(dat_Pop2, dat_Pop2_SNP[,Pop_SNPnames])


## Replace population groups by Gaussian distribution with different means
set.seed(2019)
dat_Pop1$cont.pop = rnorm(N_Pop_1,0,1)
dat_Pop2$cont.pop = rnorm(N_Pop_2,50,1)

## Replace evaluation locations by t_13 distribution with different means
dat_Pop1$cont.env = c(rt(halfN_Pop1,13), 50+rt(N_Pop_1 - halfN_Pop1,13))
dat_Pop2$cont.env = c(rt(halfN_Pop2,13), 50+rt(N_Pop_2 - halfN_Pop2,13))

## combine two dataset and rearrange columns
dat = rbind(dat_Pop1, dat_Pop2)
newdat = data.frame(X = dat$X, production = dat$production, rust = dat$rust, green = dat$green,
                   Pop = dat$cont.pop, Env = dat$cont.env, dat[,Pop_SNPnames])
write.csv(newdat,'GenDat.csv')


# draw histograms for responses
jpeg("Figures/prodhist.jpg")
hist(newdat$production, xlab="production", freq=FALSE, ylab="density",breaks=20,main = "")
dev.off()

jpeg("Figures/rusthist.jpg")
hist(newdat$rust, xlab="rust", freq=FALSE, ylab="density",breaks=20,main = "")
dev.off()

jpeg("Figures/greenhist.jpg")
hist(newdat$green, xlab="green", freq=FALSE, ylab="density",breaks=20,main = "")
dev.off()





