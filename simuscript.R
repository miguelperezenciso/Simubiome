################################################################################################
# Miguel Perez-Enciso (miguel.perez@uab.es)
# otu data from Difford
# pruned snp data from Dallace et al gen.txt.gz, 32k snps
################################################################################################

source('simubiome.R')

# folder containing snp and otu data
ddata = 'data/'

# read bacteria
B = read.biome(paste0(ddata,'journal.pgen.1007580.s005.txt.gz'))
# read archaea (not used)
A = read.biome(paste0(ddata,'journal.pgen.1007580.s006.txt.gz'))
# read genotypes
X = read.gen(paste0(ddata,'gen.txt.gz'))
# N SNPs
Nsnp = nrow(X)
# N individuals
Nind = ncol(X)
# N OTUs
Notu = nrow(B)
try (if (Nind != ncol(B)) stop('Nind in B and X does not match'))

#--> main parameters
h2 = 0.25
b2 = 0.25
Nqtl_y = 100
Notu_y = 25
Notu_y_g = 25
Nqtl_otu = 10
Nclust = 500
Nmiss = 75

# Clusters st minimum within cluster cor is > 3rd quantile average rho
Cl=hclust(dist(B),method="ward.D2")
Bclust=cutree(Cl,Nclust)

#--> simulate data
s = SimuBiome(X, B=B0, Bclust=Bclust, h2=h2, b2=b2, Nqtl_y=Nqtl_y, Notu_y=Notu_y, Notu_y_g=Notu_y_g)
print(names(s))
# simulated phenotype
y = s$y
# returns reordered B, always in log scale,
# Note: B is shuffled wrt to genotypes in every call to SimuBiome irrespective of the model
B = s$B

#--> plots
par(mfrow=c(2,2))
hist(s$gq, main='Indiv Genetic values')
hist(s$gb, main='Indiv Genetic values')
hist(abs(s$b_qtl), main='Genetic coefficients (abs)')
print(c('Nqtl ',length(s$b_qtl)))
hist(abs(s$b_otu), main='Microbiome coefficients (abs)')
print(c('Notu causative ',length(s$b_otu)))

#--> save simubiome data
save(s,file='simubiome.Rdata')

#--> scale and transpose data
X = scale(t(X))
B = scale(t(B))
y = scale(y)

#--> prediction
tst = sample(seq(length(y)),size=Nmiss)
yNA = y
yNA[tst] = NA

#--> BayesC, WARNING: can take long to run and writes big files
probin=0.001
p0=5
fm_Cgb = doBayesC(yNA, X=X, B=B, out='bayCgb_', pi1=probin, pi2=probin, p0=p0)
fm_Cg = doBayesC(yNA, X=X, out='bayCg_', pi1=probin, p0=p0)
fm_Cb = doBayesC(yNA, B=B, out='bayCb_', pi2=probin, p0=p0)
c = list('fm_Cgb'=fm_Cgb,'fm_Cg'=fm_Cg,'fm_Cb'=fm_Cb)
save(c,file='fmc.Rdata')

# postprocess Bayes C results
load('simubiome.Rdata')
load('fmc.Rdata')
# predictive accuracy in Bayes Cgb
y = s$y
# recover missing if needed
ytst = c$fm_Cgb$fm$y
tst = which(is.na(ytst))

yhat_Cgb = c$fm_Cgb$fm$yHat[tst]
rho = cor(c$fm_Cgb$fm$yHat[tst], y[tst])
# plot y yhat
plot(y[tst], c$fm_Cgb$fm$yHat[tst], xlab='Observed y', ylab='Predicted y')

# computes h2, b2 and cov(g,b)
G=readBinMat('bayCgb_ETA_1_b.bin')
C=readBinMat('bayCgb_ETA_2_b.bin')
u=G%*%t(X)
b=C%*%t(B)
Y=matrix(rep(y,nrow(G)), nrow = nrow(G), byrow=TRUE)
varU=apply(u,1,var)
varB=apply(b,1,var)
covUB=(apply(u+b,1,var) - varU - varB)*0.5
varE=apply(Y-u-b,1,var)
h2g=varU/(varU+varB+varE)
h2b=varB/(varU+varB+varE)
h2gb=covUB/(varU+varB+varE)
# plots
plot(density(h2g),xlim=c(-0.05,1),main=c('Bayes Cgb',rho),xlab='var component')
lines(density(h2b),col='blue',lty=2)
lines(density(h2gb),col='red',lty=3)
print(c(mean(h2g),mean(h2b),mean(h2gb)))


#--> GBLUP: Not used in main paper
if (1<0) {
fm_Ggb = doGBLUP(yNA, X, B, out='gblupgb_')
fm_Gg = doGBLUP(yNA, X=X, out='gblupg_')
fm_Gb = doGBLUP(yNA, B=B, out='gblupb_')
g = list('tst'=tst,'fm_Ggb'=fm_Ggb,'fm_Gg'=fm_Gg,'fm_Gb'=fm_Gb)
save(g,file='fmg.Rdata')
}
#--> This is it
