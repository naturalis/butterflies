library(ape)
library(nlme)
library(caper)
library(optparse)

setwd("I:\\Vlinders\\Phylogeny")

# Read tree
butterflies <- read.nexus("grafted.pruned.nex")
names(butterflies)
plot(butterflies)
axisPhylo()

data<-read.csv("DataGLS.csv", header=TRUE, as.is=TRUE)
attach(data)


# Analysis of Range size

# exclude rows with NA
sel <- is.na(Range.size)
drop <- Range.size[!sel]
data_rangesize <- data[!sel, ]

# Plot
rownames(data_rangesize) <- data_rangesize$Species_in_Phylo
PC.C12 <- data_rangesize$PC.C1*data_rangesize$PC.C1
PC.C22 <- data_rangesize$PC.C2*data_rangesize$PC.C2
plot(data_rangesize$PC.B1,data_rangesize$Range.size)

# Trim the tree and keep the tips for the species with trait data
sel <- butterflies$tip.label %in% data_rangesize$Species_in_Phylo 
drop <- butterflies$tip.label[!sel]
butterflies.trim <- drop.tip(butterflies,drop)
plot(butterflies.trim)

# Weight for non-ultrametric tree
# see http://blog.phytools.org/2012/04/using-nlmegls-for-phylogenetic.html
weight<-diag(vcv.phylo(butterflies.trim))


# GLS Range size

# GLS using ape

cor.matrix <- corBrownian(phy=butterflies.trim)
modelA <- gls(Range_size ~ PC.B1+PC.B2+PC.B3+PC.C1+PC.C12+PC.C2+PC.C22,data=data_rangesize,correlation=cor.matrix,weights=varFixed(~weight),method="REML")
summary(modelA)
hist(residuals(modelA))
R2A <- cor(data_rangesize$Range_size,predict(modelA))^2
R2A

modelA <- gls(Range_size ~ PC.B1+PC.B2+PC.B3+PC.C1+PC.C12+PC.C2+PC.C22,data=data_rangesize,correlation=cor.matrix,weights=varFixed(~weight))
summary(modelA)
hist(residuals(modelA))
R2A <- cor(data_rangesize$Range_size,predict(modelA))^2
R2A

modelB <- gls(Range.size ~ PC.B1+PC.B2+PC.B3+PC.C1+PC.C12+PC.C2+PC.C22,data=data_rangesize)
summary(modelB)
hist(residuals(modelB))
R2B <- cor(data_rangesize$Range.size,predict(modelB))^2
R2B



# PGLS using caper

# merge into a comparative data object
cdat <- comparative.data(
  phy=butterflies.trim,
  data=data_rangesize,
  names.col="tip.label", # column to join on
  vcv=TRUE, # create variance/covariance matrix, needed for pgls
  na.omit=FALSE, # do NOT omit species with missing data
  warn.dropped=TRUE # warn about any species not shared by tree & data
);

model.pgls.range_size<-pgls(Range_size~PC.B1+PC.B2+PC.B3+PC.C1+PC.C12+PC.C2+PC.C22,data=cdat,lambda="ML");
summary(model.pgls.range_size);





# Analysis of Habitat specificity

# exclude rows with NA
sel <- is.na(SSI)
drop <- SSI[!sel]
data_SSI <- data[!sel, ]

# Plot
rownames(data_SSI) <- data_SSI$Species_in_Phylo
PC.C12 <- data_SSI$PC.C1*data_SSI$PC.C1
PC.C22 <- data_SSI$PC.C2*data_SSI$PC.C2
plot(data_SSI$PC.B1,data_SSI$SSI)

# Trim the tree and keep the tips for the species with trait data
sel <- butterflies$tip.label %in% data_SSI$Species_in_Phylo 
drop <- butterflies$tip.label[!sel]
butterflies.trim <- drop.tip(butterflies,drop)
plot(butterflies.trim)

# GLS SSI
cor.matrix <- corBrownian(phy=butterflies.trim)
modelA <- gls(SSI ~ PC.B1+PC.B2+PC.B3+PC.C1+PC.C12+PC.C2+PC.C22,data=data_SSI,correlation=cor.matrix)
summary(modelA)
hist(residuals(modelA))
R2A <- cor(data_SSI$SSI,predict(modelA))^2
R2A

modelB <- gls(SSI ~ PC.B1+PC.B2+PC.B3+PC.C1+PC.C12+PC.C2+PC.C22,data=data_SSI)
summary(modelB)
hist(residuals(modelB))
R2B <- cor(data_SSI$SSI,predict(modelB))^2
R2B

