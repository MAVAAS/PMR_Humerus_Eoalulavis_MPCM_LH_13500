# Script PEM – Lucas Legendre
# SEP 27 – December 9, 2020
# Compiled in R version 4.0.3 (2020-10-10)
# Loading packages
library(MPSEM)
library(evobiR)
library(phytools)
library(Metrics)
library(RColorBrewer)
library(dplyr)
library(textshape)
library(remotes)
library(phylolm)
library(evobiR)
# WARNING: edit the working directory to your preferred folder
# Set the R environment language to English
Sys.setenv(LANGUAGE = "en")

# You can check the current language setting by running:
Sys.getenv("LANGUAGE")

# Phylogeny and dataset
data<-read.table("MMRdataset.txt", header=T)
tree<-read.nexus("untitled.trees.trees (1).trees.nex")
#tree<-drop.tip(tree, setdiff(tree$tip.label, data$species))
plot(tree)
data<-ReorderData(tree, data, taxa.names = 1)
data[,c(2:3)]<-log1p(data[,c(2:3)])
grloc <- getGraphLocations(tree,data[is.na(data[,"MMR"]),"species"])
sporder <- match(attr(grloc$x,"vlabel")[grloc$x$vertex$species],data[,"species"])

# Extract eigenvectors from the tree

PEMfs <- PEM.fitSimple(y=data[sporder,"MMR"],x=NULL,w=grloc$x,d="distance",sp="species",lower=0,upper=1)
# to see all phylog eigenvectors:
#In object tree there are 21 species and 20 nodes
#With data for 14 species so we expect 13 eigenvectors
PEMfs$u
# when we create a table:
b<-PEMfs$u
View(b)
write.table(b, file="eigenvectors.txt", sep="\t")
# the species column is present but the names of eigenvectors are lagged


# Build the PEM
# List of additional variables
aux <- c("Qdot")

PEMfs <- list()
  for(m in aux) {
    PEMfs[[m]] <- PEM.fitSimple(y=data[sporder,"MMR"],
                                        x=data[sporder,m],w=grloc$x,d="distance",sp="species",lower=0,upper=1)
  PEMfs[["none"]] <- PEM.fitSimple(y=data[sporder,"MMR"],
                                           x=NULL,w=grloc$x,d="distance",sp="species",lower=0,upper=1)
} ; rm(m)

for(m in c(aux,"none")) print(PEMfs[[m]]$optim$par) 
# Value of alpha parameter (very small alpha means high phylog signal)
rm(m)


# Select the best model (with one of the co-predictors or without co-predictors) based on AICc
PEMAIC <- list()
for(m in aux) {
  PEMAIC[[m]] <- lmforwardsequentialAICc(y=data[sporder,"MMR"],
                                            x=data[sporder,m,drop=FALSE],object=PEMfs[[m]])
  PEMAIC[["none"]] <- lmforwardsequentialAICc(y=data[sporder,"MMR"],object=PEMfs[["none"]])
}

for(m in c(aux,"none"))
  cat(m,summary(PEMAIC[[m]])$adj,PEMAIC[[m]]$AICc,"\n") # R-squared and AICc value for each co-predictor

rm(m)


summary(PEMAIC[["Qdot"]])

# Predicting missing values for MMR

m <- "Qdot"

atr <- data[is.na(data[,"MMR"]),m,drop=FALSE]
resultsPEM<-predict(object=PEMfs[[m]],targets=grloc,lmobject=PEMAIC[[m]],newdata=atr,interval="confidence")

# Predicted MMR, with upper and lower limits of the confidence interval for the prediction:
exp(resultsPEM$values); exp(resultsPEM$upper); exp(resultsPEM$lower)

#PLot the phylogenetic eigenvectors
#In object tree there are 21 species and 20 nodes
#With data for 14 species so we expect 13 eigenvectors
PEMfs$none$u
c<-PEMfs$none$u
View(c)
write.table(c, file="eigenvectors-bis.txt", sep="\t")
# the species column is present but the names of eigenvectors are lagged
# THEY ARE THE SAME AS ABOVE


### Plot with the phylogeny, values, and confidence intervals ###

# Objects for predicted and fitted values
if(m == "none") {
  ypred <- predict(object=PEMfs[[m]],targets=grloc,lmobject=PEMAIC[[m]],interval="confidence")
} else {
  atr <- data[is.na(data[,"MMR"]),m,drop=FALSE]
  ypred <- predict(object=PEMfs[[m]],targets=grloc,lmobject=PEMAIC[[m]],newdata=atr,interval="confidence")
}
yfit <- numeric(length(PEMAIC[[m]]$fitted.values)) ; yfit[sporder] <- PEMAIC[[m]]$fitted.values

data<-read.table("MMRdataset.txt", header=T) # Reload data to plot it with log10 conversion

# Order of tips in the tree for the plot
tree2<-ladderize(tree)
is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
tree2$tip.label[ordered_tips]

# Plot the tree and axis
layout(matrix(c(1,1,2),1L,3L))
par(mar=c(4,1,1,1))
plot(tree2,cex=1.00)
lab <- log10(range(data[,"MMR"],na.rm=TRUE)) ; lab <- c(floor(lab[1L]),ceiling(lab[2L]))
plot(NA,ylim=c(1L,length(tree2$tip.label)),xlim=lab,axes=FALSE,xlab="MMR (mLO2 h-1 g-0.93)")
lab <- lab[1L]:lab[2L]
axis(1,at=lab,label=10^lab)

# Plot the observed, fitted, and predicted values
points(y=match(data[!is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       x=log10(data[!is.na(data[,"MMR"]),"MMR"]),pch=22,bg="white",cex=2)
points(y=match(data[!is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       x=log10(exp(yfit)),pch=4,cex=2)
points(y=match(data[is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       x=log10(exp(ypred$value)),pch=22,bg="black",cex=2)
arrows(y0=match(data[is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       y1=match(data[is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       x0=log10(exp(ypred$values)),x1=log10(exp(ypred$upper)),length=0.05,angle=90)
arrows(y0=match(data[is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       y1=match(data[is.na(data[,"MMR"]),"species"],tree2$tip.label[ordered_tips]),
       x0=log10(exp(ypred$values)),x1=log10(exp(ypred$lower)),length=0.05,angle=90)


# Leave-one-out cross-validation (LOOCV)
data[,c(2:3)]<-log1p(data[,c(2:3)])

dataCV<-subset(data, !is.na(data$MMR))
dataCV2<-subset(data, !is.na(data$MMR))
m <- "Qdot"

# Estimate PEM predictions for each extant taxon
predictions<-list()
for (i in 1:nrow(dataCV)) {
  dataCV[i,2]<-NA
  treeCV<-drop.tip(tree, setdiff(tree$tip.label, dataCV$species))
  grloc <- getGraphLocations(treeCV,dataCV[is.na(dataCV[,"MMR"]),"species"])
  sporder <- match(attr(grloc$x,"vlabel")[grloc$x$vertex$species],dataCV[,"species"])
  PEMfs<-PEM.fitSimple(y=dataCV[sporder,"MMR"],
                       x=dataCV[sporder,m],w=grloc$x,d="distance",sp="species",lower=0,upper=1)
  PEMAIC<-lmforwardsequentialAICc(y=dataCV[sporder,"MMR"],x=dataCV[sporder,m,drop=FALSE],object=PEMfs)
  atr <- dataCV[is.na(dataCV[,"MMR"]),m,drop=FALSE]
  resultsPEM<-predict(object=PEMfs,targets=grloc,lmobject=PEMAIC,newdata=atr,interval="confidence")
  predictions<-c(predictions,resultsPEM$values)
  dataCV[i,2]<-dataCV2[i,2]
}


# Compile predictions for extant taxa
predata<-vector(length=length(na.omit(data$MMR)))
for (i in 1:length(predata)) {
  predata[i]<-predictions[[i]]
}

str(predictions)
length(predictions)


predata <- vector(length = length(na.omit(data$MMR)))
for (i in seq_along(predata)) {
  if (i <= length(predictions)) {
    predata[i] <- predictions[[i]]
  } else {
    predata[i] <- NA  # Assign NA or another placeholder if predictions are incomplete
  }
}


cat("Length of predictions:", length(predictions), "\n")
cat("Length of predata:", length(predata), "\n")


predata <- sapply(seq_along(na.omit(data$MMR)), function(i) {
  if (i <= length(predictions)) predictions[[i]] else NA
})

names(predata)<-dataCV2$species
dataCV2<-cbind(dataCV2,predata)


# Additional tests
# 1) Wilcoxon signed-rank test: pairwise test for difference between observed and predicted values
wilcox.test(dataCV2$MMR,dataCV2$predata, paired=T)
# Difference not significant: predictions are valid


