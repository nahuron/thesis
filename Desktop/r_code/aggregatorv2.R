#required libraries
library(viridis)
library(FactoMineR)
library(pcaMethods)
library(corrplot)
library(Hmisc)
library(maptools)
library(scales)
library(TeachingDemos)

#mac morphological
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.l/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.m/")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.noTL/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.noTL.l/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.noTL.m/")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.red/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.red.l/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.red.m/")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL/")

#mac genetic
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.l")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.m")

#linux morphological
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph/")
#linux genetic
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")


mypal_viridis <- viridis(20,1,0,1,"D")
  
#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

#morphological
mydata.filepathshort[-grep("morph_fr", mydata.filepathshort)] -> mydata.filepathshort
mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

#genetic
mydata.filepathshort[-grep("phylo_fr", mydata.filepathshort)] -> mydata.filepathshort
mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

bwholder <- NA

#morphological
for (a in 1:length(mydata.filepathshort)){
  
  #set Bandwidth for density functions that equals the average among values
  #bw <- 0.90864586  #all data
  #bw <- 0.4871598 #luzon only
  #bw <- 0.9148217 #mindanao only
  #bw <- 0.5
  #bw <- 0.8817729 #all data noTL
  #bw <- 0.8503455 #all data red
  #bw <- 0.5569555 #luzon red
  #bw <- 0.930679  #mindanao red
  bw <- 1.969154 #d^2 data
  
  #first loop to read in first dataset
  if (a == 1) {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"emp.com.mean.pvalues"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    #sigs <- density(mydata$emp.com.mean.holder[sig.hold], kernel="triangular", bw=bw)
    #nonsigs <- density(mydata$emp.com.mean.holder[nonsig.hold], kernel="triangular", bw=bw)
    
    print(a)
    print(mydata$emp.com.mean.holder[sig.hold])
    print(mydata$emp.com.mean.holder[nonsig.hold])
    
    mydata.density <- density(mydata$emp.com.mean.holder[-as.numeric(nrow(mydata))], kernel="triangular",bw=bw, from=0)
    
    plot(mydata.density, col=mypal_viridis[a], main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Average Mahalanobis D^2 By Community", ylim=c(0,0.15), xlim=c(0,50), xaxs="i", yaxs="i")
    #plot(sigs, col=mypal_viridis[a], ylim=c(0,0.40))
    #lines(nonsigs,lty=2, col=mypal_viridis[a])

    bwholder <- c(bwholder, mydata.density$bw)
    
    sig.hold <- NA
    nonsig.hold <- NA
    #print(range(c(sigs$y, nonsigs$y)))
  }
  else  {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"emp.com.mean.pvalues"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    #sigs <- density(mydata$emp.com.mean.holder[sig.hold], kernel="triangular", bw=bw)
    #nonsigs <- density(mydata$emp.com.mean.holder[nonsig.hold], kernel="triangular", bw=bw)
    
    print(a)
    print(mydata$emp.com.mean.holder[sig.hold])
    print(mydata$emp.com.mean.holder[nonsig.hold])
    
    mydata.density <- density(mydata$emp.com.mean.holder[-as.numeric(nrow(mydata))], kernel="triangular",bw=bw, from=0)
    
    lines(mydata.density, col=mypal_viridis[a])
    #lines(sigs, col=mypal_viridis[a])
    #lines(nonsigs,lty=2, col=mypal_viridis[a])
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    sig.hold <- NA
    nonsig.hold <- NA
    #print(range(c(sigs$y, nonsigs$y)))
  }
  lines(density(answer@sim.com.mdistances, kernel="triangular", bw=bw, from=0), lty=2, lwd=1.5)
}

legend(x=40, y=0.145, legend=format(round(seq(0.05, 1.00, 0.05),2), nsmall=2), col=mypal_viridis[1:length(mydata.filepathshort)], lty=1, lwd=2, xpd=TRUE)

bwholder <- NA

#PSV
for (a in 1:length(mydata.filepathshort)){
  #for (a in 1){  
  #bw <- 0.07204041  #determined from sig data
  bw <- 0.0571599773 #determined from all data
  #bw <- 0.04171355   #determined from luzon data
  #bw <- 0.05521344 #determined from all mindanao data
  
  if (a == 1) {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"psv.p"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    mydata.density <- density(mydata$psv[-as.numeric(nrow(mydata))], from=0, to=1, bw=bw)
    
    sigs <- density(mydata$psv[sig.hold], kernel="triangular", from=0, to=1, bw=bw)
    nonsigs <- density(mydata$psv[nonsig.hold], kernel="triangular", from=0, to=1, bw=bw)
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    print(a)
    print(mydata$psv[sig.hold])
    print(mydata$psv[nonsig.hold])
    print(mydata$n.coms[sig.hold])
    print(mydata$n.coms[nonsig.hold])
    #print(mydata.density$bw)
    #print(range(c(sigs$y, nonsigs$y)))
    
    plot(mydata.density, col=mypal_viridis[a],main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability", ylim=c(0,7), xaxs="i", yaxs="i")
    #plot(sigs, col=mypal_viridis[a], ylim=ceiling(range(c(sigs$y, nonsigs$y))), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability")
    #lines(nonsigs,lty=2, col=mypal_viridis[a])
    
    sig.hold <- NA
    nonsig.hold <- NA
  }
  else  {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"psv.p"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    mydata.density <- density(mydata$psv[-as.numeric(nrow(mydata))], from=0, to=1, bw=bw)
    
    sigs <- density(mydata$psv[sig.hold], kernel="triangular", from=0, to=1, bw=bw)
    nonsigs <- density(mydata$psv[nonsig.hold], kernel="triangular", from=0, to=1, bw=bw)
    
    bwholder <- c(bwholder, mydata.density$bw)
    print(a)
    print(mydata$psv[sig.hold])
    print(mydata$psv[nonsig.hold])
    print(mydata$n.coms[sig.hold])
    print(mydata$n.coms[nonsig.hold])
    #print(mydata.density$bw)
    #print(range(c(sigs$y, nonsigs$y)))
    
    lines(mydata.density, col=mypal_viridis[a],main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability")
    #lines(sigs, col=mypal_viridis[a], ylim=range(c(sigs$y, nonsigs$y)))
    #lines(nonsigs, col=mypal_viridis[a], lty=2)
    
    sig.hold <- NA
    nonsig.hold <- NA
  }
}
#all PAIC
abline(v=0.854, lty=2)
#Luzon PAIC
abline(v=0.851, lty=2)
#Mindanao PAIC
abline(v=0.739, lty=2)

legend(x=0.05, y=6.5, legend=format(round(seq(0.05, 1.00, 0.05),2), nsmall=2), col=mypal_viridis[1:length(mydata.filepathshort)], lty=1, lwd=2, xpd=TRUE)

#SES.MPD
for (a in 1:length(mydata.filepathshort)){
#for (a in 1){  
  
  bw <-  0.455850724  #determined from all data
  #bw <- 0.3110999 #determined from Luzon data
  #bw <- 0.3230618 #determined from Mindanao data
  
  if (a == 1) {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"ses.mpd.p"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    print(a)
    print(mydata$ses.mpd[sig.hold])
    print(mydata$ses.mpd[nonsig.hold])
    print(mydata$n.coms[sig.hold])
    print(mydata$n.coms[nonsig.hold])
    
    if(length(sig.hold[!is.na(sig.hold)])>=1 & length(nonsig.hold[!is.na(nonsig.hold)])>=1)  {
    
      mydata.density <- density(mydata$ses.mpd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      
    sigs <- density(mydata$ses.mpd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    nonsigs <- density(mydata$ses.mpd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    print(range(c(sigs$y, nonsigs$y)))
    print(range(c(sigs$x, nonsigs$x)))
    
    plot(mydata.density, col=mypal_viridis[a], ylim=c(0,1), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance", xaxs="i", yaxs="i")
    #plot(sigs, col=mypal_viridis[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
    #lines(nonsigs,lty=2, col=mypal_viridis[a])
    
    sig.hold <- NA
    nonsig.hold <- NA
    }
    else{
      mydata.density <- density(mydata$ses.mpd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)
      
      bwholder <- c(bwholder, mydata.density$bw)
      
      print(range(mydata.density$y))
      print(range(mydata.density$x))
      
      plot(mydata.density, col=mypal_viridis[a], ylim=c(0,0.90), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      
    }
  }
  else  {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"ses.mpd.p"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    print(a)
    print(mydata$ses.mpd[sig.hold])
    print(mydata$ses.mpd[nonsig.hold])
    print(mydata$n.coms[sig.hold])
    print(mydata$n.coms[nonsig.hold])
    
    if(length(sig.hold[!is.na(sig.hold)])>=1 & length(nonsig.hold[!is.na(nonsig.hold)])>=1)  {
      
    mydata.density <- density(mydata$ses.mpd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
    
    sigs <- density(mydata$ses.mpd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    nonsigs <- density(mydata$ses.mpd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    
    lines(mydata.density, col=mypal_viridis[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
    #lines(sigs, col=mypal_viridis[a])
    #lines(nonsigs, col=mypal_viridis[a], lty=2)
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    print(range(c(sigs$y, nonsigs$y)))
    print(range(c(sigs$x, nonsigs$x)))
    
    sig.hold <- NA
    nonsig.hold <- NA
    }
    else{
      mydata.density <- density(mydata$ses.mpd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)
      
      bwholder <- c(bwholder, mydata.density$bw)
      
      print(range(mydata.density$y))
      print(range(mydata.density$x))
      
      lines(mydata.density, col=mypal_viridis[a], ylim=c(0,0.75), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      
    }
  }
}
abline(v=0, lty=2)
legend(x=3.5, y=0.90, legend=format(round(seq(0.05, 1.00, 0.05),2), nsmall=2), col=mypal_viridis[1:length(mydata.filepathshort)], lty=1, lwd=2, xpd=TRUE)

#SES.MNTD
for (a in 1:length(mydata.filepathshort)){
  #for (a in 1){  
  
  bw <- 0.441174087 #determined from all data
  #bw <- 0.3110999 #determined from Luzon data
  #bw <- 0.3230618 #determined from Mindanao data
  
  if (a == 1) {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"ses.mntd.p"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    print(a)
    print(mydata$ses.mntd[sig.hold])
    print(mydata$ses.mntd[nonsig.hold])
    
    if(length(sig.hold[!is.na(sig.hold)])>=1 & length(nonsig.hold[!is.na(nonsig.hold)])>=1)  {
      
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      
      sigs <- density(mydata$ses.mntd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      nonsigs <- density(mydata$ses.mntd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(c(sigs$y, nonsigs$y)))
      print(range(c(sigs$x, nonsigs$x)))
      
      plot(mydata.density, col=mypal_viridis[a], ylim=c(0,1), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Nearest Taxon Distance", xaxs="i", yaxs="i")
      #plot(sigs, col=mypal_viridis[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Nearest Taxon Distance")
      #lines(nonsigs,lty=2, col=mypal_viridis[a])
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
    else{
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)
      
      bwholder <- c(bwholder, mydata.density$bw)
      
      print(range(mydata.density$y))
      print(range(mydata.density$x))
      
      plot(mydata.density, col=mypal_viridis[a], ylim=c(0,1), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Nearest Taxon Distance")
      
    }
  }
  else  {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    pvaluesorted <- unlist(lapply(mydata[-as.numeric(nrow(mydata)),"ses.mntd.p"], function (x) x <= 0.025 | x >= 0.975))
    
    sig.hold <- NA
    nonsig.hold <- NA
    
    for(b in 1:length(pvaluesorted)){
      if (pvaluesorted[b]==TRUE){
        sig.hold <- c(sig.hold, b)
        sig.hold <- sig.hold[!is.na(sig.hold)]
      } 
      else  {
        nonsig.hold <- c(nonsig.hold, b)
        nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
      }
    }
    
    print(a)
    print(mydata$ses.mntd[sig.hold])
    print(mydata$ses.mntd[nonsig.hold])
    
    if(length(sig.hold[!is.na(sig.hold)])>=1 & length(nonsig.hold[!is.na(nonsig.hold)])>=1)  {
      
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      
      sigs <- density(mydata$ses.mntd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      nonsigs <- density(mydata$ses.mntd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      
      lines(mydata.density, col=mypal_viridis[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      #lines(sigs, col=mypal_viridis[a])
      #lines(nonsigs, col=mypal_viridis[a], lty=2)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(c(sigs$y, nonsigs$y)))
      print(range(c(sigs$x, nonsigs$x)))
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
    else if(length(sig.hold)<=1){
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      nonsigs <- density(mydata$ses.mntd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      lines(mydata.density, col=mypal_viridis[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Nearest Taxon Distance")
      #lines(sigs, col=mypal_viridis[a])
      #lines(nonsigs, col=mypal_viridis[a], lty=2)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(nonsigs$y))
      print(range(nonsigs$x))
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
    else {
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      
      #lines(sigs, col=mypal_viridis[a])
      #lines(nonsigs, col=mypal_viridis[a], lty=2)
      
      bwholder <- c(bwholder, mydata.density$bw)
      
      print(range(mydata.density$y))
      print(range(mydata.density$x))
      
      lines(mydata.density, col=mypal_viridis[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
  }
}
abline(v=0, lty=2)
legend(x=3.5, y=0.90, legend=format(round(seq(0.05, 1.00, 0.05),2), nsmall=2), col=mypal_viridis[1:length(mydata.filepathshort)], lty=1, lwd=2, xpd=TRUE)

mydata.means <- rep(NA, times=length(35:(ncol(mydata)-1)))

#obtain average values of communities for each scale and store in an object
#genetic
for (a in 1:length(mydata.filepathshort)){
  
    if (a==1) {
      mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
      mydata.means <- rep(NA, times=length(35:(ncol(mydata)-1)))
      
      mydata.means.hold <- mydata[nrow(mydata),35:(ncol(mydata)-1)]
      
      mydata.means <- rbind(mydata.means, mydata.means.hold)
      mydata.means <- mydata.means[-1,]
      
    }
  else{
  mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
  mydata.means.hold <- mydata[nrow(mydata),35:(ncol(mydata)-1)]
  
 
  print(mydata.means)
  mydata.means <- rbind(mydata.means, mydata.means.hold)
  
  }  
  rm(mydata)
}

#morphology
for (a in 1:length(mydata.filepathshort)){
  
  if (a==1) {
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    mydata.means <- rep(NA, times=length(35:(ncol(mydata))))
    
    mydata.means.hold <- mydata[nrow(mydata),35:(ncol(mydata))]
    
    mydata.means <- rbind(mydata.means, mydata.means.hold)
    mydata.means <- mydata.means[-1,]
    
  }
  else{
    mydata <- read.csv(mydata.filepath[a], header=T, row.names=1)
    mydata.means.hold <- mydata[nrow(mydata),35:(ncol(mydata))]
    
    
    print(mydata.means)
    mydata.means <- rbind(mydata.means, mydata.means.hold)
    
  }  
  rm(mydata)
}


#ses.mpd
pvaluesorted <- unlist(lapply(mydata.means[,"ses.mpd.p"], function (x) x <= 0.025 | x >= 0.975))

#ses.mntd
pvaluesorted <- unlist(lapply(mydata.means[,"ses.mntd.p"], function (x) x <= 0.025 | x >= 0.975))

#psv
pvaluesorted <- unlist(lapply(mydata.means[,"psv.p"], function (x) x <= 0.025 | x >= 0.975))

#morphology
pvaluesorted <- unlist(lapply(mydata.means[,"emp.com.mean.pvalues"], function (x) x <= 0.025 | x >= 0.975))


sig.hold <- NA
nonsig.hold <- NA

for(b in 1:length(pvaluesorted)){
  if (pvaluesorted[b]==TRUE){
    sig.hold <- c(sig.hold, b)
    sig.hold <- sig.hold[!is.na(sig.hold)]
  } 
  else  {
    nonsig.hold <- c(nonsig.hold, b)
    nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
    
  }
}
if(length(sig.hold)==1 && is.na(sig.hold)==TRUE){
  sig.hold <- NULL
}

bw <- 0.07395
print(sig.hold);print(nonsig.hold)

#ses.mpd
#bw <- 0.03717#luzon only
#bw <- 0.07039#mindanao only
plot(density(mydata.means$ses.mpd, bw=bw), xlim=c(-1,1), ylim=c(0,5), main="Density Curves for Mean of Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
abline(v=0, lty=2)
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$ses.mpd[sig.hold], kernel="triangular", bw=bw), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$ses.mpd[nonsig.hold], kernel="triangular", bw=bw), col="red", lty=4)
}
legend(x=0.65, y=4.5, legend=c("All", "Significant", "N.S.", "Random \nDispersion"), lty=c(1,1,4,2), col=c("black","blue","red", "black"), cex=0.75)

#ses.mntd
#bw <- 0.05424#luzon only
#bw <- 0.08518#mindanao only
plot(density(mydata.means$ses.mntd), ylim=c(0,6), xlim=c(-1,1), main="Density Curves for Mean of Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Nearest Taxon Distance")
abline(v=0, lty=2)
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$ses.mntd[sig.hold], kernel="triangular", bw=bw), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$ses.mntd[nonsig.hold], kernel="triangular", bw=bw), col="red", lty=4)
}
legend(x=0.65, y=5.5, legend=c("All", "Significant", "N.S.", "Random \nDispersion"), lty=c(1,1,4,2), col=c("black","blue","red", "black"), cex=0.75)

#psv
#bw <- 0.004874#luzon only
#bw <- 0.0127#mindanao only
plot(density(mydata.means$psv, bw=bw), xlim=c(0.5,1), ylim=c(0,15), main="Density Curves for Mean of Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability")
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$psv[sig.hold], kernel="triangular", bw=bw), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$psv[nonsig.hold], kernel="triangular", bw=bw), col="red", lty=4)
}
abline(v=0.854, lty=2)
legend(x=0.55, y=12.5, legend=c("All", "Significant", "N.S.", "Random \nDispersion"), lty=c(1,1,4,2), col=c("black","blue","red", "black"), cex=0.75)


#morphology
plot(density(mydata.means$emp.com.mean.holder),ylim=c(0,4), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales (D^2 noTL)", xlab="Average Morphological Disparity By Community", col="black")
if(length(sig.hold[!is.na(sig.hold)]) > 0)  {
  lines(density(mydata.means$emp.com.mean.holder[sig.hold], kernel="triangular"), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$emp.com.mean.holder[nonsig.hold], kernel="triangular"), col="red", lty=4)
}
legend(x=3.1, y=3.5, legend=c("All", "Significant", "N.S.", "Random \nDispersion"), lty=c(1,1,4,2), col=c("black","blue","red", "black"), cex=0.75)
