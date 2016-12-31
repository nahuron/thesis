#mac morphological
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph/")
#mac genetic
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

#linux morphological
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph/")
#linux genetic
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")


tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "fr\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = "fr.csv$", full.names = FALSE))

bwholder <- NA

#morphological
for (a in 1:length(mydata.filepathshort)){
  
  #set Bandwidth for density functions that equals the average among values
  bw <- 0.90864586
  
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
    
    sigs <- density(mydata$emp.com.mean.holder[sig.hold], kernel="triangular", bw=bw)
    nonsigs <- density(mydata$emp.com.mean.holder[nonsig.hold], kernel="triangular", bw=bw)
    
    print(mydata$emp.com.mean.holder[sig.hold])
    print(mydata$emp.com.mean.holder[nonsig.hold])
    
    mydata.density <- density(mydata$emp.com.mean.holder[-as.numeric(nrow(mydata))], kernel="triangular", bw=bw)
    
    plot(mydata.density, col=tol21rainbow[a], main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Average Morphological Disparity By Community", ylim=c(0,0.25))
    #plot(sigs, col=tol21rainbow[a], ylim=c(0,0.40))
    #lines(nonsigs,lty=5, col=tol21rainbow[a])

    bwholder <- c(bwholder, mydata.density$bw)
    
    sig.hold <- NA
    nonsig.hold <- NA
    print(range(c(sigs$y, nonsigs$y)))
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
    
    sigs <- density(mydata$emp.com.mean.holder[sig.hold], kernel="triangular", bw=bw)
    nonsigs <- density(mydata$emp.com.mean.holder[nonsig.hold], kernel="triangular", bw=bw)
    
    print(mydata$emp.com.mean.holder[sig.hold])
    print(mydata$emp.com.mean.holder[nonsig.hold])
    
    mydata.density <- density(mydata$emp.com.mean.holder[-as.numeric(nrow(mydata))], kernel="triangular", bw=bw)
    
    lines(mydata.density, col=tol21rainbow[a])
    #lines(sigs, col=tol21rainbow[a])
    #lines(nonsigs,lty=5, col=tol21rainbow[a])
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    sig.hold <- NA
    nonsig.hold <- NA
    print(range(c(sigs$y, nonsigs$y)))
  }
  lines(density(answer@sim.com.mdistances, kernel="triangular", bw=bw), lty=5, lwd=2)
}

legend(x=12, y=0.225, legend=seq(0.10, 1.00, 0.05), col=tol21rainbow[1:length(mydata.filepathshort)], lty=1, lwd=2, xpd=TRUE)

bwholder <- NA

#PSV
for (a in 1:length(mydata.filepathshort)){
  #for (a in 1){  
  #bw <- 0.07204041  #determined from sig data
  bw <- 0.0571599773 #determined from all data
  
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
    print(mydata$n.coms[nonsig.hold])
    print(sigs$bw)
    print(range(c(sigs$y, nonsigs$y)))
    
    plot(mydata.density, col=tol21rainbow[a],main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability")
    #plot(sigs, col=tol21rainbow[a], ylim=ceiling(range(c(sigs$y, nonsigs$y))), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability")
    #lines(nonsigs,lty=5, col=tol21rainbow[a])
    
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
    print(mydata$n.coms[nonsig.hold])
    print(sigs$bw)
    print(range(c(sigs$y, nonsigs$y)))
    
    lines(mydata.density, col=tol21rainbow[a],main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Phylogenetic Species Variability")
    #lines(sigs, col=tol21rainbow[a], ylim=range(c(sigs$y, nonsigs$y)))
    #lines(nonsigs, col=tol21rainbow[a], lty=5)
    
    sig.hold <- NA
    nonsig.hold <- NA
  }
}

#SES.MPD
for (a in 1:length(mydata.filepathshort)){
#for (a in 1){  
  
  bw <-  0.455850724  #determined from all data
  
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
    
    if(length(sig.hold)>=1 & length(nonsig.hold)>=1)  {
    
      mydata.density <- density(mydata$ses.mpd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      
    sigs <- density(mydata$ses.mpd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    nonsigs <- density(mydata$ses.mpd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    print(range(c(sigs$y, nonsigs$y)))
    print(range(c(sigs$x, nonsigs$x)))
    
    plot(mydata.density, col=tol21rainbow[a], ylim=c(0,0.75), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
    #plot(sigs, col=tol21rainbow[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
    #lines(nonsigs,lty=5, col=tol21rainbow[a])
    
    sig.hold <- NA
    nonsig.hold <- NA
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
    
    if(length(sig.hold)>=1 & length(nonsig.hold)>=1)  {
      
    mydata.density <- density(mydata$ses.mpd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
    
    sigs <- density(mydata$ses.mpd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    nonsigs <- density(mydata$ses.mpd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
    
    lines(mydata.density, col=tol21rainbow[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
    #lines(sigs, col=tol21rainbow[a])
    #lines(nonsigs, col=tol21rainbow[a], lty=5)
    
    bwholder <- c(bwholder, mydata.density$bw)
    
    print(range(c(sigs$y, nonsigs$y)))
    print(range(c(sigs$x, nonsigs$x)))
    
    sig.hold <- NA
    nonsig.hold <- NA
    }
  }
}

#SES.MNTD
for (a in 1:length(mydata.filepathshort)){
  #for (a in 1){  
  
  bw <- 0.441174087 #determined from all data
  
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
    
    if(length(sig.hold)>=1 & length(nonsig.hold)>=1)  {
      
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      
      sigs <- density(mydata$ses.mntd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      nonsigs <- density(mydata$ses.mntd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(c(sigs$y, nonsigs$y)))
      print(range(c(sigs$x, nonsigs$x)))
      
      plot(mydata.density, col=tol21rainbow[a], ylim=c(0,0.75), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      #plot(sigs, col=tol21rainbow[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Nearest Taxon Distance")
      #lines(nonsigs,lty=5, col=tol21rainbow[a])
      
      sig.hold <- NA
      nonsig.hold <- NA
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
      
      lines(mydata.density, col=tol21rainbow[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      #lines(sigs, col=tol21rainbow[a])
      #lines(nonsigs, col=tol21rainbow[a], lty=5)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(c(sigs$y, nonsigs$y)))
      print(range(c(sigs$x, nonsigs$x)))
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
    else if(length(sig.hold)<=1){
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      nonsigs <- density(mydata$ses.mntd[nonsig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      lines(mydata.density, col=tol21rainbow[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      #lines(sigs, col=tol21rainbow[a])
      #lines(nonsigs, col=tol21rainbow[a], lty=5)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(nonsigs$y))
      print(range(nonsigs$x))
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
    else if(length(nonsig.hold)<=1){
      mydata.density <- density(mydata$ses.mntd[-as.numeric(nrow(mydata))], kernel="triangular", from=-6, to=6, bw=bw)  
      sigs <- density(mydata$ses.mntd[sig.hold], kernel="triangular", from=-6, to=6, bw=bw)
      lines(mydata.density, col=tol21rainbow[a], ylim=c(0,6), xlim=c(-6,6), main="Density Curves for Empirical Communities \nAcross a Range of Spatial Scales", xlab="Standardized Effect Size Mean Phylogenetic Distance")
      #lines(sigs, col=tol21rainbow[a])
      #lines(nonsigs, col=tol21rainbow[a], lty=5)
      
      bwholder <- c(bwholder, mydata.density$bw)
      print(range(sigs$y))
      print(range(sigs$x))
      
      sig.hold <- NA
      nonsig.hold <- NA
    }
  }
}

mydata.means <- rep(NA, times=length(35:(ncol(mydata)-1)))

#obtain average values of communities for each scale and store in an object
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
  } else  {
    nonsig.hold <- c(nonsig.hold, b)
    nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
  }
}
#bw <- 0.04395
print(sig.hold);print(nonsig.hold)

#ses.mpd
plot(density(mydata.means$ses.mpd), xlim=c(-1,1))
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$ses.mpd[sig.hold], kernel="triangular"), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$ses.mpd[nonsig.hold], kernel="triangular"), col="blue", lty=5)
}

#ses.mntd
plot(density(mydata.means$ses.mntd), ylim=c(0,6), xlim=c(-1,1))
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$ses.mntd[sig.hold], kernel="triangular"), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$ses.mntd[nonsig.hold], kernel="triangular"), col="blue", lty=5)
}

#psv
plot(density(mydata.means$psv), xlim=c(0.75,1))
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$psv[sig.hold], kernel="triangular"), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$psv[nonsig.hold], kernel="triangular"), col="blue", lty=5)
}

#morphology
plot(density(mydata.means$emp.com.mean.holder))
if(length(sig.hold) > 0)  {
  lines(density(mydata.means$emp.com.mean.holder[sig.hold], kernel="triangular"), col="blue")
}
if(length(nonsig.hold) > 0)  {
  lines(density(mydata.means$emp.com.mean.holder[nonsig.hold], kernel="triangular"), col="blue", lty=5)
}







#script to split the results by significance into two new objects
pvaluesorted <- lapply(mydata[,c("ses.mpd.p","ses.mntd.p","psv.p")], function (x) x <= 0.025 | x >= 0.975)

sig.hold <- NA
nonsig.hold <- NA

for(b in 1:length(pvaluesorted)){
if (pvaluesorted[b]==TRUE){
  sig.hold <- c(sig.hold, b)
  sig.hold <- sig.hold[!is.na(sig.hold)]
} else  {
    nonsig.hold <- c(nonsig.hold, b)
    nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
  }
}
