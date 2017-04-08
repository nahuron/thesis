#community summary stats GENETIC

####GENETIC DATA####
#set working directory
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.l")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.m")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.standard/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.standard.l")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.standard.m")

#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

#genetic
mydata.filepathshort[-grep("phylo_fr", mydata.filepathshort)] -> mydata.filepathshort
mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

{results.psv <- list()
  length(results.psv) <- length(mydata.filepathshort)
  results.mpd <- list()
  length(results.mpd) <- length(mydata.filepathshort)
  results.mntd <- list()
  length(results.mntd) <- length(mydata.filepathshort)}

#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)

#PSV
for (e in 1:length(mydata.filepathshort)){
  mydata.hold <- read.csv(mydata.filepath[e], header=T, row.names=1)
  names(results.psv)[e] <- paste0("PSV_",mydata.filepathshort[e])
  
  if(colnames(mydata.hold)[as.numeric(ncol(mydata.hold))]=="com.fr"){
    mydata.red <- mydata.hold[,!(names(mydata.hold) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata.hold[-as.numeric(nrow(mydata.hold)),]
  }
  
  mydata.sig <- mydata.red[mydata.red$psv.p>=0.975 | mydata.red$psv.p<=0.025,]
  mydata.rm <- c("ses.mpd", "ses.mpd.p", "ses.mntd", "ses.mntd.p")
  mydata.sig <- mydata.sig[,!(names(mydata.sig) %in% mydata.rm)]
  mydata.sig <- mydata.sig[,colSums(mydata.sig) !=0]
  mydata.sig <- unique(mydata.sig)
  
  if(nrow(mydata.sig)>0){
    com.fr <- rep(NA, times=nrow(mydata.sig))
    for (a in 1:(nrow(mydata.sig))) { #isolate row of interest
      com.match.hold <- 0
      for (b in 1:(ncol(mydata.sig)-3)) { #isolate cell within row of interest
        if (mydata.sig[a,b]==1){
          m <- match(as.character(colnames(mydata.sig[a,][b])), as.character(brach_fr_key[,1]))
          com.match.hold <- c(com.match.hold, as.character(brach_fr_key$FR.code[m]))
          com.match.hold <- com.match.hold[com.match.hold!=0]
          com.match.hold <- com.match.hold[!is.na(com.match.hold)]
          com.match.hold <- unique(com.match.hold)
          print(com.match.hold)
        }
      }
      com.match.hold <- unique(com.match.hold)
      #print(com.match.hold)
      if (length(com.match.hold) > 1 ) {
        #print(paste0("Community ", rownames(mydata.sig)[a], " spans more than 1 Faunal Regions!"))
        #print(com.match.hold)
        com.fr[a] <- as.character(length(com.match.hold))
      }
      else if (length(com.match.hold == 1)) {
        com.fr[a] <- as.character(brach_fr_key$FR.code[m])
      }
      
      mydata.revised <- cbind(mydata.sig, com.fr)
      
      rm(com.match.hold)
    }
    
    results.psv [[e]] <- mydata.revised
  }
}

#MPD
for (e in 1:length(mydata.filepathshort)){
  
  mydata.hold <- read.csv(mydata.filepath[e], header=T, row.names=1)
  names(results.mpd)[e] <- paste0("MPD_",mydata.filepathshort[e])
  
  if(colnames(mydata.hold)[as.numeric(ncol(mydata.hold))]=="com.fr"){
    mydata.red <- mydata.hold[,!(names(mydata.hold) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata.hold[-as.numeric(nrow(mydata.hold)),]
  }
  
  mydata.sig <- mydata.red[mydata.red$ses.mpd.p>=0.975 | mydata.red$ses.mpd.p<=0.025,]
  mydata.rm <- c("psv", "psv.p", "ses.mntd", "ses.mntd.p")
  mydata.sig <- mydata.sig[,!(names(mydata.sig) %in% mydata.rm)]
  mydata.sig <- mydata.sig[,colSums(mydata.sig) !=0]
  mydata.sig <- unique(mydata.sig)
  
  if(nrow(mydata.sig)>0){
    #reassigns FRs
    com.fr <- rep(NA, times=nrow(mydata.sig))
    for (a in 1:(nrow(mydata.sig))) { #isolate row of interest
      com.match.hold <- 0
      for (b in 1:(ncol(mydata.sig)-3)) { #isolate cell within row of interest
        if (mydata.sig[a,b]==1){
          m <- match(as.character(colnames(mydata.sig[a,][b])), as.character(brach_fr_key[,1]))
          com.match.hold <- c(com.match.hold, as.character(brach_fr_key$FR.code[m]))
          com.match.hold <- com.match.hold[com.match.hold!=0]
          com.match.hold <- com.match.hold[!is.na(com.match.hold)]
          print(com.match.hold)
        }
      }
      com.match.hold <- unique(com.match.hold)
      print(com.match.hold)
      if (length(com.match.hold) > 1 ) {
        print(paste0("Community ", rownames(mydata.sig)[a], " spans more than 1 Faunal Regions!"))
        print(com.match.hold)
        com.fr[a] <- as.character(length(com.match.hold))
      }
      else if (length(com.match.hold == 1)) {com.fr[a] <- as.character(brach_fr_key$FR.code[m])}
      
      mydata.revised <- cbind(mydata.sig, com.fr)
      
      rm(com.match.hold)
    }
    
    #store results
    results.mpd [[e]] <- mydata.revised
  } 
}

#MNTD
for (e in 1:length(mydata.filepathshort)){
  
  mydata.hold <- read.csv(mydata.filepath[e], header=T, row.names=1)
  names(results.mntd)[e] <- paste0("MNTD_",mydata.filepathshort[e])
  
  if(colnames(mydata.hold)[as.numeric(ncol(mydata.hold))]=="com.fr"){
    mydata.red <- mydata.hold[,!(names(mydata.hold) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata.hold[-as.numeric(nrow(mydata.hold)),]
  }
  
  mydata.sig <- mydata.red[mydata.red$ses.mntd.p>=0.975 | mydata.red$ses.mntd.p<=0.025,]
  mydata.rm <- c("psv", "psv.p", "ses.mpd", "ses.mpd.p")
  mydata.sig <- mydata.sig[,!(names(mydata.sig) %in% mydata.rm)]
  mydata.sig <- mydata.sig[,colSums(mydata.sig) !=0]
  mydata.sig <- unique(mydata.sig)
  
  if(nrow(mydata.sig)>0){
    
    #reassigns FRs
    com.fr <- rep(NA, times=nrow(mydata.sig))
    for (a in 1:(nrow(mydata.sig))) { #isolate row of interest
      com.match.hold <- 0
      for (b in 1:(ncol(mydata.sig)-3)) { #isolate cell within row of interest
        if (mydata.sig[a,b]==1){
          m <- match(as.character(colnames(mydata.sig[a,][b])), as.character(brach_fr_key[,1]))
          com.match.hold <- c(com.match.hold, as.character(brach_fr_key$FR.code[m]))
          com.match.hold <- com.match.hold[com.match.hold!=0]
          com.match.hold <- com.match.hold[!is.na(com.match.hold)]
          print(com.match.hold)
        }
      }
      com.match.hold <- unique(com.match.hold)
      print(com.match.hold)
      if (length(com.match.hold) > 1 ) {
        print(paste0("Community ", rownames(mydata.sig)[a], " spans more than 1 Faunal Regions!"))
        print(com.match.hold)
        com.fr[a] <- as.character(length(com.match.hold))
      }
      else if (length(com.match.hold == 1)) {com.fr[a] <- as.character(brach_fr_key$FR.code[m])}
      
      mydata.revised <- cbind(mydata.sig, com.fr)
      
      rm(com.match.hold)
    }
    
    #store results
    results.mntd [[e]] <- mydata.revised
    
  }
  
}


lengths(results.psv); lengths(results.mpd); lengths(results.mntd)

hist(lengths(results.psv),xlim=c(0,35))
hist(lengths(results.mpd), col="red", add=T)
hist(lengths(results.mntd), col="blue", add=T)

#obtain the means for each metric and grid size
mydata.all <- rep(NA, 7)
mydata.sesmpd <- rep(NA, 20)
mydata.sesmntd <- rep(NA, 20)
mydata.psv <- rep(NA, 20)
mydata.sesmpd2 <- rep(NA, 20)
mydata.sesmntd2 <- rep(NA, 20)
mydata.psv2 <- rep(NA, 20)

for (e in 1:length(mydata.filepath)){
  #for (e in 1){
  
  mydata <- read.csv(mydata.filepath[e], header=T, row.names=1)
  mydata <- mydata[,!(names(mydata) %in% ("com.fr"))]
  print(mydata[,(ncol(mydata)-6):ncol(mydata)])
  #mydata.psv [e] <- length(which(mydata[,(ncol(mydata)-2)]<0.85))
  #mydata.psv2 [e] <- length(which(mydata[,(ncol(mydata)-2)]>0.85))
  
  mydata.all <- rbind(mydata[nrow(mydata),(ncol(mydata)-6):ncol(mydata)], mydata.all)
  
}

#remove NA row
mydata.all <- mydata.all[-nrow(mydata.all),]
#flip the row order so that grid size increases from top to bottom
mydata.all <- mydata.all[ nrow(mydata.all):1, ]
#label row names with grid size
rownames(mydata.all) <- seq(0.05,1.00,by=0.05)

#test correlation with metrics and gridsizes
psv.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=mydata.all$psv, method = "pearson")
ses.mpd.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=mydata.all$ses.mpd, method = "pearson")
ses.mntd.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=mydata.all$ses.mntd, method = "pearson")

#test correlation with metrics and mean number of species per community
psv.coms.co <- cor.test(x=mycoms.meanrichness, y=mydata.all$psv, method = "pearson")
ses.mpd.coms.co <- cor.test(x=mycoms.meanrichness, y=mydata.all$ses.mpd, method = "pearson")
ses.mntd.coms.co <- cor.test(x=mycoms.meanrichness, y=mydata.all$ses.mntd, method = "pearson")


#################################################################################################################

#split results by PAIC 
splitbyPAIC <- function(mydata, splitrule, colofinterest="com.fr"){
  split.results<- as.list(rep(NA,length(mydata)))
  length(split.results) <- length(mydata)
  names(split.results) <- names(mydata)
  
  for (f in 1: length(mydata)){
    ind.com <- mydata[[f]]
    paic.holder <- which(colnames(ind.com)==colofinterest)
    
    com.holder <- ind.com[ind.com[,paic.holder]==splitrule,]
    
    if(length(nrow(com.holder) < 1) == 0){
      split.results[[f]] <- NA
    }
    else if(length(nrow(com.holder) > 0) > 0){
      com.holder <- ind.com[ind.com[,paic.holder]==splitrule,]
      split.results[[f]] <- com.holder
    }
    
  }
  return(split.results)
}

splitbyPAIC(results.mpd,"M")
splitbyPAIC(results.mntd,"M")
splitbyPAIC(results.psv,"M")

splitbyPAIC(results.mpd,"L")
splitbyPAIC(results.mntd,"L")
splitbyPAIC(results.psv,"L")

splitbyPAIC(results.psv,"V")

#################################################################################################################
library(viridis)

mypal_viridis <- viridis(20,1,0,1,"D")
mypal_viridis_t <- viridis(20,0.5,0,1,"D")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic/")
#individual PAICs
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.l/")
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.m/")


#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

#genetic
mydata.filepathshort[-grep("phylo_fr", mydata.filepathshort)] -> mydata.filepathshort
mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

#objects to store results
results.gen <- list()
length(results.gen) <- as.numeric(length(mydata.filepathshort))
names(results.gen) <- seq(from=0.05,to=1.00,by=0.05)

#loop to read in all grid size objects to a list
for (a in 1: length(mydata.filepathshort)){
  mydata.hold <- read.csv(mydata.filepath[a], header=T, row.names=1)
  results.gen [[a]] <- mydata.hold[colnames(mydata.hold) %in% c("ses.mpd",	"ses.mpd.p",	"ses.mntd",	"ses.mntd.p",	"psv",	"psv.p",	"n.coms")]
  
}

results.gen.means <- as.data.frame(matrix(nrow=length(results.gen),ncol=length(c("ses.mpd",	"ses.mpd.p",	"ses.mntd",	"ses.mntd.p",	"psv",	"psv.p",	"n.coms"))))
colnames(results.gen.means) <- c("ses.mpd",	"ses.mpd.p",	"ses.mntd",	"ses.mntd.p",	"psv",	"psv.p",	"n.coms")
#loop to get the means out of all objects in results.gen
for (b in 1:length(results.gen)){
  results.gen.means[b,] <- results.gen[[b]][as.numeric(nrow(results.gen[[b]])),]
  print(results.gen[[b]])
  results.gen[[b]] <- results.gen[[b]][-as.numeric(nrow(results.gen[[b]])),]
}

#if results are for a PAIC, assign to new object for plotting
results.gen.l <- results.gen
results.gen.means.l <- results.gen.means

results.gen.m <- results.gen
results.gen.means.m <- results.gen.means


#plot standard stacked MPD vs MNTD
plot(seq(from=0.05,to=1.00,by=0.05), results.gen.means$ses.mpd, pch=16, col=rgb(0,114/255,178/255,1),  ylim=c(-6,6), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Metric Value")
for(c in 1:length(results.gen)){
  boxplot(results.gen[[c]]$ses.mpd, border=rgb(0,114/255,178/255,1), boxwex=0.075, at=(seq(from=0.05,to=1.00,by=0.05)[c]), add=TRUE, axes=FALSE,col=mypal_viridis_t[c])
  
  points(x=rep(as.numeric(names(results.gen)[c]), times=length(results.gen[[c]]$ses.mpd[results.gen[[c]]$ses.mpd.p >= 0.975 | results.gen[[c]]$ses.mpd.p <= 0.025])),
         y=results.gen[[c]]$ses.mpd[results.gen[[c]]$ses.mpd.p >= 0.975 | results.gen[[c]]$ses.mpd.p <= 0.025],
         col=rgb(0,114/255,178/255,1),
         pch=18, cex=1.75)
  
  #abline(h=results.tails$ses.mpd.lower[c], col="blue", lty=2)
  #abline(h=results.tails$ses.mpd.upper[c], col="blue", lty=2)
}
points(seq(from=0.05,to=1.00,by=0.05), results.gen.means$ses.mntd, pch=16, col=rgb(213/255,94/255,0,1),  ylim=c(-6,6), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Metric Value")
for(c in 1:length(results.gen)){
  boxplot(results.gen[[c]]$ses.mntd, border=rgb(213/255,94/255,0,1), boxwex=0.075, at=(seq(from=0.05,to=1.00,by=0.05)[c]), add=TRUE, axes=FALSE,col=mypal_viridis_t[c])
  
  points(x=rep(as.numeric(names(results.gen)[c]), times=length(results.gen[[c]]$ses.mntd[results.gen[[c]]$ses.mntd.p >= 0.975 | results.gen[[c]]$ses.mntd.p <= 0.025])),
         y=results.gen[[c]]$ses.mntd[results.gen[[c]]$ses.mntd.p >= 0.975 | results.gen[[c]]$ses.mntd.p <= 0.025],
         col=rgb(213/255,94/255,0,1),
         pch=18, cex=1.75)
}

abline(h=0, lty=2)
points(seq(from=0.05,to=1.00,by=0.05), results.gen.means$ses.mntd, pch=16, col=rgb(213/255,94/255,0,1),  ylim=c(-5,5), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Metric Value")
points(seq(from=0.05,to=1.00,by=0.05), results.gen.means$ses.mpd, pch=16, col=rgb(0,114/255,178/255,1),  ylim=c(-5,5), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Metric Value")
legend(x=0.8, y=4.75, legend=c("SES.MPD", "SES.MNTD"), col="black", fill=c(rgb(0,114/255,178/255,1), rgb(213/255,94/255,0,1)), bty="n")
for(c in 1:nrow(results.gen.means)){
  if(results.gen.means$ses.mpd.p[c] <= 0.025 | results.gen.means$ses.mpd.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c], results.gen.means$ses.mpd[c], pch=23, col="black", bg=rgb(0,114/255,178/255,1))
  }
  if(results.gen.means$ses.mntd.p[c]<= 0.025 | results.gen.means$ses.mntd.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c], results.gen.means$ses.mntd[c], pch=23, col="black", bg=rgb(213/255,94/255,0,1))
  }
}

#plot standard PSV
plot(seq(from=0.05,to=1.00,by=0.05), results.gen.means$psv, pch=16, col=rgb(0,158/255,115/255,1),  ylim=c(0,1.05), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="PSV")
for(c in 1:length(results.gen)){
  boxplot(results.gen[[c]]$psv, border=rgb(0,158/255,115/255,1), boxwex=0.075, at=(seq(from=0.05,to=1.00,by=0.05)[c]), add=TRUE, axes=FALSE,col=mypal_viridis_t[c])
  abline(h=0.854, lty=2)
  
  points(x=rep(as.numeric(names(results.gen)[c]), times=length(results.gen[[c]]$psv[results.gen[[c]]$psv.p >= 0.975 | results.gen[[c]]$psv.p <= 0.025])),
         y=results.gen[[c]]$psv[results.gen[[c]]$psv.p >= 0.975 | results.gen[[c]]$psv.p <= 0.025],
         col=rgb(0,158/255,115/255,1),
         pch=18, cex=1.75)
}
for(c in 1:nrow(results.gen.means)){
  if(results.gen.means$psv.p[c] <= 0.025 | results.gen.means$psv.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c], results.gen.means$psv[c], pch=23, col="black", bg=rgb(0,158/255,115/255,1))
  }
}


#plot PAICs against one another

#MPD
#first plot the mean by grid size points
#luzon
plot(seq(from=0.05,to=1.00,by=0.05)-0.01, results.gen.means.l$ses.mpd, pch=16, col=rgb(0,114/255,178/255,1),  ylim=c(-5,4), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Phylogenetic Dispersion")
abline(h=0, lty=2)
#mindanao
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.gen.means.m$ses.mpd, pch=16, col=rgb(0,158/255,115/255,1))
for(c in 1:length(results.gen.l)){
  boxplot(results.gen.l[[c]]$ses.mpd, border="black", lwd=2, boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]-0.01), add=TRUE, axes=FALSE, col=rgb(0,114/255,178/255,0.5), outcol=rgb(0,114/255,178/255,1))
  boxplot(results.gen.m[[c]]$ses.mpd, border="black", lwd=2, boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]+0.01), add=TRUE, axes=FALSE, col=rgb(0,158/255,115/255,0.5), outcol=rgb(0,158/255,115/255,1))

  points(x=rep(as.numeric(names(results.gen.l)[c])-0.01, times=length(results.gen.l[[c]]$ses.mpd[results.gen.l[[c]]$ses.mpd.p >= 0.975 | results.gen.l[[c]]$ses.mpd.p <= 0.025])),
         y=results.gen.l[[c]]$ses.mpd[results.gen.l[[c]]$ses.mpd.p >= 0.975 | results.gen.l[[c]]$ses.mpd.p <= 0.025],
         col="black",
           bg=rgb(0,114/255,178/255,1),
         pch=23, cex=1.5, lwd=2)
  points(x=rep(as.numeric(names(results.gen.m)[c])+0.01, times=length(results.gen.m[[c]]$ses.mpd[results.gen.m[[c]]$ses.mpd.p >= 0.975 | results.gen.m[[c]]$ses.mpd.p <= 0.025])),
         y=results.gen.m[[c]]$ses.mpd[results.gen.m[[c]]$ses.mpd.p >= 0.975 | results.gen.m[[c]]$ses.mpd.p <= 0.025],
         col="black",
           bg=rgb(0,158/255,115/255,1),
         pch=23, cex=1.5, lwd=2)
}

#MNTD on same figure
for(c in 1:length(results.gen.l)){
  boxplot(results.gen.l[[c]]$ses.mntd, border=rgb(0,114/255,178/255,0.75), lwd=2, boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]-0.01), add=TRUE, axes=FALSE, col=rgb(0,0,0,0.35), outcol=rgb(0,0,0,1))
  boxplot(results.gen.m[[c]]$ses.mntd, border=rgb(0,158/255,115/255,0.75), lwd=2, boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]+0.01), add=TRUE, axes=FALSE, col=rgb(0,0,0,0.35), outcol=rgb(0,0,0,1))
  
  points(x=rep(as.numeric(names(results.gen.l)[c])-0.01, times=length(results.gen.l[[c]]$ses.mntd[results.gen.l[[c]]$ses.mntd.p >= 0.975 | results.gen.l[[c]]$ses.mntd.p <= 0.025])),
         y=results.gen.l[[c]]$ses.mntd[results.gen.l[[c]]$ses.mntd.p >= 0.975 | results.gen.l[[c]]$ses.mntd.p <= 0.025],
         bg="black",
         col=rgb(0,114/255,178/255,1),
         pch=23, cex=1.5, lwd=2)
  points(x=rep(as.numeric(names(results.gen.m)[c])+0.01, times=length(results.gen.m[[c]]$ses.mntd[results.gen.m[[c]]$ses.mntd.p >= 0.975 | results.gen.m[[c]]$ses.mntd.p <= 0.025])),
         y=results.gen.m[[c]]$ses.mntd[results.gen.m[[c]]$ses.mntd.p >= 0.975 | results.gen.m[[c]]$ses.mntd.p <= 0.025],
         bg="black",
          col=rgb(0,158/255,115/255,1),
         pch=23, cex=1.5, lwd=2)
}

#MPD
points(seq(from=0.05,to=1.00,by=0.05)-0.01, results.gen.means.l$ses.mpd, pch=21, col="black", bg=rgb(0,114/255,178/255,1), lwd=2, cex=1.5)
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.gen.means.m$ses.mpd, pch=21, col="black", bg=rgb(0,158/255,115/255,1), lwd=2, cex=1.5)
#MNTD
points(seq(from=0.05,to=1.00,by=0.05)-0.01, results.gen.means.l$ses.mntd, pch=21, col=rgb(0,114/255,178/255,1), bg="black", lwd=2, cex=1.5)
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.gen.means.m$ses.mntd, pch=21, col=rgb(0,158/255,115/255,1), bg="black", lwd=2, cex=1.5)

#community means
for(c in 1:nrow(results.gen.means)){
  if(results.gen.means.l$ses.mpd.p[c] <= 0.025 | results.gen.means.l$ses.mpd.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c]-0.01, results.gen.means.l$ses.mpd[c], pch=23, col=rgb(1,1,1,0.65), bg=rgb(0,114/255,178/255,1), cex=1.25)
  }
  if(results.gen.means.m$ses.mpd.p[c] <= 0.025 | results.gen.means.m$ses.mpd.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c]+0.01, results.gen.means.m$ses.mpd[c], pch=23, col=rgb(1,1,1,0.65), bg=rgb(0,158/255,115/255,1), cex=1.25)
  }
  
  if(results.gen.means.l$ses.mntd.p[c] <= 0.025 | results.gen.means.l$ses.mntd.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c]-0.01, results.gen.means.l$ses.mntd[c], pch=23, bg=rgb(1,1,1,0.65), col=rgb(0,114/255,178/255,1), cex=1.25)
  }
  if(results.gen.means.m$ses.mntd.p[c] <= 0.025 | results.gen.means.m$ses.mntd.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c]+0.01, results.gen.means.m$ses.mntd[c], pch=23, bg=rgb(1,1,1,0.65), col=rgb(0,158/255,115/255,1), cex=1.25)
  }
  
}



legend(x="topright", 
       legend=c(expression('Luzon SES'[MPD]),expression('Luzon SES'[MNTD]),expression('Mindanao SES'[MPD]),expression('Mindanao SES'[MNTD])),
       pt.bg= c(rgb(0,114/255,178/255,1), "black", rgb(0,158/255,115/255,1),"black"), 
       col=c("black",rgb(0,114/255,178/255,1),"black", rgb(0,158/255,115/255,1)), 
       pt.lwd = 2,
       pt.cex = 2,
       pch=rep(22,4),
       bty="n", cex=0.75, y.intersp = 0.65)


#PSV
plot(seq(from=0.05,to=1.00,by=0.05)-0.01, results.gen.means.l$psv, pch=16, col=rgb(0,0,0,0.65),  ylim=c(0,1.05), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="PSV")
abline(h=0.854, lty=2)
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.gen.means.m$psv, pch=16, col=rgb(0,0,0,0.65))
for(c in 1:length(results.gen.l)){
  boxplot(results.gen.l[[c]]$psv, col=rgb(0,0,0,0.65), lwd=2, boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]-0.01), add=TRUE, axes=FALSE, border=rgb(0,114/255,178/255,0.5), outcol=rgb(0,114/255,178/255,1))
  boxplot(results.gen.m[[c]]$psv, col=rgb(0,0,0,0.65), lwd=2, boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]+0.01), add=TRUE, axes=FALSE, border=rgb(0,158/255,115/255,0.5), outcol=rgb(0,158/255,115/255,1))
  
  points(x=rep(as.numeric(names(results.gen.l)[c])-0.01, times=length(results.gen.l[[c]]$psv[results.gen.l[[c]]$psv.p >= 0.975 | results.gen.l[[c]]$psv.p <= 0.025])),
         y=results.gen.l[[c]]$psv[results.gen.l[[c]]$psv.p >= 0.975 | results.gen.l[[c]]$psv.p <= 0.025],
         bg=rgb(0,0,0,0.65),
         col=rgb(0,114/255,178/255,1),
         pch=23, cex=1.5, lwd=2)
  points(x=rep(as.numeric(names(results.gen.m)[c])+0.01, times=length(results.gen.m[[c]]$psv[results.gen.m[[c]]$psv.p >= 0.975 | results.gen.m[[c]]$psv.p <= 0.025])),
         y=results.gen.m[[c]]$psv[results.gen.m[[c]]$psv.p >= 0.975 | results.gen.m[[c]]$psv.p <= 0.025],
         bg=rgb(0,0,0,0.65),
         col=rgb(0,158/255,115/255,1),
         pch=23, cex=1.5, lwd=2)
}
points(seq(from=0.05,to=1.00,by=0.05)-0.01, results.gen.means.l$psv, pch=21, bg=rgb(0,0,0,0.65),col=rgb(0,114/255,178/255,1))
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.gen.means.m$psv, pch=21, bg=rgb(0,0,0,0.65), col=rgb(0,158/255,115/255,1))
for(c in 1:nrow(results.gen.means)){
  if(results.gen.means.l$psv.p[c] <= 0.025 | results.gen.means.l$psv.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c]-0.01, results.gen.means.l$psv[c], pch=23, col=rgb(1,1,1,0.65), bg=rgb(0,114/255,178/255,1), cex=1.25)
  }
  if(results.gen.means.m$psv.p[c] <= 0.025 | results.gen.means.m$psv.p[c] >= 0.975){
    points(seq(from=0.05,to=1.00,by=0.05)[c]+0.01, results.gen.means.m$psv[c], pch=23, col=rgb(1,1,1,0.65), bg=rgb(0,158/255,115/255,1), cex=1.25)
  }
  
}


legend(x="bottomleft", 
       legend=c('Luzon PSV','Mindanao PSV'),
       pt.bg= c(rgb(0,0,0,0.65), rgb(0,0,0,0.65)), 
       col=c(rgb(0,114/255,178/255,1), rgb(0,158/255,115/255,1)), 
       pt.lwd = 2,
       pt.cex = 2,
       pch=rep(22,4),
       bty="n", cex=0.75, y.intersp = 0.65)






#plot for standard vs new method?
plot(seq(0.05,1,0.05), results.gen.means$ses.mpd, ylim=c(-1,1), cex=1.25, ylab = "Metric Mean", xlab = "Grid Size", pch=21, bg=rgb(0,114/255,178/255,1))
points(seq(0.05,1,0.05), results.gen.means.s$ses.mpd, pch=23, cex=1.25, bg=rgb(0,114/255,178/255,1))
points(seq(0.05,1,0.05), results.gen.means$ses.mntd, bg=rgb(213/255,94/255,0,1), pch=21, cex=1.25)
points(seq(0.05,1,0.05), results.gen.means.s$ses.mntd, pch=23, bg=rgb(213/255,94/255,0,1), cex=1.25)
points(seq(0.05,1,0.05), results.gen.means.s$psv, bg=rgb(0,158/255,115/255,1), pch=23, cex=1.35)
points(seq(0.05,1,0.05), results.gen.means$psv, bg=rgb(0,158/255,115/255,1), pch=21, cex=1.20)

       legend=c("SES.MPD","SES.MPD.S", "SES.MNTD", "SES.MNTD.S", "PSV", "PSV.S"), 
       pt.bg=c(rgb(0,114/255,178/255,1),rgb(0,114/255,178/255,1),rgb(213/255,94/255,0,1),rgb(213/255,94/255,0,1),rgb(0,158/255,115/255,1), rgb(0,158/255,115/255,1)), 
       col=rep("black", 6),
       pch=c(21,23,21,23,21,23),
       cex=rep(1.25,6),
       bty="n")

plot(seq(0.05,1,0.05), results.gen.means$ses.mpd, ylim=c(-4,2), cex=1.25, ylab = "Metric Mean", xlab = "Grid Size", pch=21, bg=rgb(0,114/255,178/255,1))
points(seq(0.05,1,0.05), results.gen.means.s$ses.mpd, pch=23, cex=1.25, bg=rgb(0,114/255,178/255,1))
points(seq(0.05,1,0.05), results.gen.means$ses.mntd, bg=rgb(213/255,94/255,0,1), pch=21, cex=1.25)
points(seq(0.05,1,0.05), results.gen.means.s$ses.mntd, pch=23, bg=rgb(213/255,94/255,0,1), cex=1.25)
points(seq(0.05,1,0.05), results.gen.means.s$psv, bg=rgb(0,158/255,115/255,1), pch=23, cex=1.35)
points(seq(0.05,1,0.05), results.gen.means$psv, bg=rgb(0,158/255,115/255,1), pch=21, cex=1.20)

for(c in 1:length(results.gen)){
  points(rep(as.numeric(names(results.gen)[c])+0.01, times=length(results.gen[[c]]$ses.mpd)), results.gen[[c]]$ses.mpd, pch=1, cex=0.9, col=rgb(0,114/255,178/255,1))
  points(rep(as.numeric(names(results.gen)[c])+0.01, times=length(results.gen.s[[c]]$ses.mpd)), results.gen.s[[c]]$ses.mpd, pch=5, cex=0.9, col=rgb(0,114/255,178/255,1))

  points(rep(as.numeric(names(results.gen)[c])+0.01, times=length(results.gen[[c]]$ses.mntd)), results.gen[[c]]$ses.mntd, pch=1, cex=0.9, col=rgb(213/255,94/255,0,1))
  points(rep(as.numeric(names(results.gen))[c]+0.01, times=length(results.gen.s[[c]]$ses.mntd)), results.gen.s[[c]]$ses.mntd, pch=5, cex=0.9, col=rgb(213/255,94/255,0,1))

  points(rep(as.numeric(names(results.gen)[c])+0.01, times=length(results.gen[[c]]$psv)), results.gen[[c]]$psv, pch=1, cex=0.9, col=rgb(0,158/255,115/255,1))
  points(rep(as.numeric(names(results.gen))[c]+0.01, times=length(results.gen.s[[c]]$psv)), results.gen.s[[c]]$psv, pch=5, cex=0.9, col=rgb(0,158/255,115/255,1))
  
}

