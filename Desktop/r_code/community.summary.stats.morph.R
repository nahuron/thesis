#community summary stats MORPHOLOGY
#determine the mean of mean morphological disparity by community

morph.short <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL", pattern="[fr.csv]*$")
morph.long <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL", pattern="[fr.csv]*$", full.names = TRUE)


morph.short[-grep("_fr.csv", morph.short)] -> morph.short

morph.long[-grep("_fr.csv", morph.long)] -> morph.long

morph.comdist.means <- rep(NA, times=length(morph.short))

for (a in 1:length(morph.long)){
  read.csv(morph.long[a], header=T, row.names=1)->mymorph
  morph.comdist.means[a] <- tail(mymorph$emp.com.mean.holder, 1)
  
}


#sort out significant morphology communities

results.morph <- list()
length(results.morph) <- length(morph.short)

results.morphd <- list()
length(results.morph) <- length(morph.short)

#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)

#isolate the significant communities
for(z in 1:20){
  
  read.csv(morph.long[z], header=T, row.names=1)->mymorph
  if(colnames(mymorph)[ncol(mymorph)]=="com.fr"){
    mymorph.red <- mymorph[-ncol(mymorph),]
    mymorph.red <- mymorph.red[-nrow(mymorph.red),]
  }
  else{
    mymorph.red <- mymorph[-nrow(mymorph),]
    
  }
  mymorph.red[mymorph.red$emp.com.mean.pvalues>=0.975 | mymorph.red$emp.com.mean.pvalues<=0.025,] -> mymorph.sig
  mymorph.sig <- mymorph.sig[,colSums(mymorph.sig) !=0]
  mymorph.sig <- unique(mymorph.sig)
  names(results.morph)[z] <- morph.short[z]
  
  com.fr <- rep(NA, times=nrow(mymorph.sig))
  for (a in 1:(nrow(mymorph.sig))) { #isolate row of interest
    #for (a in 20)  {
    com.match.hold <- 0
    for (b in 1:(ncol(mymorph.sig)-2)) { #isolate cell within row of interest
      if (mymorph.sig[a,b]==1){
        m <- match(as.character(colnames(mymorph.sig[a,][b])), as.character(brach_fr_key[,1]))
        com.match.hold <- c(com.match.hold, as.character(brach_fr_key$FR.code[m]))
        com.match.hold <- com.match.hold[com.match.hold!=0]
        com.match.hold <- com.match.hold[!is.na(com.match.hold)]
        print(com.match.hold)
      }
    }
    com.match.hold <- unique(com.match.hold)
    print(com.match.hold)
    if (length(com.match.hold) > 1 ) {
      print(paste0("Community ", rownames(mymorph.sig)[a], " spans more than 1 Faunal Regions!"))
      print(com.match.hold)
      com.fr[a] <- as.character(length(com.match.hold))
    }
    else if (length(com.match.hold == 1)) {com.fr[a] <- as.character(brach_fr_key$FR.code[m])}
    
    mymorph.revised <- cbind(mymorph.sig, com.fr)
    
    
    rm(com.match.hold)
  }
  
  
  print(z)
  print(mymorph.sig)
  
  results.morph [[z]] <- mymorph.revised
  results.morphd [[z]] <- mymorph.sig[,-ncol(mymorph.sig)]
  results.morphd [[z]] <- results.morphd[[z]][,-ncol(results.morphd[[z]])]
}

#################################################################################################################
#####MORPHOLOGY VERSION
morph.short <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL", pattern="[fr.csv]*$")
morph.long <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL", pattern="[fr.csv]*$", full.names = TRUE)

results.d2 <- list()
length(results.d2) <- length(morph.long)

morph.short[--grep("_fr.csv", morph.short)] -> morph.short

morph.long[-grep("_fr.csv", morph.long)] -> morph.long


for (e in 1:length(mydata.filepathshort)){
  mydata <- read.csv(morph.long[e], header=T, row.names=1)
  names(results.d2)[e] <- paste0("D2_",mydata.filepathshort[e])
  
  if(colnames(mydata)[as.numeric(ncol(mydata))]=="com.fr"){
    mydata.red <- mydata[,!(names(mydata) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata[-as.numeric(nrow(mydata)),]
  }
  
  mydata.sig <- mydata.red[mydata.red$emp.com.mean.pvalues>=0.975 | mydata.red$emp.com.mean.pvalues<=0.025,]
  mydata.rm <- c("emp.com.mean.holder", "emp.com.mean.pvalues")
  mydata.sig <- mydata.sig[,!(names(mydata.sig) %in% mydata.rm)]
  mydata.sig <- mydata.sig[,colSums(mydata.sig) !=0]
  mydata.sig <- unique(mydata.sig)
  
  if(nrow(mydata.sig)>0){
    com.fr <- rep(NA, times=nrow(mydata.sig))
    for (a in 1:(nrow(mydata.sig))) { #isolate row of interest
      com.match.hold <- 0
      for (b in 1:(ncol(mydata.sig))) { #isolate cell within row of interest
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
      
      mydata.revised <- cbind(mydata.sig,mydata.red[rownames(mydata.red) %in% rownames(mydata.sig),c("emp.com.mean.holder",  "emp.com.mean.pvalues")], com.fr)
      
      
      rm(com.match.hold)
    }
    print(mydata.revised)
    results.d2 [[e]] <- mydata.revised
    
  }
}

results.d2 <- results.d2[lengths(results.d2)>0]


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

splitbyPAIC(results.d2,"M")
splitbyPAIC(results.d2,"L")


#################################################################################################################
library(viridis)

mypal_viridis <- viridis(20,1,0,1,"D")
mypal_viridis_t <- viridis(20,0.5,0,1,"D")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL/")

mydata.filepath.short <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL", pattern="[fr.csv]*$")
mydata.filepath <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL", pattern="[fr.csv]*$", full.names = TRUE)

mydata.filepath.short[-grep("_fr.csv", mydata.filepath.short)] -> mydata.filepath.short
mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

#objects to store results
results.morphd <- list()
length(results.morphd) <- as.numeric(length(mydata.filepathshort))
names(results.morphd) <- seq(from=0.05,to=1.00,by=0.05)



#loop to read in all grid size objects to a list
for (a in 1: length(mydata.filepathshort)){
  mydata.hold <- read.csv(mydata.filepath[a], header=T, row.names=1)
  results.morphd [[a]] <- mydata.hold[colnames(mydata.hold) %in% c("emp.com.mean.holder",	"emp.com.mean.pvalues")]
  
}

results.morphd.means <- as.data.frame(matrix(nrow=length(results.morphd),ncol=length(c("emp.com.mean.holder",	"emp.com.mean.pvalues"))))
colnames(results.morphd.means) <- c("emp.com.mean.holder",	"emp.com.mean.pvalues")
#loop to get the means out of all objects in results.morphd
for (b in 1:length(results.morphd)){
  results.morphd.means[b,] <- results.morphd[[b]][as.numeric(nrow(results.morphd[[b]])),]
  print(results.morphd[[b]])
  results.morphd[[b]] <- results.morphd[[b]][-as.numeric(nrow(results.morphd[[b]])),]
}

plot(seq(from=0.05,to=1.00,by=0.05), results.morphd.means$emp.com.mean.holder, pch=16, col=rgb(0,0,0,1),  ylim=c(-0.1,15), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Mahalanobis D^2 Distance")
for(c in 1:length(results.gen)){
  boxplot(results.morphd[[c]]$emp.com.mean.holder, border=rgb(0,0,0,1), boxwex=0.075, at=(seq(from=0.05,to=1.00,by=0.05)[c]), add=TRUE, axes=FALSE,col=mypal_viridis_t[c])
  
  points(x=rep(as.numeric(names(results.morphd)[c]), times=length(results.morphd[[c]]$emp.com.mean.holder[results.morphd[[c]]$emp.com.mean.pvalues >= 0.975 | results.morphd[[c]]$emp.com.mean.pvalues <= 0.025])),
         y=results.morphd[[c]]$emp.com.mean.holder[results.morphd[[c]]$emp.com.mean.pvalues >= 0.975 | results.morphd[[c]]$emp.com.mean.pvalues <= 0.025],
         col=rgb(0,0,0,1),
         pch=18, cex=1.75)
  abline(h=4,lty=2)
}

#plot grid size and mean morphological disparity
plot(seq(from=0.05,to=1.00,by=0.05), results.morphd.means$emp.com.mean.holder, pch=16, col="black", xlab="Grid Size (Decimal Degrees)", xlim=c(0,1.05), ylim=c(3,5), main="Mean Morphological Disparityper\nCommunity Membership Grid Size")

cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=results.morphd.means$emp.com.mean.holder, method = "pearson")
