# find out the average community size across grid sizes
library(fields)


test <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern="[fr.csv]*$")
test2 <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern="[fr.csv]*$", full.names=TRUE)

test[-grep("_fr.csv", test)] -> test.short

test2[-grep("_fr.csv", test2)] -> test.long


for (a in 1:length(test.long)){
  read.csv(test.long[a], header=T, row.names=1)->mycoms
  
  if(a==1){
    mycoms.ncoms <- NA
    mycoms.meanrichness <- NA
    simcoms.meanrichness <- list()
    length(simcoms.meanrichness) <- length(test.long)
    names(simcoms.meanrichness) <- seq(from=0.05,to=1.00,by=0.05)
    
    print(test.short[a])
    print(nrow(mycoms))
    mycoms.ncoms <- c(mycoms.ncoms, nrow(mycoms))
    mycoms.ncoms <- mycoms.ncoms[-is.na(mycoms.ncoms)]
    
    print(rowSums(mycoms))

    print(mean(rowSums(mycoms)))
    mycoms.meanrichness <- c(mycoms.meanrichness, mean(rowSums(mycoms)))
    mycoms.meanrichness <- mycoms.meanrichness[-is.na(mycoms.meanrichness)]
    
    sim.com.hold <- replicate(1000, as.matrix(com.simulator(2, 5, nrow(mycoms), sort(unique(brach_loc$Species)))))
    simcoms.meanrichness[[a]] <- apply(X=sim.com.hold, MARGIN = 3, FUN = function(x) mean(rowSums(x)))
    
  }
  else{
    print(test.short[a])
    print(nrow(mycoms))
    mycoms.ncoms <- c(mycoms.ncoms, nrow(mycoms))
    
    print(rowSums(mycoms))
    
    print(mean(rowSums(mycoms))) 
    mycoms.meanrichness <- c(mycoms.meanrichness, mean(rowSums(mycoms)))
    
    sim.com.hold <- replicate(1000, as.matrix(com.simulator(2, 5, nrow(mycoms), sort(unique(brach_loc$Species)))))
    simcoms.meanrichness[[a]] <- apply(X=sim.com.hold, MARGIN = 3, FUN = function(x) mean(rowSums(x)))
    
  }
}

simcoms.meanrichness.gen <- simcoms.meanrichness

#simcom means of species richness

simcoms.meanrichness.means <- rep(NA, times=length(simcoms.meanrichness))
simcoms.meanrichness.sd <- rep(NA, times=length(simcoms.meanrichness))

for(b in 1: length(simcoms.meanrichness)){
  print(mean(simcoms.meanrichness[[b]]))
  print(sd(simcoms.meanrichness[[b]]))
  simcoms.meanrichness.means[b] <- mean(simcoms.meanrichness[[b]])
  simcoms.meanrichness.sd[b] <- sd(simcoms.meanrichness[[b]])
}

#determine the mean of mean morphological disparity by community

morph.short <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph", pattern="[fr.csv]*$")
morph.long <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph", pattern="[fr.csv]*$", full.names = TRUE)


morph.short[-grep("_fr.csv", morph.short)] -> morph.short

morph.long[-grep("_fr.csv", morph.long)] -> morph.long

morph.comdist.means <- rep(NA, times=length(morph.short))

for (a in 1:length(morph.long)){
  read.csv(morph.long[a], header=T, row.names=1)->mymorph
  morph.comdist.means[a] <- tail(mymorph$emp.com.mean.holder, 1)
  
}


#plot some results
plot(seq(from=0.05,to=1.00,by=0.05), mycoms.meanrichness, pch=16, col="blue",  ylim=c(2,5), xlim=c(0.00,1.05))
bplot(simcoms.meanrichness, pos=seq(from=0.05,to=1.00,by=0.05), add=TRUE, axes=FALSE)

plot(seq(from=0.05,to=1.00,by=0.05), mycoms.ncoms, pch=17, col="red", xlim=c(0,1.05), ylim=c(30,45), ylab="Number of Communities", xlab="Grid Size (Decimal Degrees)")

results.morph <- list()
length(results.morph) <- length(morph.short)

#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)

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
}

####GENETIC DATA####
#set working directory
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

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
  
  mydata <- read.csv(mydata.filepath[e], header=T, row.names=1)
  names(results.psv)[e] <- paste0("PSV_",mydata.filepathshort[e])
  
  if(colnames(mydata)[as.numeric(ncol(mydata))]=="com.fr"){
    mydata.red <- mydata[,!(names(mydata) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata[-as.numeric(nrow(mydata)),]
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
      for (b in 1:(ncol(mydata.sig)-2)) { #isolate cell within row of interest
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
    
    results.psv [[e]] <- mydata.revised
  }
}

#MPD
for (e in 1:length(mydata.filepathshort)){
  
  mydata <- read.csv(mydata.filepath[e], header=T, row.names=1)
  names(results.mpd)[e] <- paste0("MPD_",mydata.filepathshort[e])
  
  if(colnames(mydata)[as.numeric(ncol(mydata))]=="com.fr"){
    mydata.red <- mydata[,!(names(mydata) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata[-as.numeric(nrow(mydata)),]
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
      for (b in 1:(ncol(mydata.sig)-2)) { #isolate cell within row of interest
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
  
  mydata <- read.csv(mydata.filepath[e], header=T, row.names=1)
  names(results.mntd)[e] <- paste0("MNTD_",mydata.filepathshort[e])
  
  if(colnames(mydata)[as.numeric(ncol(mydata))]=="com.fr"){
    mydata.red <- mydata[,!(names(mydata) %in% c("com.fr", "n.coms"))]
    mydata.red <- mydata.red[-as.numeric(nrow(mydata.red)),]
  }
  else{
    mydata.red <- mydata[-as.numeric(nrow(mydata)),]
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
      for (b in 1:(ncol(mydata.sig)-2)) { #isolate cell within row of interest
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


mydata.all <- rep(NA, 7)
mydata.sesmpd <- rep(NA, 20)
mydata.sesmntd <- rep(NA, 20)
mydata.psv <- rep(NA, 20)
mydata.sesmpd2 <- rep(NA, 20)
mydata.sesmntd2 <- rep(NA, 20)
mydata.psv2 <- rep(NA, 20)

for (e in 1:length(mydata.filepath)){
  
  mydata <- read.csv(mydata.filepath[e], header=T, row.names=1)
  mydata <- mydata[,!(names(mydata) %in% ("com.fr"))]
  #print(mydata[,(ncol(mydata)-6):ncol(mydata)])
  mydata.psv [e] <- length(which(mydata[,(ncol(mydata)-2)]<0.85))
  mydata.psv2 [e] <- length(which(mydata[,(ncol(mydata)-2)]>0.85))
  
  mydata.all <- rbind(mydata[nrow(mydata),(ncol(mydata)-6):ncol(mydata)], mydata.all)
  
}

mydata.all <- mydata.all[-nrow(mydata.all),]




####ENM DATA####
setwd("/Users/nicholashuron/Desktop/niche_overlap_tests/test2")

list.files(full.names = TRUE) -> folders
list.files() -> folders.short

folders <- folders[-1]
folders <- folders[-32]
folders <- folders[-18]
folders <- sort(folders)

folders.short <- folders.short[-1]
folders.short <- folders.short[-32]
folders.short <- folders.short[-18]
folders.short <- sort(folders.short)

enms.auc <- rep(NA, times=length(folders))
names(enms.auc) <- folders.short
enms.om <- rep(NA, times=length(folders))
names(enms.om) <- folders.short
enms.alt <- rep(NA, times=length(folders))
names(enms.alt) <- folders.short
enms.bio02 <- rep(NA, times=length(folders))
names(enms.bio02) <- folders.short
enms.bio03 <- rep(NA, times=length(folders))
names(enms.bio03) <- folders.short
enms.bio12 <- rep(NA, times=length(folders))
names(enms.bio12) <- folders.short
enms.bio15 <- rep(NA, times=length(folders))
names(enms.bio15) <- folders.short

for (a in 1:length(folders)){
  mydata <- read.csv(list.files(folders[a], pattern="maxentResults.csv", full.names = TRUE), header=TRUE, row.names=1)
  print(folders.short[a])
  print(mydata$Test.AUC)
  print(mydata$Minimum.training.presence.test.omission)
  enms.auc[a]<-mydata$Test.AUC[nrow(mydata)]
  enms.om[a]<-mydata$Minimum.training.presence.test.omission[nrow(mydata)]
  enms.alt[a]<-mydata$alt.contribution[nrow(mydata)]
  enms.bio02[a]<-mydata$bio_02.contribution[nrow(mydata)]
  enms.bio03[a]<-mydata$bio_03.contribution[nrow(mydata)]
  enms.bio12[a]<-mydata$bio_12.contribution[nrow(mydata)]
  enms.bio15[a]<-mydata$bio_15.contribution[nrow(mydata)]
  
}

tester <- cbind(enms.auc,enms.om, enms.alt, enms.bio02, enms.bio03, enms.bio12, enms.bio15)

rmer <- c("brevidactylus", "cobos", "dalawangdaliri", "elerae", "isangdaliri", "ligtas", "minimus", "pathfinderi", "suluensis", "tiboliorum", "vermis", "vindumi", "wrighti")

tester.reg <- tester[!rownames(tester) %in% rmer,]


###raw communities
# find out the average community size across grid sizes

mydata.filepathshort <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern="[fr.csv]*$")
mydata.filepath <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern="[fr.csv]*$", full.names=TRUE)

mydata.filepathshort[-grep("_fr.csv", mydata.filepathshort)] -> mydata.filepathshort

mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

com.spnum <- rep(NA,20)

for (f in 1:length(mydata.filepath)){
  mydata <- read.csv(mydata.filepath[f], header=T, row.names=1)
  
  mydata <- mydata[colnames(mydata)!="species2"]
  mydata <- mydata[colnames(mydata)!="species3"]
  mydata <- mydata[colnames(mydata)!="c.f._bonitae"]
  
  sp.com <- NA
  for (e in 1:nrow(mydata)) {
    
    sp.com <- c(sp.com, colnames(mydata[e,mydata[e,]==1]))
    #print(sp.com)
  }
  #print(length(unique(sp.com[-is.na(sp.com)])))
  print(nrow(mydata))
  #print(rowSums(mydata))
  #print(mean(rowSums(mydata)))
  com.spnum[f] <- length(unique(sp.com[-is.na(sp.com)]))
  sp.com <- NA
}



##genetic data
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

mydata.filepath <- list.files(path=getwd(), pattern = "\\_fr.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = "_fr.csv$", full.names = FALSE))

metrics <- as.data.frame(matrix(ncol=3, nrow=length(mydata.filepathshort)))


for (e in 1:length(mydata.filepathshort)){
  
  mydata <- read.csv(mydata.filepath[e], header=T, row.names=1)
  print(mydata[nrow(mydata),])
  cols <- c("psv.p", "ses.mpd.p","ses.mntd.p")
  colnames(metrics) <- cols
  print(mydata[nrow(mydata),cols])
  metrics[e,] <- mydata[nrow(mydata),cols]
  
}
