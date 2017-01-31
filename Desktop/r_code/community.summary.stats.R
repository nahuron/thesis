# find out the average community size across grid sizes
library(fields)

test <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv3", pattern="[fr.csv]*$")
test2 <- list.files(path="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv3", pattern="[fr.csv]*$", full.names=TRUE)

test[-grep("_fr.csv", test)] -> test.short
test2[-grep("_fr.csv", test2)] -> test.long

comlist <- list()
length(comlist) <- length(test.short)

#total species across communities
for (c in 1:length(test.long)){
  read.csv(test.long[c], header=T, row.names=1)->comlist[[c]]
  comlist[[c]] <- comlist[[c]][colnames(comlist[[c]])!="species2"]
  comlist[[c]] <- comlist[[c]][colnames(comlist[[c]])!="species3"]
  comlist[[c]] <- comlist[[c]][colnames(comlist[[c]])!="c.f._bonitae"]
}

comlister <- lapply(comlist, colSums)
comlister.num <- NA

for(a in 1:length(comlister)){
  comlister.num <- c(comlister.num, (40-length(which(comlister[[a]]==0))))
  print(comlister.num)
  comlister.num <- comlister.num[!is.na(comlister.num)]
  }

#obtain mean sp richness per community per grid size and number of communities per grid size
for (a in 1:length(test.long)){
  read.csv(test.long[a], header=T, row.names=1)->mycoms
  mycoms <- mycoms[colnames(mycoms)!="species2"]
  mycoms <- mycoms[colnames(mycoms)!="species3"]
  mycoms <- mycoms[colnames(mycoms)!="c.f._bonitae"]
  
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
    
    #sim.com.hold <- replicate(1000, as.matrix(com.simulator(2, 5, nrow(mycoms), sort(unique(brach_loc$Species)))))
    #simcoms.meanrichness[[a]] <- apply(X=sim.com.hold, MARGIN = 3, FUN = function(x) mean(rowSums(x)))
    
  }
  else{
    print(test.short[a])
    print(nrow(mycoms))
    mycoms.ncoms <- c(mycoms.ncoms, nrow(mycoms))
    
    print(rowSums(mycoms))
    
    print(mean(rowSums(mycoms))) 
    mycoms.meanrichness <- c(mycoms.meanrichness, mean(rowSums(mycoms)))
    
    #sim.com.hold <- replicate(1000, as.matrix(com.simulator(2, 5, nrow(mycoms), sort(unique(brach_loc$Species)))))
    #simcoms.meanrichness[[a]] <- apply(X=sim.com.hold, MARGIN = 3, FUN = function(x) mean(rowSums(x)))
    
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


##plot some results
#plot mean number of species per community for empirical and simulated communities
plot(seq(from=0.05,to=1.00,by=0.05), mycoms.meanrichness, pch=16, col="blue",  ylim=c(2,5), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Mean Number of Species Per Community", main="Mean Number of Species Per Community per\nCommunity Membership Grid Size")
bplot(simcoms.meanrichness, pos=seq(from=0.05,to=1.00,by=0.05), add=TRUE, axes=FALSE)

#plot number of communities per grid size
plot(seq(from=0.05,to=1.00,by=0.05), mycoms.ncoms, pch=16, col="black", xlim=c(0,1.05), ylim=c(30,40), ylab="Number of Communities", xlab="Grid Size (Decimal Degrees)", main="Number of Unique Communities per\nCommunity Membership Grid Size")

#plot total number of species per grid size
plot(seq(from=0.05,to=1.00,by=0.05), comlister.num, pch=16, col="blue", xlab="Grid Size (Decimal Degrees)", ylab="Total Number of Species", xlim=c(0,1.05), ylim=c(30,40), main="Number of Unique Species per\nCommunity Membership Grid Size")

#calculate Pearson's Product Moment Correlation Coefficient
mycoms.meanrichness.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=mycoms.meanrichness, method = "pearson")
simcoms.meanrichness.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=simcoms.meanrichness.means, method = "pearson")
mycoms.ncoms.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=mycoms.ncoms, method = "pearson")
mycoms.comlister.num.co <- cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=comlister.num, method = "pearson")
mycoms.meanrichness.co; simcoms.meanrichness.co; mycoms.ncoms.co; mycoms.comlister.num.co

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


####GENETIC DATA####
#set working directory
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

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
splitbyPAIC(results.mpd,"M")
splitbyPAIC(results.mntd,"M")
splitbyPAIC(results.psv,"M")

splitbyPAIC(results.d2,"L")
splitbyPAIC(results.mpd,"L")
splitbyPAIC(results.mntd,"L")
splitbyPAIC(results.psv,"L")

splitbyPAIC(results.psv,"V")

#################################################################################################################
library(viridis)

mypal_viridis <- viridis(20,1,0,1,"D")
mypal_viridis_t <- viridis(20,0.5,0,1,"D")

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

#genetic
mydata.filepathshort[-grep("phylo_fr", mydata.filepathshort)] -> mydata.filepathshort
mydata.filepath[-grep("_fr.csv", mydata.filepath)] -> mydata.filepath

#Tails
load(file="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/results.tails.rda")

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
  
  #abline(h=results.tails$ses.mntd.lower[c], col="red", lty=2)
  #abline(h=results.tails$ses.mntd.upper[c], col="red", lty=2)
}

abline(h=0, lty=2)
points(seq(from=0.05,to=1.00,by=0.05), results.gen.means$ses.mntd, pch=16, col=rgb(213/255,94/255,0,1),  ylim=c(-5,5), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Metric Value")
points(seq(from=0.05,to=1.00,by=0.05), results.gen.means$ses.mpd, pch=16, col=rgb(0,114/255,178/255,1),  ylim=c(-5,5), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Metric Value")
legend(x=0.8, y=4.75, legend=c("SES.MPD", "SES.MNTD"), col="black", fill=c(rgb(0,114/255,178/255,1), rgb(213/255,94/255,0,1)), bty="n")


plot(seq(from=0.05,to=1.00,by=0.05), results.gen.means$psv, pch=16, col=rgb(0,158/255,115/255,1),  ylim=c(0,1.05), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="PSV")
for(c in 1:length(results.gen)){
  boxplot(results.gen[[c]]$psv, border=rgb(0,158/255,115/255,1), boxwex=0.075, at=(seq(from=0.05,to=1.00,by=0.05)[c]), add=TRUE, axes=FALSE,col=mypal_viridis_t[c])
  abline(h=0.854, lty=2)
  
  points(x=rep(as.numeric(names(results.gen)[c]), times=length(results.gen[[c]]$psv[results.gen[[c]]$psv.p >= 0.975 | results.gen[[c]]$psv.p <= 0.025])),
       y=results.gen[[c]]$psv[results.gen[[c]]$psv.p >= 0.975 | results.gen[[c]]$psv.p <= 0.025],
       col=rgb(0,158/255,115/255,1),
       pch=18, cex=1.75)
  
  #abline(h=results.tails$psv.lower[c], col="black", lty=2)
  #abline(h=results.tails$psv.upper[c], col="black", lty=2)
}


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
  mydata.hold <- read.csv(list.files(folders[a], pattern="maxentResults.csv", full.names = TRUE), header=TRUE, row.names=1)
  print(folders.short[a])
  print(mydata.hold$Test.AUC)
  print(mydata.hold$Minimum.training.presence.test.omission)
  enms.auc[a]<-mydata.hold$Test.AUC[nrow(mydata.hold)]
  enms.om[a]<-mydata.hold$Minimum.training.presence.test.omission[nrow(mydata.hold)]
  enms.alt[a]<-mydata.hold$alt.contribution[nrow(mydata.hold)]
  enms.bio02[a]<-mydata.hold$bio_02.contribution[nrow(mydata.hold)]
  enms.bio03[a]<-mydata.hold$bio_03.contribution[nrow(mydata.hold)]
  enms.bio12[a]<-mydata.hold$bio_12.contribution[nrow(mydata.hold)]
  enms.bio15[a]<-mydata.hold$bio_15.contribution[nrow(mydata.hold)]
  
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
  mydata.hold <- read.csv(mydata.filepath[f], header=T, row.names=1)
  
  mydata.hold <- mydata.hold[colnames(mydata.hold)!="species2"]
  mydata.hold <- mydata.hold[colnames(mydata.hold)!="species3"]
  mydata.hold <- mydata.hold[colnames(mydata.hold)!="c.f._bonitae"]
  
  sp.com <- NA
  for (e in 1:nrow(mydata.hold)) {
    
    sp.com <- c(sp.com, colnames(mydata.hold[e,mydata.hold[e,]==1]))
    #print(sp.com)
  }
  #print(length(unique(sp.com[-is.na(sp.com)])))
  print(nrow(mydata.hold))
  #print(rowSums(mydata.hold))
  #print(mean(rowSums(mydata.hold)))
  com.spnum[f] <- length(unique(sp.com[-is.na(sp.com)]))
  sp.com <- NA
}


