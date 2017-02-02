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


