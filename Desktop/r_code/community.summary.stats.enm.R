####ENM DATA####
#set working directory
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/com.NET/")

#obtain file paths
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))


#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)

library(viridis)

mypal_viridis <- viridis(20,1,0,1,"D")
mypal_viridis_t <- viridis(20,0.5,0,1,"D")

#objects to store results
results.enm <- list()
length(results.enm) <- as.numeric(length(mydata.filepathshort))
names(results.enm) <- seq(from=0.05,to=1.00,by=0.05)

#objects to store results for Luzon PAIC
results.enm.l <- list()
length(results.enm.l) <- as.numeric(length(mydata.filepathshort))
names(results.enm.l) <- seq(from=0.05,to=1.00,by=0.05)
#objects to store results for Mindanao PAIC
results.enm.m <- list()
length(results.enm.m) <- as.numeric(length(mydata.filepathshort))
names(results.enm.m) <- seq(from=0.05,to=1.00,by=0.05)

for (a in 1: length(mydata.filepathshort)){
  mydata.hold <- read.csv(mydata.filepath[a], header=T, row.names=1)
  results.enm [[a]] <- mydata.hold[colnames(mydata.hold) %in% c("n.coms",	"NETS",	"com.fr")]
  results.enm.l[[a]] <- results.enm[[a]][results.enm[[a]]$com.fr %in% "L",]
  results.enm.m[[a]] <- results.enm[[a]][results.enm[[a]]$com.fr %in% "M",]
}

#object to store means
results.enm.means <- as.data.frame(matrix(nrow=length(results.enm),ncol=length(c("n.coms",	"NETS",	"com.fr"))))
colnames(results.enm.means) <- c("n.coms",	"NETS",	"com.fr")

#object to store Luzon PAIC means
results.enm.means.l <- as.data.frame(matrix(nrow=length(results.enm.l),ncol=length(c("n.coms",	"NETS",	"com.fr"))))
colnames(results.enm.means.l) <- c("n.coms",	"NETS",	"com.fr")
#object to store Mindanao PAIC means
results.enm.means.m <- as.data.frame(matrix(nrow=length(results.enm.m),ncol=length(c("n.coms",	"NETS",	"com.fr"))))
colnames(results.enm.means.m) <- c("n.coms",	"NETS",	"com.fr")


#loop to get the means out of all objects in results.enm
for (b in 1:length(results.enm)){
  results.enm.means[b,] <- results.enm[[b]][as.numeric(nrow(results.enm[[b]])),]
  print(results.enm[[b]])
  results.enm[[b]] <- results.enm[[b]][-as.numeric(nrow(results.enm[[b]])),]
  
  #separate PAICs
  results.enm.means.l[b,]$n.coms <- mean(results.enm.l[[b]]$n.coms, na.rm = T)
  results.enm.means.l[b,]$NETS <- mean(results.enm.l[[b]]$NETS, na.rm = T)
  
  results.enm.means.m[b,]$n.coms <- mean(results.enm.m[[b]]$n.coms, na.rm = T)
  results.enm.means.m[b,]$NETS <- mean(results.enm.m[[b]]$NETS, na.rm = T)
  
}

plot(seq(from=0.05,to=1.00,by=0.05), results.enm.means$NETS, pch=16, col="black",  ylim=c(0,100), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Percent of Equivalent ENMs per Community")
for(c in 1:length(results.enm)){
  boxplot(results.enm[[c]]$NETS, border="black", boxwex=0.075, at=(seq(from=0.05,to=1.00,by=0.05)[c]), add=TRUE, axes=FALSE,col=mypal_viridis_t[c])
}
points(seq(from=0.05,to=1.00,by=0.05), results.enm.means$NETS, pch=16, col="black")


plot(seq(from=0.05,to=1.00,by=0.05)-0.01, results.enm.means.l$NETS, pch=16, col=rgb(0,114/255,178/255,1),  ylim=c(0,100), xlim=c(0.00,1.05), xlab="Grid Size (Decimal Degrees)", ylab="Percent of Equivalent ENMs per Community")
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.enm.means.m$NETS, pch=16, col=rgb(0,158/255,115/255,1))
for(c in 1:length(results.enm)){
  boxplot(results.enm.l[[c]]$NETS, border="black", boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]-0.01), add=TRUE, axes=FALSE, col=rgb(0,114/255,178/255,0.5), outcol=rgb(0,114/255,178/255,1), lwd=2)
  boxplot(results.enm.m[[c]]$NETS, border="black", boxwex=0.075/2, at=(seq(from=0.05,to=1.00,by=0.05)[c]+0.01), add=TRUE, axes=FALSE, col=rgb(0,158/255,115/255,0.5), outcol=rgb(0,158/255,115/255,1), lwd=2)
  
}
points(seq(from=0.05,to=1.00,by=0.05)-0.01, results.enm.means.l$NETS, pch=21, col="black", bg=rgb(0,114/255,178/255,1), lwd=2, cex=1.5)
points(seq(from=0.05,to=1.00,by=0.05)+0.01, results.enm.means.m$NETS, pch=21, col="black", bg=rgb(0,158/255,115/255,1), lwd=2, cex=1.5)

legend(x=0, y=15, 
       legend=c('Luzon PSV','Mindanao PSV'),
       col= c(rgb(0,0,0,1), rgb(0,0,0,1)), 
       pt.bg=c(rgb(0,114/255,178/255,1), rgb(0,158/255,115/255,1)), 
       pt.lwd = 2,
       pt.cex = 2,
       pch=rep(22,4),
       bty="n", cex=0.75, y.intersp = 0.65)


##############################################################
#plot grid size and mean percent of species with ENM overlap
plot(seq(from=0.05,to=1.00,by=0.05), results.enm.means$NETS, pch=16, col="black", xlab="Grid Size (Decimal Degrees)", xlim=c(0,1.05), ylim=c(0,100), main="Mean ENM Overlap per\nCommunity Membership Grid Size")

cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=results.enm.means$NETS, method = "pearson")
cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=results.enm.means.l$NETS, method = "pearson")
cor.test(x=seq(from=0.05,to=1.00,by=0.05), y=results.enm.means.m$NETS, method = "pearson")

