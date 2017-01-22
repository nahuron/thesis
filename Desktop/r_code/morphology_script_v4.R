#Revised Morphology Script (MS Thesis)
#Nicholas Huron
#-----------------------------------------------------------------------------------------
library(FactoMineR)
library(pcaMethods)
library(corrplot)
library(Hmisc)

#set working space
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology")
#-----------------------------------------------------------------------------------------


#load in dataset 1 (here it is the morphometric data provided via CDS per individual with species identifier)
read.table("huron_2016_brachymeles_morph_v5_rawmeans.txt", header=T, fill=T, nrows=530) -> brach_morph

#Exclude any entries that are missing data
brach_morph <- na.omit(brach_morph)

#drop species that aren't considered in the CA structure portion
brach_morph <- brach_morph[!brach_morph[,1]=="species2",]
brach_morph <- brach_morph[!brach_morph[,1]=="species3",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f._bonitae",]
brach_morph <- brach_morph[!brach_morph[,1]=="vermis",]
brach_morph <- brach_morph[!brach_morph[,1]=="vindumi",]
brach_morph <- brach_morph[!brach_morph[,1]=="wrighti",]

#check number of entries that remain
length(brach_morph$Species)	#518 #here it is 515, why?

unique(brach_morph$Species) #45 levels, should be 42 since dalawangdaliri and suluensis are omitted; libayani dropped
brach_morph$Species <- droplevels(brach_morph$Species)

#now check to ensure levels were dropped
unique(brach_morph$Species) #45 levels, should be 43 since dalawangdaliri and suluensis are omitted; libayani dropped
#-----------------------------------------------------------------------------------------

#determine number of species in based on the number of unique values in the scaled dataset
spp.num <- length(unique(brach_morph[,"Species"]))
#create palette
sp.palette <- c("#DD422B","#78E63A","#626EDC","#2B352A","#5FD3C3","#E5B832","#D13D91","#D1937B","#48628D","#748F2F","#D1A7C9","#D44BD6","#DEDA97","#914A1A","#497972","#C3DC3C","#DC3F5C","#549F69","#501F21","#69E379","#8A325A","#6B3F82","#826A72","#3B2946","#385521","#D4CCC5","#9D8025","#7691D7","#D37EC0","#66E4AA","#9BA27D","#98B5D0","#D46F7C","#DA7A35","#BAE2BB","#ABDD7C","#CECD5F","#D4AB65","#704D2B","#A265CE","#993530","#4EBB3A","#7E7A46","#69C4D6","#4A942E")
#create set of colors for points
pts.pick <- sp.palette[1:spp.num]
names(pts.pick)<-unique(as.character(brach_morph[,"Species"]))

#store colors and species
species.legend <- cbind(unique(as.character(as.factor(brach_morph[,"Species"]))), pts.pick)
#-----------------------------------------------------------------------------------------

#alternative function
pca.alt <- PCA(brach_morph, scale.unit=T, ncp=(ncol(brach_morph)-1), graph=T, quali.sup=1)
plot.PCA(pca.alt, axes=c(1,2), choix="ind", col.quali=sp.palette, habillage="none", label="quali", cex=0.8)
#plotellipses(pca.alt, label="none", pch.means=16, pch=17)

summary(pca.alt)


#loadings<-sweep(pca.alt$var$coord,2,sqrt(pca.alt$eig[,1]),FUN="/")

#obtain all loadings regardless of community to establish approximate min and max for figure
pca.all <- prcomp(brach_morph[,2:18], scale=T)

plot(pca.all$x[,1], pca.all$x[,2], pch=16)

summary(pca.all)

plot(PC.values$MEANS[,1],PC.values$MEANS[,2], pch=16, xlab="PC 1", ylab= "PC 2", col=species.legend[,2],xlim=c(min(PC.values$ALL[,1]*1.25),max(PC.values$ALL[,1]*1.25)), ylim=c(min(PC.values$ALL[,2]*1.25),max(PC.values$ALL[,2]*1.25)), main=expression(Morphometric~PCA~of~Species~of~italic(Brachymeles)))
text(PC.values$MEANS[,1],PC.values$MEANS[,2],rownames(PC.values$MEANS), cex=0.6,pos=4?)
minor.tick(nx=5, ny=5, tick.ratio=0.5)
axis(side=1, at=seq(-10,7,1), labels=seq(-10,7,1)) 
#-----------------------------------------------------------------------------------------

#pull PC loadings and store in an object to make plotting easier

PC.storage <- function(pca_object, factor.names, subjectIDs)	{
		factor.names <- as.factor(factor.names)
		print(length(factor.names))
		#subjectIDs <- droplevels(as.factor(subjectIDs))
		#print(subjectIDs)

		PC.matrix <- matrix(data= NA, nrow= length(pca_object$x[,1]), ncol= length(pca_object$sdev))
		PC.mean.matrix <- matrix(data= NA, nrow= length(factor.names), ncol= length(pca_object$sdev))									
		colnames(PC.matrix) <- paste("PC",1:length(pca_object$sdev), sep="")
		rownames(PC.matrix) <- subjectIDs
		colnames(PC.mean.matrix) <- paste("PC",1:length(pca_object$sdev), sep="")
		rownames(PC.mean.matrix) <- factor.names
		
		#obtain matrix of all individual coordinates
		for (i in 1:length(pca_object$sdev))	{
			PC.matrix[,i] <- pca_object$x[,i]
	                                      	}
		#print(PC.matrix)
										
		for (j in 1:length(pca_object$sdev))	{ #for each of the PC's
			#PC.mean.holder <- as.vector(tapply(PC.matrix[,j],subjectIDs, mean, na.rm=1))
			PC.mean.holder <- tapply(PC.matrix[,j],subjectIDs, mean, na.rm=1)
			PC.mean.matrix[,j] <- PC.mean.holder
									                      	}
										
		PC.BOTH <- new.env()
		PC.BOTH$ALL <- PC.matrix
		PC.BOTH$MEANS <- PC.mean.matrix
		as.list(PC.BOTH)
		
											}
										
PC.storage(pca.all, as.character(rownames(species.legend)), brach_morph[,"Species"]) -> PC.values

#-----------------------------------------------------------------------------------------

#empirical communities in binary form
empir_coms <- read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brachymeles.communitymatrix.csv", header=T, sep=",", row.names=1)

#empir_coms <- empir_coms[-35,];empir_coms <- empir_coms[!colnames(empir_coms)=="species3"]



setClass(Class="com_distances", representation(sim.com.distances="list", emp.com.distances="list", sim.com.mdistances="vector", emp.com.mdistances="vector"))	#use @ notation to get at individual pieces


#-----------------------------------------------------------------------------------------
#takes empirical community matrix membership, species PCA results, and com.simulator function input
#and calculates intracommunity distances that are averaged to get by-community "disparity"
#which are then plotted by empirical vs. simulated communities

com.morph.disp <- function(empirical, simulated, spp.pca, plot=F)	{
	
	#prepare object with PC1 and PC2 values for reference
	PC2D <- cbind(spp.pca[,1], spp.pca[,2])
	names(PC2D) <- rownames(spp.pca)
	
	distance.holder <- {}
	species.holder <- {}
	M <- numeric(1)
	
	#create list objects of the intracommunity distances and mean intracommunity distances
	sim.com.distances <- as.list(rep(NA, (nrow(simulated))))
	names(sim.com.distances) <- paste("SIM", rownames(simulated), sep=".")
	sim.com.mdistances <- as.list(rep(NA, (nrow(simulated))))
	names(sim.com.mdistances) <- paste("SIM", rownames(simulated), sep=".")

	emp.com.distances <- as.list(rep(NA, (nrow(empirical)+1)))
	names(emp.com.distances) <- paste("EMP", rownames(empirical), sep=".")
	emp.com.mdistances <- as.list(rep(NA, (nrow(empirical))))
	names(emp.com.mdistances) <- paste("EMP", rownames(empirical), sep=".")
	
	# loop to get the distances for simulated communities
	for (r in 1:nrow(simulated))	{
	
		species.holder <- names(simulated[r,simulated[r,]==1])
		distance.holder <- as.data.frame(matrix(numeric(3*length(species.holder)), ncol=3))
		colnames(distance.holder) <- c("PC 1", "PC 2", "Species")
			
			for (s in 1: nrow(distance.holder))	{
			
				distance.holder[s,3] <- species.holder[s]
				distance.holder[s,1:2] <- PC2D[rownames(PC2D)==species.holder[s]]
				
												}
												
			M <- dist(as.matrix(cbind(distance.holder[,1],distance.holder[,2])))
			M <- as.vector(M)
			sim.com.distances[[r]] <- M
			sim.com.mdistances[[r]] <- mean(M)
			distance.holder <- {}
			M <- {}									
		
											}
											
	# loop to get the distances for empirical communities
	for (t in 1:nrow(empirical))	{
	
		species.holder <- names(empirical[t,empirical[t,]==1])
		distance.holder <- as.data.frame(matrix(numeric(3*length(species.holder)), ncol=3))
		colnames(distance.holder) <- c("PC 1", "PC 2", "Species")
			
			for (u in 1: nrow(distance.holder))	{
			
				distance.holder[u,3] <- species.holder[s]
				distance.holder[u,1:2] <- PC2D[rownames(PC2D)==species.holder[u]]

												}
												
			M <- dist(as.matrix(cbind(distance.holder[,1],distance.holder[,2])))
			M <- as.vector(M)
			emp.com.distances[[t]] <- M
			emp.com.mdistances[[t]] <- mean(M)
			distance.holder <- {}
			M <- {}									
	
											}
	
	if(plot==T) {
	#dev.new()
	#hist(as.numeric(sim.com.mdistances), xlim=c(0,12), col="blue", probability=T)
	#dev.new()
	hist(as.numeric(emp.com.mdistances), col="red", xlim=c(0,12), probability = T, ylim=c(0,0.3))
	lines(density(as.numeric(sim.com.mdistances)), col="blue")	
	
	#dev.new()
	plot(density(as.numeric(sim.com.mdistances)), col="blue")
	rug(as.numeric(emp.com.mdistances), col = "red", ticksize = 0.20)
	lines(density(as.numeric(emp.com.mdistances)),col="red")										
	            }										
	
	invisible(new("com_distances", sim.com.distances=sim.com.distances, emp.com.distances=emp.com.distances, sim.com.mdistances=as.vector(sim.com.mdistances, mode="numeric"), emp.com.mdistances=as.vector(emp.com.mdistances, mode="numeric")))													
															}
															
com.morph.disp(empir_coms, null.coms, PC.values$MEANS, plot=F)->answer


#need to develop methods to obtain values for comparison against the null distribution

emp.com.files <- list.files(path = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communities", pattern = ".csv$", full.names = TRUE)
emp.com.files.short <- gsub("\\.csv$","",list.files(path = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communities", pattern = ".csv$", full.names = FALSE))

for(a in 1:length(emp.com.files)){
#for(a in 1){
  #read csv file for emp_com
  emp.comm <- read.csv(file=paste0(emp.com.files[a]), header=T, row.names=1)
  emp.com.short <- gsub("\\.csv$","",emp.comm)
  
  #edit the community matrix
  print(paste0("There are ", nrow(emp.comm), " original empirical communities"))
  emp.comm <- emp.comm[colnames(emp.comm)!="species2"]
  emp.comm <- emp.comm[colnames(emp.comm)!="dalawangdaliri"]
  emp.comm <- emp.comm[colnames(emp.comm)!="c.f..bonitae"]
  emp.comm <- emp.comm[colnames(emp.comm)!="vermis"]
  emp.comm <- emp.comm[colnames(emp.comm)!="vindumi"]
  emp.comm <- emp.comm[colnames(emp.comm)!="wrighti"]
  emp.comm <- emp.comm[colnames(emp.comm)!="suluensis"]
  emp.comm <- emp.comm[colnames(emp.comm)!="libayani"]
  emp.comm <- emp.comm[rowSums(emp.comm)>=2,]
  print(paste0("There are ", nrow(emp.comm), " empirical communities"))
    
  #com morph function
  com.morph.holder <- com.morph.disp(emp.comm, null.coms, PC.values$MEANS, plot=F)
  
  #store the mean emp.com morphological distances
  emp.com.mean.holder <- com.morph.holder@emp.com.mdistances
  
  #create object for p-values per each emp.com
  emp.com.mean.pvalues <- rep(NA, length(emp.com.mean.holder))
  
  for (b in 1:nrow(emp.comm)){
    #(length(answer@sim.com.mdistances[answer@sim.com.mdistances<=answer@emp.com.mdistances[1]])/length(answer@sim.com.mdistances))
    emp.com.mean.pvalues[b] <- 1-mean(com.morph.holder@sim.com.mdistances > emp.com.mean.holder[b])
  }
  emp.com.newdata <- cbind(emp.com.mean.holder, emp.com.mean.pvalues)
  emp.com.newdata <- cbind(emp.comm, emp.com.newdata)
  print(emp.com.newdata)
  
  #write to new csv in the new morphology folder
  
  write.csv(emp.com.newdata, paste0("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph/", emp.com.files.short[a], "_morph.csv"))
  
  rm(list=c("emp.com.newdata", "emp.com.mean.pvalues", "emp.com.mean.holder", "com.morph.holder", "emp.comm"))
}


