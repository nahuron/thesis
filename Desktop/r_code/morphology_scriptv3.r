#Revised Morphology Script (MS Thesis)
#Nicholas Huron
#-----------------------------------------------------------------------------------------
library(FactoMineR)
library(pcaMethods)
library(corrplot)

#set working space
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/brachymeles.morphology/")
#-----------------------------------------------------------------------------------------


#SKIP THIS SECTION IF YOU HAVE YOUR DATA TRANSFORMED!
#-----------------------------------------------------------------------------------------
#untransformed data
#read.table("huron_brachymeles_morph_raw.txt", header=T, fill=T, nrows=555) -> brachymeles

#transformed data, no means
read.table("huron_2015_brachymeles_morph.txt", header=T, fill=T, nrows=555) -> brachymeles

	brach.char.means <- as.data.frame(matrix(data=NA, ncol= ncol(brachymeles), nrow=length(unique(brachymeles$Species))))
	brach.char.means[,1] <- unique(brachymeles$Species)
	colnames(brach.char.means) <- colnames(brachymeles)
	brach.char.vars <- as.data.frame(matrix(data=NA, ncol= ncol(brachymeles), nrow=length(unique(brachymeles$Species))))
	brach.char.vars[,1] <- unique(brachymeles$Species)
	colnames(brach.char.vars) <- colnames(brachymeles)


for (z in 2:ncol(brachymeles))	{
	brach.char.means[,z] <- tapply(brachymeles[,z], brachymeles$Species, mean)
	brach.char.vars[,z] <- tapply(brachymeles[,z], brachymeles$Species, var)
	
											}
save.xlsx <- function (file, ...)
  {
      require(xlsx, quietly = TRUE)
      objects <- list(...)
      fargs <- as.list(match.call(expand.dots = TRUE))
      objnames <- as.character(fargs)[-c(1, 2)]
      nobjects <- length(objects)
      for (i in 1:nobjects) {
          if (i == 1)
              write.xlsx(objects[[i]], file, sheetName = objnames[i])
          else write.xlsx(objects[[i]], file, sheetName = objnames[i],
              append = TRUE)
      }
      print(paste("Workbook", file, "has", nobjects, "worksheets."))
}

save.xlsx("brachymeles.standardized.means.vars.xlsx", brach.char.means, brach.char.vars )
#-----------------------------------------------------------------------------------------




#load in dataset 1 (here it is the morphometric data provided via CDS per individual with species identifier)
read.table("huron_2016_brachymeles_morph_v4_rawmeans.txt", header=T, fill=T, nrows=530) -> brachymeles



#alternative morphometric data (no TL)
#read.table("huron_2016_brachymeles_morph_v3_rawmeans_noTL.txt", header=T, fill=T, nrows=530) -> brachymeles


#grab the variances and means for all of the characters



#Exclude any entries that are missing data
brachymeles <- na.omit(brachymeles)
#check number of entries that remain
length(brachymeles$Species)	#518

unique(brachymeles$Species) #45 levels, should be 42 since dalawangdaliri and suluensis are omitted; libayani dropped
brachymeles$Species <- droplevels(brachymeles$Species)
#now check to ensure levels were dropped
unique(brachymeles$Species) #45 levels, should be 43 since dalawangdaliri and suluensis are omitted; libayani dropped
#-----------------------------------------------------------------------------------------


#summary of object brachymeles
summary(brachymeles)
cor.mat <- round(cor(brachymeles[,2:18]),2)
dev.new()
corrplot(cor.mat, type="upper", order="hclust", 
         tl.col="black", tl.srt=45)
#-----------------------------------------------------------------------------------------

#determine number of species in based on the number of unique values in the scaled dataset
spp.num <- length(unique(brachymeles[,"Species"]))
#create palette
sp.palette <- c("#DD422B","#78E63A","#626EDC","#2B352A","#5FD3C3","#E5B832","#D13D91","#D1937B","#48628D","#748F2F","#D1A7C9","#D44BD6","#DEDA97","#914A1A","#497972","#C3DC3C","#DC3F5C","#549F69","#501F21","#69E379","#8A325A","#6B3F82","#826A72","#3B2946","#385521","#D4CCC5","#9D8025","#7691D7","#D37EC0","#66E4AA","#9BA27D","#98B5D0","#D46F7C","#DA7A35","#BAE2BB","#ABDD7C","#CECD5F","#D4AB65","#704D2B","#A265CE","#993530","#4EBB3A","#7E7A46","#69C4D6","#4A942E")
#create set of colors for points
pts.pick <- sp.palette[1:spp.num]
names(pts.pick)<-unique(as.character(brachymeles[,"Species"]))

#store colors and species
species.legend <- cbind(unique(as.character(as.factor(brachymeles[,"Species"]))), pts.pick)
#-----------------------------------------------------------------------------------------


#obtain all loadings regardless of community to establish approximate min and max for figure
pca.all <- prcomp(brachymeles[,2:18], scale=T)

#alternative function
PCA <- PCA(brachymeles, scale.unit=T, ncp=(ncol(brachymeles)-1), graph=F, quali.sup=1)
dev.new()
plot.PCA(PCA, axes=c(1,2), choix="ind", col.quali=sp.palette)

summary(pca.all)

#-----------------------------------------------------------------------------------------





#pull PC loadings and store in an object to make plotting easier

PC.storage <- function(PCA, factor.names, subjectIDs)	{
		factor.names <- as.factor(factor.names)
		subjectIDs <- droplevels(as.factor(subjectIDs))

		PC.matrix <- matrix(data= NA, nrow= length(PCA$x[,1]), ncol= length(PCA$sdev))
		PC.mean.matrix <- matrix(data= NA, nrow= length(factor.names), ncol= length(PCA$sdev))									
		colnames(PC.matrix) <- paste("PC",1:length(PCA$sdev), sep="")
		rownames(PC.matrix) <- subjectIDs
		colnames(PC.mean.matrix) <- paste("PC",1:length(PCA$sdev), sep="")
		rownames(PC.mean.matrix) <- factor.names
		
		for (i in 1:length(PCA$sdev))	{
			PC.matrix[,i] <- PCA$x[,i]
										}
										
		for (j in 1:length(PCA$sdev))	{
			PC.mean.holder <- as.vector(tapply(PC.matrix[,j],subjectIDs, mean, na.rm=1))
			PC.mean.matrix[,j] <- PC.mean.holder
										}
										
		PC.BOTH <- new.env()
		PC.BOTH$ALL <- PC.matrix
		PC.BOTH$MEANS <- PC.mean.matrix
		as.list(PC.BOTH)
		
											}
										
PC.storage(pca.all, species.legend[,1], brachymeles[,"Species"]) -> PC.values

#-----------------------------------------------------------------------------------------

dev.new()
plot(PC.values$MEANS[,1],PC.values$MEANS[,2], pch=16, xlab="PC 1", ylab= "PC 2", col=species.legend[,2],xlim=c(min(PC.values$ALL[,1]*1.25),max(PC.values$ALL[,1]*1.25)), ylim=c(min(PC.values$ALL[,2]*1.25),max(PC.values$ALL[,2]*1.25)), main=expression(Morphometric~PCA~of~Species~of~italic(Brachymeles)))
text(PC.values$MEANS[,1],PC.values$MEANS[,2],as.numeric(1:length(species.legend[,1])), cex=0.7,pos=4)

dev.new()
#Plot first community and store the results in the list
plot(PC.values$ALL[,1],PC.values$ALL[,2], col=brachymeles$Species, pch=16, xlab="PC 1", ylab= "PC2", main=expression(Morphometric~PCA~of~Species~of~italic(Brachymeles)))
text(PC.values$ALL[,1],PC.values$ALL[,2],as.numeric(brachymeles[,1]), cex=0.7,pos=4)
#legend(legend=as.numeric(as.factor(species.legend[,1])), pch=16, col=species.legend[,2], x="topright")

dev.new()
#Plot first community and store the results in the list
plot(PC.values$ALL[,1],PC.values$ALL[,2],xlim=c(min(PC.values$ALL[,1]*1.25),max(PC.values$ALL[,1]*1.25)), ylim=c(min(PC.values$ALL[,2]*1.25),max(PC.values$ALL[,2]*1.25)), col=brachymeles$Species, pch=16, xlab="PC 1", ylab= "PC2", main=expression(Morphometric~PCA~of~Species~of~italic(Brachymeles)))
text(PC.values$ALL[,1],PC.values$ALL[,2],as.numeric(brachymeles[,1]), cex=0.7,pos=4)
legend(legend=as.numeric(as.factor(species.legend[,1])), pch=16, col=species.legend[,2], x="topright")

#-----------------------------------------------------------------------------------------

#load in communities (community identifiers, species identifiers)
#communities <- read.table("community_assignments.txt", header=T)

#empirical communities in binary form
empir_coms <- read.table("brachymeles.communitymatrix_2016.csv", header=T, sep=",", row.names=1)
colnames(empir_coms)[9]<- "c.f._bonitae"
empir_coms <- empir_coms[!colnames(empir_coms)=="apus"];empir_coms <- empir_coms[!colnames(empir_coms)=="miriamae"];empir_coms <- empir_coms[!colnames(empir_coms)=="libayani"]
empir_coms <- empir_coms[empir_coms[,"dalawangdaliri"]!=1,]
#empir_coms <- empir_coms[-35,];empir_coms <- empir_coms[!colnames(empir_coms)=="species3"]



#-----------------------------------------------------------------------------------------

#Optional code to prepare a community membership matrix in binary form for future reference (looks like null community matrix)
communities<- communities[order(communities$Species),]
communities <- cbind(communities,order(communities[,1]))


communities[,ncol(communities)] <- letters[communities[,ncol(communities)]]
alphacommunities <- communities[order(communities[,ncol(communities)]),]
colnames(alphacommunities) <- c(colnames(communities[1:3]),"species")

binary.communities <- matrix(ncol=length(unique(alphacommunities$species)), nrow=length(unique(alphacommunities$Community)))

for (i in 1:length(unique(alphacommunities$Community)))	{
	for (j in 1:length(unique(alphacommunities$species)))	{
		
		if(alphacommunities[j,2]== i)	{
			binary.communities[i,j]<-1
										}
																}
														}
		binary.communities[is.na(binary.communities)] <- 0
		colnames(binary.communities) <- alphacommunities$species
		rownames(binary.communities) <- seq(1:nrow(binary.communities))
		colnames(binary.communities) <- communities$Species

#-----------------------------------------------------------------------------------------

#function that creates a list containing PCs for each species (grouped by community)
#input data should be a:
#	binary community matrix (species x community)
#	output object from PC.storage function, specifically the PC.BOTH$MEANS
#resulting object should be list

spp.by.com.charsorter <- function(communities, PC.results, means.as.data.frame=FALSE)	{
	sxc <<- vector("list", length(1:nrow(communities)))
	names(sxc) <<- as.character(rownames(communities))
	
	
	
	
	results.holder <- matrix(NA)
	
	if(means.as.data.frame==FALSE)	{
		results.holder <- as.data.frame(PC.results$MEANS)
		species.names <- rownames(PC.results$MEANS)
		for (k in 1:nrow(communities))	{
			for (l in 1: ncol(communities))	{
				if (communities[k,l]==1)		{
					spp.hold <- species.names[l]
					sxc[[k]] <<- c(sxc[[k]], results.holder[rownames(results.holder)==spp.hold,])
			
												}
											}	
									
										}
		
									}	
	
	if(means.as.data.frame==TRUE){
								 }
																						}

spp.by.com.charsorter(empir_coms, PC.values)
#-----------------------------------------------------------------------------------------

#new function that takes empirical community matrix membership, species PCA results, and com.simulator function input
#and calculates intracommunity distances that are averaged to get by-community "disparity"
#which are then plotted by empirical vs. simulated communities

com.morph.disp <- function(empirical, simulated, spp.pca)	{
	PC1 <- spp.pca[,1]
	names(PC1) <- names(spp.pca[,1])
	PC2 <- spp.pca[,2]
	names(PC2) <- names(spp.pca[,2])
	PC2D <- cbind(PC1, PC2)
	
	#initialize storage objects
	
	DISTANCES <- new.env()
	DISTANCES <- as.list(DISTANCES)
	length(DISTANCES) <- sum(c(nrow(simulated),nrow(empirical)))
	
	for (o in 1:length(DISTANCES))	{		#consider setting length of each larger element of the list so that placement is easier below
		DISTANCES[[o]] <- NA
									}
	
	names(DISTANCES) <- c(paste("SIM", rownames(simulated), sep="."),rownames(empirical))
	print(DISTANCES)
	distance.holder <<- NA
	species.holder <<- NA

	
	#transfer simulated community information to new list object
	
	for (m in 1:nrow(simulated))	{			
		DISTANCES[[m]] <- rep(NA, times=(2*sum(simulated[m,])))
		for (n in 1:length(simulated[m,]))	{
			if (simulated[m,n]==1)				{
				distance.holder <<- PC2D[rownames(PC2D)==colnames(simulated)[n],]
				species.holder <<- c(species.holder, colnames(simulated)[n])
				species.holder <<- na.omit(species.holder)
				DISTANCES [[m]] <- c(DISTANCES[[m]], distance.holder)					
												}
											}
		DISTANCES[[m]] <- na.omit(DISTANCES[[m]])
		attr(DISTANCES[[m]], 'na.action') <- NULL
		attr(DISTANCES[[m]], 'Species') <- species.holder
		species.holder <<- NA							
									}
	
	
	print(DISTANCES)
	
	
	#transfer the empirical community information to the same list object
	
	for (p in 1:nrow(empirical))	{			
		DISTANCES[[nrow(simulated)+p]] <- rep(NA, times=(2*sum(empirical[p,])))
		for (q in 1:length(empirical[p,]))	{
			if (empirical[p,q]==1)				{
				distance.holder <<- PC2D[rownames(PC2D)==colnames(empirical)[q],]
				species.holder <<- c(species.holder, colnames(empirical)[q])
				species.holder <<- na.omit(species.holder)
				DISTANCES [[nrow(simulated)+p]] <- c(DISTANCES[[nrow(simulated)+p]], distance.holder)					
												}
											}
		DISTANCES[[nrow(simulated)+p]] <- na.omit(DISTANCES[[nrow(simulated)+p]])
		attr(DISTANCES[[nrow(simulated)+p]], 'na.action') <- NULL
		attr(DISTANCES[[nrow(simulated)+p]], 'Species') <- species.holder
		species.holder <<- NA							
									}
									
	print(DISTANCES)								
									
									
	#segregate com1 values into columns
	#spp.holder <- as.data.frame(matrix(numeric(3*length(attr(DISTANCES[[1]], 'Species'))),ncol=3))
	#colnames(spp.holder) <- c("PC 1", "PC 2", "Species")
	
	
	#create list object of the intracommunity distances and mean intracommunity distances
	
	com.distances <- new.env()
	com.distances <- as.list(com.distances)
	length(com.distances) <- length(DISTANCES)
	com.distances.means <<- numeric(length(DISTANCES))
	M <- numeric(1)
	
	
	for (r in 1:length(DISTANCES))	{
		spp.holder <- as.data.frame(matrix(numeric(3*length(attr(DISTANCES[[r]], 'Species'))),ncol=3))
		colnames(spp.holder) <- c("PC 1", "PC 2", "Species")

						spp.holder [,3] <- attr(DISTANCES[[r]], 'Species')
						pc1.holder <- NA
						pc2.holder <- NA
						
						for (s in seq(from=1, to=length(DISTANCES[[r]]), by=2))	{
							pc1.holder <- c(pc1.holder, DISTANCES[[r]][s])
							pc1.holder <- na.omit(pc1.holder)
																				}
						spp.holder [,1] <- pc1.holder
						
						for (t in seq(from=2,to=length(DISTANCES[[r]]), by=2))	{
							pc2.holder <- c(pc2.holder, DISTANCES[[r]][t])
							pc2.holder <- na.omit(pc2.holder)
																				}
						
						spp.holder [,2] <- pc2.holder
						
						M <- dist(as.matrix(cbind(spp.holder[,1],spp.holder[,2])))
						M <- as.vector(M)
						com.distances[[r]] <- M
						com.distances.means[r] <<- mean(M) 

									}
						print(com.distances.means)			
	
	
	
	#seperate the mean community values by simulated vs. empirical
	
	#simulated means
	
	simcom.distances.means <<- com.distances.means[1:nrow(simulated)]
	
	#empirical means
	
	empcom.distances.means <<- com.distances.means[(1+nrow(simulated)):length(com.distances.means)]
	
	
	dev.new()
	hist(simcom.distances.means, xlim=c(0,12), col="blue")
	dev.new()
	hist(empcom.distances.means, col="red", xlim=c(0,12))
	
	dev.new()
	plot(density(simcom.distances.means), col="blue")
	lines(density(empcom.distances.means),col="red")
	

	DISTANCE <<- list(PCs=DISTANCES, AllDistances=com.distances, MeanDistances=com.distances.means)
	
	#output.all <- new.env()
	#output.all$PCs <- DISTANCE[[1]]
	#output.all$AllDistances <- DISTANCE[[2]]
	#output.all$MeanDistances <- DISTANCE[[3]]
	
	invisible(DISTANCE)

	
															}
															
com.morph.disp(empir_coms, null.coms, PC.values$MEANS)->answer
