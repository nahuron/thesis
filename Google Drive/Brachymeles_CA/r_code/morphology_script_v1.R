#Final R Project 
#Nicholas Huron
#-----------------------------------------------------------------------------------------
#Script to create PCA Figures for individual communities and calculate corresponding 
#average distances between the species PC averages in each community

#set working space
setwd("/Users/nicholashuron/Desktop/brachymeles.morphology/")
#-----------------------------------------------------------------------------------------


#load in dataset 1 (here it is the morphometric data provided via CDS per individual with species identifier)
read.table("brachymeles_reduced.txt", header=T, fill=T, nrows=555) -> brachymeles

#load in communities (community identifiers, species identifiers)
communities <- read.table("community_assignments.txt", header=T)
#-----------------------------------------------------------------------------------------

#prepare a community membership matrix in binary form for future reference (looks like null community matrix)
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
#-----------------------------------------------------------------------------------------
						

#Exclude any entries that are missing data
brachymeles <- na.exclude(brachymeles)
#check number of entries that remain
length(brachymeles$Species)	#550
#-----------------------------------------------------------------------------------------


#give each individual in the loaded dataset a unique identifier
individuals <- with(brachymeles, paste(Species,SPPSPEC,sep="."))

#add the individual identifiers in to the original file
brachymeles <- as.data.frame(cbind(individuals,brachymeles))

#summary of object brachymeles
summary(brachymeles)
#-----------------------------------------------------------------------------------------


#Sort out communities by species membership and place in new dataframes with community appropriate names
community.sorter <- function(X,Y)	{

#Use the merge() function to merge by species to create a composite table
#that only includes the species in the communities (Y) object
#but now includes community number as well
merge(X,Y, by="Species") -> both.comm.species

return(both.comm.species)
									}
community.sorter(brachymeles,communities) -> brach

summary(brach)
#find out remaining number of entries
length(brach[,1])	#230
#-----------------------------------------------------------------------------------------


#scale the continuous measurements to be used for a given community (scale function)
st.brach.meas <- as.data.frame(scale(brach[,4:18]))

#attach identifiers to standardized vars (species and community treated as factors in cols 1 & 2
st.brach.meas <- cbind(brach[,20],brach[,19],st.brach.meas)
as.numeric(st.brach.meas[,1]) -> st.brach.meas[,1]

#assign names for use in later script
names(c(brach[20],brach[19],brach[4:18])) ->N
names(st.brach.meas) <- N
#-----------------------------------------------------------------------------------------


#determine number of communities based on number of unique values in community assignment matrix
com.num <- length(unique(communities[,"Community"]))
#determine number of species in those communities based on the number of unique values in the community assignment matrix
spp.num <- length(unique(communities[,"Species"]))
#create set of colors for points
pts.pick <- rainbow(com.num)
#-----------------------------------------------------------------------------------------


#Create object to store the pca loadings for reference outside of the function
elements.list <- 1:com.num
pca.data <- as.list(rep(NA, length(1:com.num)))
names(pca.data) <- paste("Community", elements.list, sep=".")

#create matrix to store means of coordinates for each community
comspp.means<-matrix(numeric(4*spp.num),ncol=4)
colnames(comspp.means) <- c("Community","Species","X Means","Y Means")


#obtain all loadings regardless of community to establish approximate min and max for figure
pca.all <- prcomp(st.brach.meas[,3:17])

#initialize for saving created image
dev.new()
tiff(filename="FigureA.tiff", units="px", height= 4800, width=4800, compression="none", res=600)


#perform the first community PCA
com <- st.brach.meas[st.brach.meas[,"Community"]==1,]
pca <- prcomp(st.brach.meas[st.brach.meas[,"Community"]==1,3:17])


#Plot first community and store the results in the list
plot(pca$x[,1],pca$x[,2],xlim=c(min(pca.all$x[,1]*1.25),max(pca.all$x[,1]*1.25)), ylim=c(min(pca.all$x[,2]*1.25),max(pca.all$x[,2]*1.25)), pch=16, col=pts.pick[1], xlab="PC 1", ylab= "PC2", main=expression(Morphometric~PCA~of~Species~of~italic(Brachymeles)~By~Community))
text(pca$x[,1],pca$x[,2],com$Spp, cex=0.7,pos=4,col=pts.pick[1])
legend(paste("Community", elements.list, sep="."), pch=16, col=pts.pick, x="topright")
pca.data[[1]] <- pca


#perform rest of community pca's and plot
for (i in 2:com.num)	{
				
							pca <- prcomp(st.brach.meas[st.brach.meas[,"Community"]==i,3:17])
							com <- st.brach.meas[st.brach.meas[,"Community"]==i,]
							points(pca$x[,1], pca$x[,2], pch=16, col= pts.pick[i])
							text(pca$x[,1], pca$x[,2], com$Spp, cex=0.7, pos=4, col=pts.pick[i])
							pca.data[[i]] <- pca
							
							}
dev.off()

#display summary and list of stored PCA results 
#Individual PCA loadings can be pulled for later use by assigning a new object to one of the community elements in the list
summary(pca.data)
pca.data


#Calculate the means for all of the species in each community
for (i in 1:com.num)	{
				
							pca <- prcomp(st.brach.meas[st.brach.meas[,"Community"]==i,3:17])
							com <- st.brach.meas[st.brach.meas[,"Community"]==i,]
							pca.holder <- cbind(st.brach.meas[st.brach.meas[,"Community"]==i,1:2],pca$x[,1],pca$x[,2])
							
									for(j in as.numeric(levels(as.factor(pca.holder[,1]))))		{
											
																									comspp.means[j,1] <- i
																									comspp.means[j,2] <- j
																									comspp.means [j,3] <- mean(pca.holder[pca.holder[,1]==j,3])
																									comspp.means [j,4] <- mean(pca.holder[pca.holder[,1]==j,4])
											
																								}
								
							}
		
		
#initialize for saving created image
dev.new()
plot.new()
tiff(filename="FigureB.tiff", units="px", height= 4800, width=4800, compression="none", res=600)
							
#Plot means by community with species 1 to initialize the figure
plot(comspp.means[1,3],comspp.means[1,4], pch=16, xlab="PC 1", ylab= "PC2", main=expression(Morphometric~PCA~of~Means~of~Species~of~italic(Brachymeles)~By~Community), xlim=c(min(pca.all$x[,1]*1.25),max(pca.all$x[,1]*1.25)),ylim=c(min(pca.all$x[,2]*1.25), max(pca.all$x[,2]*1.25)), col=pts.pick[comspp.means[1,1]])
text(comspp.means[1,3],comspp.means[1,4],comspp.means[1,2], cex=0.7, pos=4, col=pts.pick[comspp.means[1,1]])
legend(paste("Community", elements.list, sep="."), pch=16, col=pts.pick, x="topright")

#Plot the remaining points
for (i in 2:spp.num)	{
							points(comspp.means[i,3],comspp.means[i,4], pch=16, col=pts.pick[comspp.means[i,1]])
							text(comspp.means[i,3],comspp.means[i,4],comspp.means[i,2], cex=0.7, pos=4, col=pts.pick[comspp.means[i,1]])
						}			

dev.off()

#calculate the distances between each species' mean and the other species' means and store in a dataframe



#initialize a list to store the distances
distances <- as.list(rep(NA, length(1:com.num)))

#Object specified for the means of the distances for each community
distances.means <- numeric(length(1:com.num))
names(distances.means) <- paste("Community", elements.list, sep=".")

#use the dist() function and a for() loop to acquire the distances within each community

comspp.means<- as.data.frame(comspp.means)
alphacomspp.means <- comspp.means						
alphacomspp.means[,2]<- as.character(letters[alphacomspp.means[,2]])

for (i in 1:com.num)	{

						spp.holder <- comspp.means[comspp.means[,1]==i,]
						M <- dist(as.matrix(cbind(spp.holder[,3],spp.holder[,4])))
						M <- as.vector(M)
						distances[[i]] <- M
						distances.means[i] <- mean(M) 
						print(M)
						print(spp.holder)

						}
#for(i in 1:com.num)	{
#	spp.holder <- comspp.means[comspp.means[,1]==i,]
#	for (j in 1:length(unique(spp.holder[,2])))	{
#		segments(spp.holder[3,j],spp.holder[3,j+1]
		