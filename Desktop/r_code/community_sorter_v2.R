

#takes empirical coordinate data and generates communities based on a grid size set by the user

com.sorter <- function(specieswithcoords, gridsize=0.5, digits=9, writeCSV=FALSE, writepath=paste(getwd(),"sim.com.matrix.csv",sep="/"))	{

require(gdistance)

colnames(specieswithcoords) <- c("species", "longitude", "latitude")

options(digits=as.numeric(digits))

#round to 3 digits
specieswithcoords[,2] <- format(round(specieswithcoords[,2], 3), nsmall=3)
specieswithcoords[,3] <- format(round(specieswithcoords[,3], 3),nsmall=3)

specieswithcoords[,2]<- as.numeric(specieswithcoords[,2])
specieswithcoords[,3]<- as.numeric(specieswithcoords[,3])

#obtain absolute bounding box information (min/max for long/lat)

bounding <- matrix(0,2,2)
rownames(bounding) <- c("longitude", "latitude")
colnames(bounding) <- c("min", "max")

bounding[1,1] <- min(specieswithcoords$longitude)
bounding[1,2] <- max(specieswithcoords$longitude)
bounding[2,1] <- min(specieswithcoords$latitude)
bounding[2,2] <- max(specieswithcoords$latitude)

print(bounding)


#obtain appropriate grid size information to the nearest whole number

coord.limits <- {}
coord.limits.rd <- {}

coord.limits[1] <- (bounding[1,2]-bounding[1,1])/(as.numeric(gridsize))
coord.limits[2] <- (bounding[2,2]-bounding[2,1])/(as.numeric(gridsize))

print(coord.limits)


coord.limits.rd <- ceiling(coord.limits)
print(coord.limits.rd)


#build reference lists for grid cells

long <- seq(from=bounding[1,1], to=(bounding[1,2]+(as.numeric(gridsize))), by=(as.numeric(gridsize)))	#set grid size
lat <- seq(from=bounding[2,1], to=(bounding[2,2]+(as.numeric(gridsize))), by=(as.numeric(gridsize)))	#set grid size

print(long)
print(lat)


#final community storage object

coms.names <- paste("com", 1:((length(long)-1)*(length(lat)-1)), sep=".")
coms <- as.list(rep(NA, length(coms.names)))
names(coms) <- coms.names


com<-{}
extent.stack <- {}
com.extent <- as.data.frame(matrix(nrow=length(coms.names),ncol=5)); colnames(com.extent)<- c("Community", "xmin", "xmax", "ymin", "ymax")
com.extent[,1]<- coms.names

#build for loop that picks a cell and tests for points within those ranges

for (i in 1:(length(lat)-1))	{

	for (j in 1:(length(long)-1))	{

		cell.bounding <- matrix(0,2,2)
		rownames(cell.bounding) <- c("longitude", "latitude")
		colnames(cell.bounding) <- c("min", "max")

		cell.bounding[1,1] <- long[j]	#long min (x)
		cell.bounding[1,2] <- long[j+1]	#long max
		cell.bounding[2,1] <- lat[i]	#lat min	(y)
		cell.bounding[2,2] <- lat[i+1]	#lat max

		#print(cell.bounding)



		# printout counter that tests the sequence of grids to take
		#print(c(i,j))

		# test all points for fitting the range of points

		for (m in 1: length(specieswithcoords[,2]))	{

			if ((specieswithcoords[m,2] >= cell.bounding[1,1] & specieswithcoords[m,2] <= cell.bounding[1,2]) & (specieswithcoords[m,3]>= cell.bounding[2,1] & specieswithcoords[m,3] <= cell.bounding[2,2]))	{

				com <- c(com, as.character(specieswithcoords[m,1]))
																																								}	else	{
																																										com <- c(com, NA)
																																										#print(com)
																																											}




												}

		#isolates the unique species names in com
		com <- unique(com)
				#print(com)
				#print(cell.bounding)

			#store extents
		extent.stack <- c(extent.stack, extent(cell.bounding[1,1], xmax=cell.bounding[1,2], ymin=cell.bounding[2,1], ymax=cell.bounding[2,2]))
			for (p in 1:length(extent.stack))	{
				com.extent[p,2:5]<- extent.stack[[p]][1:4]
												}


				#removes all NA values in com
				if (is.na(com)&length(com)==1){com <- NA}	else{com <- na.exclude(com)}
		#holder <- paste("com", (j+((i-1)*length(long))), sep=".")


		coms[[(j+((i-1)*(length(long)-1)))]] <- com

		com<- {}

								}
									}
#print(length(coms))
#print(coms)

#removes communities of NA only
for (n in length(coms):1)	{
	holder <- coms[[n]]
	tfholder <- is.na(holder)
	if (tfholder == TRUE)	{coms[[n]]<-NULL;com.extent<-com.extent[-n,]}

							}


#removes communities of 1 species
for (o in length(coms):1)	{
	if(length(coms[[o]])==1)	{coms[[o]]<-NULL;com.extent<-com.extent[-o,]}

							}

#limit communities to unique groupings
#print(length(coms))
#print(length(coms[!duplicated(coms)]))
coms <- coms[!duplicated(coms)]

#set extents for those species
coms.names.unique <- names(coms)
print(coms.names.unique)
com.extent.unique <-as.data.frame(matrix(nrow=length(coms.names.unique),ncol=5))

#com.extent.unique <- com.extent[union(names(coms.names.unique), com.extent[,1]),]
#print(com.extent.unique)
#print(coms.names.unique[1])


#match community IDs in output object and working object
for (q in 1:nrow(com.extent))	{
		for (u in 1:length(coms))	{
			if(as.character(com.extent[q,1])==coms.names.unique[u])	{com.extent.unique[u,1:5]<- com.extent[q,1:5] }
									}
								}

	colnames(com.extent.unique)<- c("Community", "xmin", "xmax", "ymin", "ymax")
#turn results into binary matrix

sim.com.matrix <- matrix(data= NA, nrow= 1, ncol= length(unique(specieswithcoords$species)))
colnames(sim.com.matrix) <- unique(specieswithcoords$species)

for (k in 1:length(coms))	{
coms[[k]]->holder.sing
row.holder <- rep(NA, times=length(colnames(sim.com.matrix)))
for (j in 1:length(holder.sing))	{
for(i in 1:length(colnames(sim.com.matrix)))	{
	if((colnames(sim.com.matrix)[i])==as.character(holder.sing[j]))	{
	row.holder[i] <- 1
																		}
												}

								}
								#print(row.holder)
								sim.com.matrix<- rbind(sim.com.matrix,row.holder)
								}
								sim.com.matrix <- sim.com.matrix[2:(length(coms)+1),]

								sim.com.matrix[is.na(sim.com.matrix)] <- 0

								if(length(coms)==1) {
								sim.com.matrix <- t(as.matrix(sim.com.matrix))
								rownames(sim.com.matrix) <- "1"
																}	else							{
								rownames(sim.com.matrix) <- 1:length(coms)
																									}


	if(writeCSV==TRUE)	{write.csv(sim.com.matrix, file=as.character(writepath))}

	output <- list(sim.com.matrix=sim.com.matrix, extent.stack=extent.stack, coms=coms, com.extent=com.extent.unique)

	#invisible(sim.com.matrix)
	#invisible(extent.stack)
	#invisible(new("CommunitiesStack", sim.com.matrix=sim.com.matrix, extent.stack=extent.stack, coms=coms, com.extent=com.extent.unique))
	invisible(output)

																	}
