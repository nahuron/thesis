com.simulator <- function(min, max, n, Species, writeCSV=FALSE, writepath=paste(getwd(),"sim.com.matrix.csv",sep="/"))		{	#function to simulate communities based on a range of possible community sizes set by the user (biologically informed)

if(min>max)		{
	stop("Minimum value is greater than maximum")
				}
else			{				

if( min==max)	{
	null.com.sizes <- rep(min, times=n)			#community size selector set for equal min:max comm sizes
				}
else			{
null.com.sizes <- sample(min:max, n, replace=T)	#community size selector set for unequal min:max comm sizes
				}

Species <- Species

sim.com <- as.list(rep(NA, length(1:n)))

	for (i in 1: n)	{
		sim.com[[i]] <- sample(as.factor(Species), null.com.sizes[i], replace=F)
					}
					

	sim.com.matrix <- matrix(data= NA, nrow= 1, ncol= length(Species))									
colnames(sim.com.matrix) <- Species	
							
for (k in 1:length(sim.com))	{
sim.com[[k]]->holder.sing
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
								sim.com.matrix <- sim.com.matrix[2:(length(sim.com)+1),]
									
								sim.com.matrix[is.na(sim.com.matrix)] <- 0
								
								if(length(null.com.sizes)==1) {
								sim.com.matrix <- t(as.matrix(sim.com.matrix))
								rownames(sim.com.matrix) <- "1"
																}
								else							{
								rownames(sim.com.matrix) <- 1:length(null.com.sizes)
																}

#print(null.com.sizes)
#print(sim.com)
										}
	if(writeCSV==TRUE)	{write.csv(sim.com.matrix, file=as.character(writepath))}
	#print(sim.com.matrix)
	invisible(sim.com.matrix)									
										}

#mac
read.csv("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Geographic/locality_data/Brachymeles.unique.locality.winter2016.csv", header=T)->brach_loc

#linux
read.csv("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/Brachymeles.unique.locality.winter2016.csv", header=T)->brach_loc

#brach_loc <- brach_loc[brach_loc[,1]!="libayani",]
brach_loc <- brach_loc[brach_loc[,1]!="species2",]
brach_loc <- brach_loc[brach_loc[,1]!="species3",]
brach_loc <- brach_loc[brach_loc[,1]!="c.f._bonitae",]

#remove for genetics
brach_loc <- brach_loc[brach_loc[,1]!="vermis",]
brach_loc <- brach_loc[brach_loc[,1]!="vindumi",]
brach_loc <- brach_loc[brach_loc[,1]!="wrighti",]
brach_loc <- brach_loc[brach_loc[,1]!="dalawangdaliri",]
brach_loc <- brach_loc[brach_loc[,1]!="suluensis",]
brach_loc$Species <- droplevels(brach_loc$Species)

#luzon only
brach_loc_l <- subset(brach_loc, brach_loc$Species %in% brach_fr_key_l$Species)
brach_loc_l$Species <- droplevels(brach_loc_l$Species)

#mindanao only
brach_loc_m <- subset(brach_loc, brach_loc$Species %in% brach_fr_key_m$Species)
brach_loc_m$Species <- droplevels(brach_loc_m$Species)


										
com.simulator(2,5,999,sort(unique(brach_loc[,1])), writeCSV=F) -> null.coms
#system.time(com.simulator(2,4,1000, species.legend[,1]) -> tester)