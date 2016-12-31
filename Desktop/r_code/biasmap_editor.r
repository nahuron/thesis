#randomized pathway script to validate communities

#1. load in communities/extents for each community

#2. load in corresponding unique localities

#3. load in bias file for trimming

#4. function (community binary matrix, community extents, bias file)
	#loop that ish by community
	#for each community:
		#clip the bias file to extent of community
		#generate a transition file for that community's extent
		#trim to presences within extent of community
		#run a pairwise test of pathways between all points in the community
		#run a counter for failed pairs?	
		
		
comchecker <- function(com, comextents, biasfile, specieswithcoords)	{

		require(gdistance)
		
		{#round occurrences to 3 digits
		specieswithcoords[,2] <- format(round(specieswithcoords[,2], 3), nsmall=3)
		specieswithcoords[,3] <- format(round(specieswithcoords[,3], 3),nsmall=3)

		specieswithcoords[,2]<- as.numeric(specieswithcoords[,2])
		specieswithcoords[,3]<- as.numeric(specieswithcoords[,3])
		
		for (a in 1:nrow(com))	{	#start the main loop
		
		
		#pull the extent data for each community and store as an extent object
		nthcom.extent <- extent(as.numeric(matrix(comextents [a, 2:5], nrow=2)))
		
		
		#clip the biasfile to the aforementioned extent?
		#There is an extent rounding error here that is accounted for later by ignoring the original extentfile 
		nthcom.biasfile <- crop(biasfile, nthcom.extent, snap='out')
		
		
		#generate a transition file from clipped; normalize matrix
		nthcom.transition <- transition(nthcom.biasfile, transitionFunction=mean, 16)
		#nthcom.transition <- normalize(nthcom.transition, method="symm")
		nthcom.transition <- geoCorrection(nthcom.transition, type="c")
		
		#trim presence points to species in the community 
		#consider making the nthcom.localities file a vector to be turned into a matrix later
		
		nthcom.localities <- as.data.frame(matrix(NA, ncol=3))
		colnames(specieswithcoords) <- c("species", "longitude", "latitude")
		colnames(nthcom.localities) <- c("species", "longitude", "latitude")}
		
		#takes row a and looks for the coordinates for species in the community and adds them to a subset object
		for (b in 1:ncol(com))	{
			if(com[a,b]==1)	{
				sp.hold <- colnames(com)[b]
				
				nthcom.localities <- rbind(nthcom.localities, specieswithcoords[as.character(specieswithcoords[,1])==as.character(sp.hold),])
				
					if(is.na(nthcom.localities[1,])==TRUE)	{nthcom.localities <- nthcom.localities[-1,]}	#removes any NA's from new object
							}
				
								}				
		
		#trim presence points to those for nth community's extent
		nthcom.localities.final <- as.data.frame(matrix(NA, ncol=3))	#final object for coordinate data for community
		colnames(nthcom.localities.final) <- c("species", "longitude", "latitude")	
		
		for (c in 1: nrow(nthcom.localities))	{
		
			if (((nthcom.localities[c,2] >= extent(nthcom.biasfile)[1] & nthcom.localities[c,2] <= extent(nthcom.biasfile)[2]) & (nthcom.localities[c,3]>= extent(nthcom.biasfile)[3] & nthcom.localities[c,3] <= extent(nthcom.biasfile)[4])))	{
				
				nthcom.localities.final <-rbind(nthcom.localities.final, nthcom.localities[c,])
				if(is.na(nthcom.localities.final[1,])==TRUE)	{nthcom.localities.final <- nthcom.localities.final[-1,]}		#removes NAs
						
				
																																																		}

												}
					print(nthcom.localities.final)
																																
			
				#Run the pairwise path test for nth community
				#score each location based on the number of non-infinity connections with all other points
				
				#matrix of connected pathways
				#failure to connect is represented by infinity (Inf)
				nthcom.connections <- as.data.frame(matrix(NA, nrow=nrow(nthcom.localities.final), ncol=nrow(nthcom.localities.final)))
				colnames(nthcom.connections) <- rownames(nthcom.localities.final)
				rownames(nthcom.connections) <- paste(rownames(nthcom.localities.final),nthcom.localities.final$species, sep=".")
				
				#start loop based on number of rows in localities
				for (d in 1: nrow(nthcom.localities.final))		{	#starts d loop
					start <- c(nthcom.localities.final[d,2], nthcom.localities.final[d,3])
						for (e in 1: nrow(nthcom.localities.final))	{	#starts e loop
								end <- c(nthcom.localities.final[e,2], nthcom.localities.final[e,3])
								nthcom.connections[d,e] <- costDistance(nthcom.transition, as.numeric(start), as.numeric(end))
								end<-{}
																	}	#ends e loop
					start<-{}	
																}	#ends d loop
				print(nthcom.connections)
																			
							
				#use the connections object to find the greatest grouping of localities					
				nthcom.connections <- as.data.frame(as.matrix(nthcom.connections))
				com.confirmed <- {}
				com.confirmed.sp <- {}
				
				splitter <- rowSums(is.finite(as.matrix(nthcom.connections)))	#counts the number of finite values in each row
				names(splitter) <- colnames(nthcom.connections)
				
				print(splitter)
				
				
				if(any(duplicated(splitter))==FALSE)	{	#scenario 1: a single max value
					
					print(paste0("Community ", a))
					splitter.max <- which.max(splitter)	#finds the value with the most connections
					com.confirmed <- as.matrix(nthcom.connections[splitter.max,])
					print(com.confirmed)
					com.confirmed[!is.finite(com.confirmed)]<-NA
					print(com.confirmed)
					com.confirmed <- na.exclude(t(as.matrix(com.confirmed)))
					com.confirmed <- rownames(com.confirmed)
					print(com.confirmed)
				
				
					for(e in 1:length(com.confirmed))	{
						com.confirmed.sp <- c(com.confirmed.sp,nthcom.localities.final[com.confirmed[e],1])
														}
						com.confirmed.sp <- unique(com.confirmed.sp)						
												
					#create a new version of row a
					com.rev <- rep(NA, length(com[a,]))
					names(com.rev) <- colnames(com)
				
					for (f in 1:length(com.confirmed.sp))	{
					for(g in 1:length(com.rev))	{
						if((names(com.rev)[g])==as.character(com.confirmed.sp[f]))	{
							com.rev[g] <- 1
																					}
																			
												}						
														}
					com.rev[is.na(com.rev)] <- 0
					com[a,] <- com.rev														
														}
				
				
				else if(any(duplicated(splitter))==TRUE)	{	#scenario 2: multiple max values of redundant species
					print(paste0("Community ", a))
					coms.dups <- splitter[splitter==max(splitter)];
					names(coms.dups) <- names(splitter[splitter==max(splitter)])
					coms.list.dups <-as.list(rep(NA, (length(coms.dups))))
					
					for (h in 1:length(coms.dups))	{
						com.confirmed <- as.matrix(nthcom.connections[h,])
						com.confirmed[!is.finite(com.confirmed)]<-NA
						com.confirmed <- na.exclude(t(as.matrix(com.confirmed)))
						com.confirmed <- rownames(com.confirmed)
							
							for(i in 1:length(com.confirmed))	{
						com.confirmed.sp <- c(com.confirmed.sp,nthcom.localities.final[com.confirmed[i],1])
						
																}
						com.confirmed.sp <- unique(com.confirmed.sp)
						coms.list.dups[[h]] <- com.confirmed.sp
						#coms.list.dups[[h]] <- unique(coms.list.dups[[h]])
						com.confirmed.sp <- {}

													}
						coms.list.dups <- unique(coms.list.dups)
						
							#trim any cases of one species
							coms.list.dups<- coms.list.dups[lengths(coms.list.dups)>=2]
							#coms.list.dups <- list(c("boulengeri", "ligtas"), c("ligtas", "mindorensis"), c("vermis", "vindumi"))			#extra line in to check the oddball possible scenarios
							
															
					#create a new version of row a
				com.rev <- rep(NA, length(com[a,]))
				names(com.rev) <- colnames(com)
				com.confirmed.list.sp <- as.vector(coms.list.dups[[1]])
				
				for (f in 1:length(com.confirmed.list.sp))	{
					for(g in 1:length(com.rev))	{
						if((names(com.rev)[g])==as.character(com.confirmed.list.sp[f]))	{
							com.rev[g] <- 1
																					}
																			
												}						
														}
				com.rev[is.na(com.rev)] <- 0
				com[a,] <- com.rev																		
																
					if(length(coms.list.dups)>1)	{			#scenario 3: multiple unique max values
						#create a new version of row a to use for adding to the end of the output matrix
						com.rev.extra <- rep(NA, length(com[a,]))
						names(com.rev.extra) <- colnames(com)
						extra.coms.holder <- {}
						ident.tester <- {}
					
					for (l in 2:length(coms.list.dups))	{
							extra.coms.holder <- coms.list.dups[[l]]
							
						for (j in 1:length(extra.coms.holder))	{
							for(k in 1:length(com.rev.extra))	{
								if((names(com.rev.extra)[k])==as.character(extra.coms.holder[j]))	{com.rev.extra[k] <- 1} #ends the 1 assigner if
																			
														}	#ends the k for					
																} #ends the j for
						com.rev.extra[is.na(com.rev.extra)] <- 0	#fills in the rest of the new object
						
						
						#test to see if it matches any existing communities in the matrix
						for (m in 1:nrow(com))	{
							ident.tester <- c(ident.tester, identical(com[m,], com.rev.extra))
												}	#ends the m for
								if (any(ident.tester)==TRUE)	{
									print(ident.tester)
									print(length(ident.tester))
									com[as.numeric(which(ident.tester)),] <- com.rev.extra
																}	else	{	#ends the if statement
									com <- rbind(com, com.rev.extra)
									rownames(com) <- 1:nrow(com)
																			}	#ends the else statement
																			
							com.rev.extra <- rep(NA, length(com[a,]))
							names(com.rev.extra) <- colnames(com)															
														}	#ends the l for 
													} #ends the scenario 3 if															
					
															}	#ends the elseif
													
													
													
							
														
														}	#ends the for loop using a
				
				com <- unique(com)		
				print(paste("Number of Unique Communities:", nrow(com), sep=" "))								
				return(com)

																}	#ends the function



notes for the PAIC BIAS module of the community_sorter Script
- obtain a raster file of the accepted model of island connectivity
- prepare that raster file for use in R
- obtain community member data with a way to reference geo coordinates for a community
- run a pairwise comparison that obtains distances for a single community and stores the lines as a SpatialLines object (one at a time? or do all of them all at once and see if you get any that fall under the threshold later)
- use cellFromLine to obtain the cell numbers in the raster and getValues to get the values for elevation for those cells (Raster package)
- compare cell values and decide if PAIC crossing occurs
- if it happens, repeat steps in a one-at-a-time manner, this time finding the two sets of points that cross
- now look at each of those points as rally points and look at which other points do not cross, grouping them with one or the other based on that information (use a counter to figure out which of the rally points should be removed?)
