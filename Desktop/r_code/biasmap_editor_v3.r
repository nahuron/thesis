#Bias Map Editor v2.0

comchecker <- function(com, comextents, biasfile, specieswithcoords, subset=1:nrow(com), writeCSV=FALSE, writepath=paste(getwd(),"rev.com.matrix.csv",sep="/"))	{	#start of function

		#Required Packages
			require(gdistance)
		
		#Format Input Data
			#round occurrences to 3 digits and ensure numeric class
			{specieswithcoords[,2] <- format(round(specieswithcoords[,2], 3), nsmall=3)
			specieswithcoords[,3] <- format(round(specieswithcoords[,3], 3),nsmall=3)

			specieswithcoords[,2]<- as.numeric(specieswithcoords[,2])
			specieswithcoords[,3]<- as.numeric(specieswithcoords[,3])}
		
		#initialize final output object
			#sp.by.com.final <- as.list(rep(NA, nrow(com)))
			com.final <- com
		
		#Pull out each community individually
		for (a in subset)	{	#starts a loop
			
			#pull the community into a new object
			nthcom <- com[a,]
			names(nthcom) <- colnames(com)
		
			#pull the extent data for each community and store as an extent object
			nthcom.extent <- extent(as.numeric(matrix(comextents [a, 2:5], nrow=2)))
		
		
			#clip the biasfile to the aforementioned extent
			#There is an extent rounding error here that is accounted for later by ignoring the original extentfile 
			nthcom.biasfile <- crop(biasfile, nthcom.extent, snap='out')
		
		
			#generate a transition file from clipped; normalize matrix
			nthcom.transition <- transition(nthcom.biasfile, transitionFunction=mean, 16)
			nthcom.transition <- geoCorrection(nthcom.transition, type="c")
		
			#create a community-specific coordinate object
			nthcom.localities <- as.data.frame(matrix(NA, ncol=3))
			colnames(specieswithcoords) <- c("species", "longitude", "latitude")
			colnames(nthcom.localities) <- c("species", "longitude", "latitude")
			
			#trim presence points to species in the community
			for (b in 1:length(nthcom))	{	#starts b loop
				if(nthcom[b]==1)	{	#starts if in b
					sp.hold <- colnames(com)[b]
				
					nthcom.localities <- rbind(nthcom.localities, specieswithcoords[as.character(specieswithcoords[,1])==as.character(sp.hold),])
				
						if(is.na(nthcom.localities[1,])==TRUE)	{nthcom.localities <- nthcom.localities[-1,]}	#removes any NA's from new object
									}	#ends if in b
				
									}	#ends b loop		
										
			#Create final locality object for coordinate data for community by extent
			nthcom.localities.final <- as.data.frame(matrix(NA, ncol=3))
			colnames(nthcom.localities.final) <- c("species", "longitude", "latitude")
			
			for (c in 1: nrow(nthcom.localities))	{	#starts c loop
		
			if (((nthcom.localities[c,2] >= extent(nthcom.biasfile)[1] & nthcom.localities[c,2] <= extent(nthcom.biasfile)[2]) & (nthcom.localities[c,3]>= extent(nthcom.biasfile)[3] & nthcom.localities[c,3] <= extent(nthcom.biasfile)[4])))	{
				
				nthcom.localities.final <-rbind(nthcom.localities.final, nthcom.localities[c,])
				if(is.na(nthcom.localities.final[1,])==TRUE)	{nthcom.localities.final <- nthcom.localities.final[-1,]}		#removes NAs
						
				
																																																												}
													}	#ends c loop
			
			#Create Pathway Matrix
				#develop results matrix
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
									end <- {}
																	}	#ends e loop
						start <- {}	
																}	#ends d loop	

			#Evaluate Resultant Pathway Matrix										
				#initialize storage object for groupings to compare
					sp.by.com <- as.list(rep(NA, nrow(nthcom.connections))) 
					names(sp.by.com) <- rownames(nthcom.connections)
				#obtain matrix as true==finite distance, false==infinite value (not connected to row species)	
					nthcom.finite <- is.finite(as.matrix(nthcom.connections))
					colnames(nthcom.finite) <- colnames(nthcom.connections)
					rownames(nthcom.finite) <- rownames(nthcom.connections)
				#evaluate each row for connections
					for (f in 1:nrow(nthcom.connections))	{	#starts f loop

						as.numeric(colnames(nthcom.finite[,nthcom.finite[f,]==TRUE])) -> row.nthcom.finite	#list of all columns in row
						row.sps <- {}
						
						for (g in 1:length(row.nthcom.finite))	{	#starts g loop
	
							holder <- row.nthcom.finite[g]	#individual column in seqence for the row of interest
							sp.holder <- specieswithcoords[rownames(specieswithcoords)==holder,1]	#hold the species
							row.sps <- c(row.sps, as.character(sp.holder))	#add the species to the object
							row.sps <- (unique(row.sps))	#reduce to unique species only
	
							
																}	#ends g loop
						sp.by.com[[f]] <- row.sps						
							
															}	#ends f loop
				#reduce possible subcommunities down to only unique species groupings
					sp.by.com <- unique(sp.by.com)
				
				#pathway for evaluation of community/subcommunity grouping
					#remove elements that contain a single species
						for (h in 1: length(sp.by.com))	{	#starts h loop
							if(length(sp.by.com[[h]])<2)	{	#starts if in h
								sp.by.com[[h]] <- NA								#replaces single species entries with NA
															}	#ends if in h
															}	#ends h loop
															
								sp.by.com <- sp.by.com[!is.na(sp.by.com)]			#removes list elements of NA
								print(sp.by.com)
								#print(nthcom.connections)				
											
					#scenario 1: sp.by.com is now empty (only entries of 1 sp)
						if(length(sp.by.com) == 0)	{ #starts if in scenario 1
							print("Scenario 1")
							print(paste0("Community ", a, " has 0 connected species"))
							com.final[a,] <- rep(NA, times=as.numeric(length(com[a,])))
							print(nrow(com.final))
													}	#ends if in scenario 1
						
					#scenario 2: sp.by.com now contains a single entry (one max sp richness	entry)
						else if(length(sp.by.com) == 1)	{	#starts if in scenario 2
							print("Scenario 2")
							print(paste0("Community ", a, " has ", length(sp.by.com[[1]]), " connected species"))
							nthcom.finalrow <- rep(NA, as.numeric(length(com[a,])))	#initialize the row holder object to be inserted later
							names(nthcom.finalrow) <- names(com[a,])
														
								for (i in 1:length(sp.by.com[[1]]))	{	#starts i loop
									for(j in 1:length(nthcom.finalrow))	{	#starts j loop
										if((names(nthcom.finalrow)[j]) == as.character(sp.by.com[[1]][i]))	{nthcom.finalrow[j] <- 1}

																	}	#ends j loop					
																		}	#ends i loop
							nthcom.finalrow[is.na(nthcom.finalrow)] <- 0
							com.final[a,] <- nthcom.finalrow
							print(nrow(com.final))																														
							
							#objects to be emptied
								nthcom.finalrow <- {}
								sp.by.com <- {}							
														}	#ends if in scenario 2
						
					#scenario 3: sp.by.com now has more than one entry	(multiple entries, more than one max sp richness entry)
						else if (length(sp.by.com) > 1 & length(which(lengths(sp.by.com) == max(lengths(sp.by.com)))) > 1)	{	#starts if in scenario 3
							print("Scenario 3")
							print(paste0("This Community can be split into ", length(which(lengths(sp.by.com) == max(lengths(sp.by.com)))), " smaller communities with the same max number of species of ", max(lengths(sp.by.com))))
							nthcom.finalrow <- rep(NA, as.numeric(length(com[a,])))	#initialize the row holder object to be inserted later
							names(nthcom.finalrow) <- names(com[a,])
							multimax <- which(lengths(sp.by.com)==max(lengths(sp.by.com)))
							
							nthcom.maxcom <- sp.by.com[[multimax[1]]]
							nthcom.maxcom <- unlist(nthcom.maxcom)
							print(sp.by.com)
							print(multimax)
							
							for (i in 1:length(nthcom.maxcom))	{	#starts i loop
									for(j in 1:length(nthcom.finalrow))	{	#starts j loop
										if((names(nthcom.finalrow)[j]) == as.character(nthcom.maxcom[i]))	{nthcom.finalrow[j] <- 1}

																	}	#ends j loop					
																		}	#ends i loop
							nthcom.finalrow[is.na(nthcom.finalrow)] <- 0
							com.final[a,] <- nthcom.finalrow
							print(nrow(com.final))
							
							for(k in multimax[-1])	{	#starts k loop
								nthcom.maxcom <- sp.by.com[[k]]
								nthcom.maxcom <- unlist(nthcom.maxcom)
								nthcom.finalrow <- rep(NA, as.numeric(length(com[a,])))	#initialize the row holder object to be inserted later
								names(nthcom.finalrow) <- names(com[a,])
								
								for (i in 1:length(nthcom.maxcom))	{	#starts i loop
									for(j in 1:length(nthcom.finalrow))	{	#starts j loop
										if((names(nthcom.finalrow)[j]) == as.character(nthcom.maxcom[i]))	{nthcom.finalrow[j] <- 1}

																	}	#ends j loop					
																		}	#ends i loop
								nthcom.finalrow[is.na(nthcom.finalrow)] <- 0
								com.final <- rbind(com.final, nthcom.finalrow)
								rownames(com.final) <- 1:nrow(com.final)
								print(nrow(com.final))
								
								#objects to be emptied
								nthcom.maxcom <- {}
								nthcom.finalrow <- {}
													}	#ends k loop
													
								#remove the max from the sp.by.com object since we are done with it					
								sp.by.com[which(lengths(sp.by.com)==max(lengths(sp.by.com)))] <- NA
								print(sp.by.com)
								sp.by.com <- sp.by.com[!is.na(sp.by.com)]			#removes list elements of NA
								print(sp.by.com)	
									
								if(length(sp.by.com) > 0)	{	
									#add any remaining subcommunities to the end of the final output object
									for (l in 1:length(sp.by.com))	{	# starts l loop
									print(paste0("Extra Community ", l, " from Original Community ", a, " has ", length(sp.by.com[[l]]), " connected species"))
									nthcom.finalrow <- rep(NA, as.numeric(length(com[a,])))	#initialize the row holder object to be inserted later
									names(nthcom.finalrow) <- names(com[a,])
									
									for (i in 1:length(sp.by.com[[l]]))	{	#starts i loop
									for(j in 1:length(nthcom.finalrow))	{	#starts j loop
										if((names(nthcom.finalrow)[j]) == as.character(sp.by.com[[l]][i]))	{nthcom.finalrow[j] <- 1}

																		}	#ends j loop					
																		}	#ends i loop
									nthcom.finalrow[is.na(nthcom.finalrow)] <- 0
									com.final <- rbind(com.final, nthcom.finalrow)
									rownames(com.final) <- 1:nrow(com.final)
									print(nrow(com.final))
									
									#objects to be emptied
									nthcom.finalrow <- {}
									
																	}	#ends l loop
															}		
								#objects to be emptied
								nthcom.finalrow <- {}
								sp.by.com <- {}									
							
																															}	#ends if in scenario 3	
																															
					#scenario 4: sp.by.com now has more than one entry	(multiple entries, one max sp richness entry)
						else if (length(sp.by.com) > 1 & length(which(lengths(sp.by.com) == max(lengths(sp.by.com)))) == 1)	{	#starts if in scenario 4
							print("Scenario 4")
							print(paste0("Community ", a, " has ", max(lengths(sp.by.com)), " connected species"))
							nthcom.finalrow <- rep(NA, as.numeric(length(com[a,])))	#initialize the row holder object to be inserted later
							names(nthcom.finalrow) <- names(com[a,])
							
							nthcom.maxcom <- sp.by.com[[which(lengths(sp.by.com) == max(lengths(sp.by.com)))]]	#object that holds species in max richness subset
							nthcom.maxcom <- unlist(nthcom.maxcom)
							
							for (i in 1:length(nthcom.maxcom))	{	#starts i loop
									for(j in 1:length(nthcom.finalrow))	{	#starts j loop
										if((names(nthcom.finalrow)[j]) == as.character(nthcom.maxcom[i]))	{nthcom.finalrow[j] <- 1}

																	}	#ends j loop					
																		}	#ends i loop
							nthcom.finalrow[is.na(nthcom.finalrow)] <- 0
							com.final[a,] <- nthcom.finalrow
							print(nrow(com.final))
							
							#remove the max from the sp.by.com object since we are done with it
								#print(sp.by.com[[which(lengths(sp.by.com)==max(lengths(sp.by.com)))]])
								sp.by.com[[which(lengths(sp.by.com)==max(lengths(sp.by.com)))]] <- NA
								sp.by.com <- sp.by.com[!is.na(sp.by.com)]			#removes list elements of NA
								#print(sp.by.com)
							
							#add any remaining subcommunities to the end of the final output object
								for (m in 1:length(sp.by.com))	{	# starts m loop
									print(paste0("Extra Community ", m, " from ", "Original Community ", a, " has ", length(sp.by.com[[m]]), " connected species"))
									nthcom.finalrow <- rep(NA, as.numeric(length(com[a,])))	#initialize the row holder object to be inserted later
									names(nthcom.finalrow) <- names(com[a,])
									
									for (i in 1:length(sp.by.com[[m]]))	{	#starts i loop
									for(j in 1:length(nthcom.finalrow))	{	#starts j loop
										if((names(nthcom.finalrow)[j]) == as.character(sp.by.com[[m]][i]))	{nthcom.finalrow[j] <- 1}

																		}	#ends j loop					
																		}	#ends i loop
									nthcom.finalrow[is.na(nthcom.finalrow)] <- 0
									com.final <- rbind(com.final, nthcom.finalrow)
									rownames(com.final) <- 1:nrow(com.final)
									print(nrow(com.final))
									
																}	#ends m loop
							
							#objects to be emptied
								nthcom.finalrow <- {}
								sp.by.com <- {}	
																															}	#ends if in scenario 4						
									
				
																		}	#ends a loop
			#Evaluate Resultant Com Matrix and Remove Duplicates
				
			#Return Final Output
				
				com.final <- com.final[!is.na(com.final)]			#removes rows of elements of NA
				com.final <- matrix(com.final, ncol=ncol(com))
				colnames(com.final) <- colnames(com)
				rownames(com.final) <- 1:nrow(com.final)
				
				if(writeCSV==TRUE)	{write.csv(com.final, file=as.character(writepath))}	
				
				return(com.final)
																		
																			}	#ends function
			
