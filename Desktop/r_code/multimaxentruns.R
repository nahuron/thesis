# r script to run maxent quickly for 25 species

#must obtain objects from revised_enmtools text file!

#baseline of function
#try to run maxent function through dismo
# bicolandia.maxent <- maxent(x=env, p=bicolandia$presence.points, a=bicolandia$background.points, nbg=10000, args=c('pictures=TRUE', 'plots=TRUE', 'replicatetype=crossvalidate', 'replicates=5', 'applythresholdrule=minimum training presence', 'visible=TRUE', 'outputformat=logistic' ), removeDuplicates=TRUE, path= "/Users/nicholashuron/Desktop/niche_overlap_tests/test2/") 
#  
#  #use the predict function through dismo
#  bicolandia.predict <- predict(object=bicolandia.maxent, x=env, progress='text', args=c('pictures=TRUE', 'responsecurves=TRUE', 'outputformat=logistic', 'projectionlayers=/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/', 'doclamp=TRUE', 'plots=TRUE','applythresholdrule=minimum training presence'))

jack <- c(15,31,38)            
             
#for (a in 1: length(sp.list))	{
  for (a in 19)	{
    sp.list[[a]]$presence.points <- unique(sp.list[[a]]$presence.points)

	if(nrow(sp.list[[a]]$presence.points) > 9)	{
	
	species.maxent <-   maxent(x=env, p= sp.list[[a]]$presence.points, a= sp.list[[a]]$background.points, args=c('pictures=TRUE', 'plots=TRUE', 'replicatetype=crossvalidate', 'replicates=5', 'applythresholdrule=minimum training presence', 'visible=TRUE', 'outputformat=logistic' ), removeDuplicates=TRUE, path= paste0("/Users/nicholashuron/Desktop/niche_overlap_tests/test2/", as.character(sp.files2[a])))
	
	model.select <- species.maxent@results[8,1:5]
	names(model.select) <- colnames(species.maxent@results)[1:5]
	
	#selects the best model from choices based on highest
	species.model <- which.max(model.select)
	
	species.predict <-  predict(object=species.maxent@models[[as.numeric(species.model)]], x=env, progress='text', args=c('pictures=TRUE', 'responsecurves=TRUE', 'outputformat=logistic', 'projectionlayers=/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/', 'doclamp=TRUE', 'plots=TRUE','applythresholdrule=minimum training presence'))
	writeRaster(species.predict, paste0("/Users/nicholashuron/Desktop/niche_overlap_tests/test2/", as.character(sp.files2[a]), "/", as.character(sp.files2[a]), ".asc"), format="ascii", overwrite=TRUE)
	
	species.maxent <- {}
	model.select <- {}
	species.model <- {}
	species.predict <- {}											
												}
	else if(nrow(sp.list[[a]]$presence.points) > 4 & nrow(sp.list[[a]]$presence.points) < 10)	{
	    args_mx <- c('pictures=TRUE', 'plots=TRUE', 'replicatetype=crossvalidate', 'applythresholdrule=minimum training presence', 'visible=TRUE', 'outputformat=logistic', 'removeduplicates=TRUE', as.character(paste0('replicates=',as.character(nrow(sp.list[[a]]$presence.points)))))
			
	    species.maxent <-   maxent(x=env, p= sp.list[[a]]$presence.points, a= sp.list[[a]]$background.points, args=args_mx, removeDuplicates=TRUE, path= paste0("/Users/nicholashuron/Desktop/niche_overlap_tests/test2/", as.character(sp.files2[a])))
	    
	    #species.maxent <-   maxent(x=env, p= sp.list[[a]]$presence.points, a= sp.list[[a]]$background.points, args=c('pictures=TRUE', 'plots=TRUE', 'replicatetype=crossvalidate', 'replicates=1', 'applythresholdrule=minimum training presence', 'visible=TRUE', 'outputformat=logistic' ), removeDuplicates=TRUE, path= paste0("/Users/nicholashuron/Desktop/niche_overlap_tests/test2/", as.character(sp.files2[a])))
			
	    model.select <- species.maxent@results[8,1:5]
	    names(model.select) <- colnames(species.maxent@results)[1:5]
	    
	    #selects the best model from choices based on highest
	    species.model <- which.max(model.select)
	    
	    
	    
	    species.predict <-  predict(object=species.maxent@models[[as.numeric(species.model)]], x=env, progress='text', args=c('pictures=TRUE', 'responsecurves=TRUE', 'outputformat=logistic', 'projectionlayers=/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/', 'doclamp=TRUE', 'plots=TRUE','applythresholdrule=minimum training presence'))
			writeRaster(species.predict, paste0("/Users/nicholashuron/Desktop/niche_overlap_tests/test2/", as.character(sp.files2[a]), "/", as.character(sp.files2[a]), ".asc"), format="ascii", overwite=TRUE)
	
			species.maxent <- {}
			model.select <- {}
			species.model <- {}
			species.predict <- {}																																
																								}

    #rm(list=c("species.maxent", "model.select", "species.model", "species.predict"))											
								}
	
	
