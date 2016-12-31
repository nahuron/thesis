#REVISED ENMTOOLS Runs for Brachymeles


#install.packages("devtools")
library(devtools)
#install_github("danlwarren/ENMTools")
#install_github("danlwarren/ENMTools", ref="apply")
library(ENMTools)
library(raster)

#set working directory for environmental layers

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/")

setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/")

#create environmental layer stack

env.files <- list.files(path = "/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/climate_data", pattern = "\\.asc$", full.names = TRUE)
env <- stack(env.files)
crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"

#read in test species (bicolandia)

#{bicolandia <- enmtools.species()
#bicolandia$species.name <- "bicolandia"
#bicolandia$presence.points <- read.csv("input/bicolandia.csv")[,2:3]
#bicolandia$range <- background.raster.buffer(bicolandia$presence.points, 111000, mask = env)
#bicolandia$background.points <- background.points.buffer(points = bicolandia$presence.points,radius = 111000, n = 10000, mask = env[[1]])
#}

#big loop to run the function for multiple species

sp.files <- list.files(path = "/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/coords_ENMS", pattern = ".csv$", full.names = TRUE)
sp.files <- sp.files[grep('background', sp.files, invert=1)]
sp.files <- sort(gsub("micro.","",sp.files))

#brach.list <- c("brevidactylus", "cobos", "dalawangdaliri", "elerae", "isangdaliri", "libayani", "ligtas", "minimus", "pathfinderi", "suluensis", "tiboliorum", "vermis", "vindumi", "wrighti")
brach.list <- c("brevidactylus", "cobos", "dalawangdaliri", "elerae", "isangdaliri", "ligtas", "minimus", "pathfinderi", "suluensis", "tiboliorum", "vermis", "vindumi", "wrighti")


for(c in 1:length(sp.files)){
  for (d in 1:length(brach.list)){
    if(grepl(brach.list[d], sp.files[c])==TRUE) {
      sp.files[c] <- gsub(brach.list[d], as.character(paste0("micro.",brach.list[d])),sp.files[c])
      print(sp.files[c])
    }
  }
}

#sp.files <- sp.files[-4]   	#remove bonitae                                          
sp.files2 <- list.files(path = "/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/coords_ENMS", pattern = ".csv$", full.names = FALSE)
sp.files2 <- sp.files2[grep('background', sp.files2, invert=1)] 
sp.files2 <- gsub("\\.csv$","",sp.files2)
sp.files2 <- sort(gsub("micro.","",sp.files2))
#sp.files2 <- sp.files2[-4] #remove bonitae    

 sp.list <-as.list(rep(NA, length(sp.files2)))
 names(sp.list) <- sp.files2
 
for (a in 1:length(sp.files2))	{
#for (a in 1)	{
	temp <- enmtools.species()
	temp$species.name <- sp.files2[a]
	location <- sp.files[a]
	temp$presence.points <- read.csv(sp.files[a])[,2:3]
	temp$range <- background.raster.buffer(temp$presence.points, 111000, mask = env)
	temp$background.points <- background.points.buffer(points = temp$presence.points,
                                                   radius = 111000, n = cel.rd[a], mask = env[[1]])
    
    sp.list[[a]] <- temp                                             
	
}
 
#save(sp.list, file=paste0("/home/nicholas/Desktop/","sp.list.v5.rda"))
								
load("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/sp.list.v5.rda")
#hold <- as.data.frame(matrix(NA, ncol=length(sp.list), nrow=length(sp.list)))

#ident.list <- as.list(rep(NA, ((length(sp.files2)*(length(sp.files2)+1))/2)))

#loop to run pairwise NET for all species included


for (b in 13: nrow(hold))	{
#for (b in 17)	{

			#for (c in 1:(b-1)){
			for (c in 10){
			print(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])))
			#ident.hold <- identity.test(species.1=sp.list[[b]], species.2=sp.list[[c]], env=env, type="mx", nreps=3, cores=1)
			#ident.list <- as.list(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])))
			#assign(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])), ident.hold)
		  #save(list=as.character(ident.list[1]), file=paste0("/home/nicholas/Desktop/",ident.list[1],".rda"), compress="xz")
		
			#rm(ident.hold)
			#rm(list=as.character(ident.list[1]))
			#gc(verbose=T)
										}
						}; #rm(ident.list)
							
	
#bg.asym.list <- as.list(rep(NA, ((length(sp.files2)*(length(sp.files2)+1))/2)))

r<-1

#for (d in 1: nrow(hold))	{
  for (d in 1)	{
	for (e in 2:2)		{
		if(names(sp.list[d]) == names(sp.list[e]))	{
											print("Skip due to same species")
										}	else	{
												print(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])))
												bg.asym.hold <- background.test(species.1 = sp.list[[d]], species.2 = sp.list[[e]], env = env, type = "mx", nreps = 4, test.type = "asymmetric" )
												bg.asym.list <- as.list(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])))
												assign(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])), bg.asym.hold)
												save(list=as.character(bg.asym.list[1]), file=paste0("/home/nicholas/Desktop/BAST/",bg.asym.list[1],".rda"), compress="xz")
												rm(bg.asym.hold)
												rm(list=as.character(bg.asym.list[1]))
												gc()
													}
								}
							}
	
	
												                                                                                                  
                                                   
#run empirical maxent run!

#enmtools version of maxent command
bicolandia.mx <- enmtools.maxent(bicolandia, env, nbg=10000, path="/Users/nicholashuron/Desktop/niche_overlap_tests/test/", args=c('pictures=TRUE', 'plots=TRUE', 'replicatetype=crossvalidate', 'replicates=5', 'applythresholdrule=minimum training presence', 'visible=TRUE', 'outputformat=logistic' ))                                 
    
 #try to run maxent function through dismo
 bicolandia.maxent <- maxent(x=env, p=bicolandia$presence.points, a=bicolandia$background.points, nbg=10000, args=c('pictures=TRUE', 'plots=TRUE', 'replicatetype=crossvalidate', 'replicates=5', 'applythresholdrule=minimum training presence', 'visible=TRUE', 'outputformat=logistic' ), removeDuplicates=TRUE, path= "/Users/nicholashuron/Desktop/niche_overlap_tests/test2/") 
  
  #use the predict function through dismo
  bicolandia.predict <- predict(object=bicolandia.maxent, x=env, progress='text', args=c('pictures=TRUE', 'responsecurves=TRUE', 'outputformat=logistic', 'projectionlayers=/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/', 'doclamp=TRUE', 'plots=TRUE','applythresholdrule=minimum training presence'))
