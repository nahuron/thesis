#REVISED ENMTOOLS Runs for Brachymeles

#install.packages("devtools")
library(devtools)
#install_github("danlwarren/ENMTools")
#install_github("danlwarren/ENMTools", ref="apply")
library(ENMTools)
library(raster)

#set working directory for environmental layers
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/")

#create environmental layer stack

env.files <- list.files(path = "/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/climate_data", pattern = "\\.asc$", full.names = TRUE)
env <- stack(env.files)
crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"


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

sp.files2 <- list.files(path = "/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/coords_ENMS", pattern = ".csv$", full.names = FALSE)
sp.files2 <- sp.files2[grep('background', sp.files2, invert=1)] 
sp.files2 <- gsub("\\.csv$","",sp.files2)
sp.files2 <- sort(gsub("micro.","",sp.files2))
#sp.files2 <- sp.files2[-4] #remove bonitae    

 sp.list <-as.list(rep(NA, length(sp.files2)))
 names(sp.list) <- sp.files2

#build species data  
for (a in 1:length(sp.files2))	{
#for (a in 10)	{
	temp <- enmtools.species()
	temp$species.name <- sp.files2[a]
	location <- sp.files[a]
	temp$presence.points <- read.csv(sp.files[a])[,2:3]
	temp$range <- background.raster.buffer(temp$presence.points, 111000, mask = env)
	temp$background.points <- background.points.buffer(points = temp$presence.points,
                                                   radius = 111000, n = cel.rd[a], mask = env[[1]])
    
    sp.list[[a]] <- temp                                             
	
}
 
#save(sp.list, file=paste0("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/","sp.list.v5.rda"))
load("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/sp.list.v5.rda")

#check data to ensure that there aren't too many NA's for bg points
for (a in 1: length(sp.list)){
  extract(env, sp.list[[a]]$presence.points, na.rm=FALSE) -> test
  test <- test[is.na(test)]
  print(sp.list[[a]]$species.name)
  print(head(test))
  extract(env, sp.list[[a]]$background.points, na.rm=FALSE) -> test
  test <- test[is.na(test)]
  print(head(test))
}

#plot points of both types for all species against elevation
for (a in 1:length(sp.list)){
  plot(env$alt, col=rev(viridis(10)), xlim=c(min(sp.list[[a]]$background.points$longitude),max(sp.list[[a]]$background.points$longitude)), ylim=c(min(sp.list[[a]]$background.points$latitude),max(sp.list[[a]]$background.points$latitude)), main=sp.list[[a]]$species.name)
  points(sp.list[[a]]$background.points$longitude,sp.list[[a]]$background.points$latitude, pch=22, col="black", bg="red")
  points(sp.list[[a]]$presence.points$longitude,sp.list[[a]]$presence.points$latitude, pch=22, col="black", bg="blue")
}

#ident.list <- as.list(rep(NA, ((length(sp.files2)*(length(sp.files2)+1))/2)))

#loop to run pairwise NET for all species included

brach.list.index <- c(4,5,6, 9, 10,11,15,17,18,23,27,31,34,37,38,40)

brach.notlist.index <- seq(1,40,1)
brach.notlist.index <- setdiff(brach.notlist.index, brach.list.index)

brach.bothlist.index <- seq(1,40,1)

#set external as writing source
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/NIT/")

for (b in sp.files2[32])	{
  
  for (c in sp.files2[6]){
    if(c>b)	{print("Skip due to later comparison")}
    else if(names(sp.list[b]) == names(sp.list[c]))	{print("Skip due to same species")}
    else	{
      print(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])))
      ident.hold <- identity.test(species.1=sp.list[[b]], species.2=sp.list[[c]], env=env, type="mx", nreps=99, cores=1)
      ident.list <- as.list(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])))
      assign(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])), ident.hold)
      save(list=as.character(ident.list[1]), file=paste0("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/NIT/",ident.list[1],".rda"), compress="xz")
      
      rm(ident.hold)
      rm(list=as.character(ident.list[1]))
      gc(verbose=T)
    }
  }
}
rm(ident.list)							
	
#bg.asym.list <- as.list(rep(NA, ((length(sp.files2)*(length(sp.files2)+1))/2)))

#forward loop
for (d in brach.notlist.index[3:24])	{
  for (e in brach.list.index)		{
    if(names(sp.list[d]) == names(sp.list[e]))	{
      print("Skip due to same species")
    }	else	{
      print(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])))
      bg.asym.hold <- background.test(species.1 = sp.list[[d]], species.2 = sp.list[[e]], env = env, type = "mx", nreps = 99, test.type = "asymmetric" )
      bg.asym.list <- as.list(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])))
      assign(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])), bg.asym.hold)
      save(list=as.character(bg.asym.list[1]), file=paste0(getwd(),"/",bg.asym.list[1],".rda"), compress="xz")
      rm(bg.asym.hold)
      rm(list=as.character(bg.asym.list[1]))
      gc()
    }
  }
}

#reverse loop	
for (d in 17:1)	{
  for (e in ncol(hold):1)		{
    if(names(sp.list[d]) == names(sp.list[e]))	{
      print("Skip due to same species")
    }	else	{
      print(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])))
      bg.asym.hold <- background.test(species.1 = sp.list[[d]], species.2 = sp.list[[e]], env = env, type = "mx", nreps = 99, test.type = "asymmetric" )
      bg.asym.list <- as.list(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])))
      assign(paste0("BAST_",names(sp.list[d]), "_vs_bg_", names(sp.list[e])), bg.asym.hold)
      save(list=as.character(bg.asym.list[1]), file=paste0(getwd(),"/BAST/",bg.asym.list[1],".rda"), compress="xz")
      rm(bg.asym.hold)
      rm(list=as.character(bg.asym.list[1]))
      gc()
    }
  }
}	

#Identity test
for (b in brach.notlist.index[19:24])	{
  
  for (c in brach.list.index){
    if(c>b)	{print("Skip due to later comparison")}
    else if(names(sp.list[b]) == names(sp.list[c]))	{print("Skip due to same species")}
    else	{
      print(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])))
      ident.hold <- identity.test(species.1=sp.list[[b]], species.2=sp.list[[c]], env=env, type="mx", nreps=99, cores=1)
      ident.list <- as.list(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])))
      assign(paste0("NET_",names(sp.list[b]), "_vs_", names(sp.list[c])), ident.hold)
      save(list=as.character(ident.list[1]), file=paste0(getwd(),"/",ident.list[1],".rda"), compress="xz")
      
      rm(ident.hold)
      rm(list=as.character(ident.list[1]))
      gc(verbose=T)
    }
  }
}
rm(ident.list)	

