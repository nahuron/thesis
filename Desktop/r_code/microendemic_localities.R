#microendemic species extra presence points finder

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/microendemics/")

#pull in raster stack of interest

env.files <- list.files(path = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/climate_data", pattern = "\\.asc$", full.names = TRUE)
env <- stack(env.files)
crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"

#pull in coordinates of interest
brach_loc <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/Brachymeles.unique.locality.fall2016.csv", header=T)

{brach_loc <- brach_loc[brach_loc[,1]=="isangdaliri",]

#identify the cell that corresponds to the point(s) of interest

cellFromXY(env, brach_loc[,2:3])

#identify adjacent pixels

adj <- adjacent(env, cells=432427, directions=16, id=TRUE)

#isolate cell numbers for adjacent pixels
adj.cells <- adj[,3]

#obtain XY data for adjacent pixels
new.coords <- xyFromCell(env, adj.cells)

#build an output style version of the XY data based on input
new.coords <- cbind(rep("isangdaliri", times=nrow(new.coords)), as.data.frame(new.coords))
colnames(new.coords) <- colnames(brach_loc)

#final output
final.coords <- rbind(new.coords, brach_loc)}



adjcoords <- function(rasterlayer, coords, adjacent.type=8, snumber=8, writeCSV=FALSE, writepath=paste(getwd(),"sampled.coordinates.csv",sep="/"))						{

if (snumber>adjacent.type*nrow(coords))	{stop("Sampled Greater Than Max Adjacent Points")} else	{
require(raster)
	#identify the cell(s) that correspond(s) to the point(s) of interest
		focal.cells <- cellFromXY(rasterlayer, coords)
	
	#identify adjacent pixels
		adjacent.cells <- adjacent(rasterlayer, cells=focal.cells, directions=adjacent.type, id=TRUE)
	#isolate cell numbers for adjacent pixels
		adjacent.cells <- adjacent.cells[,3]
	
	#if n is less than number of cells, sample for n
		if(length(adjacent.cells)>snumber)	{adjacent.cells <- sample(adjacent.cells, size=snumber, replace=FALSE)}	
		
	#obtain XY data for adjacent pixels
		new.coords <- xyFromCell(rasterlayer, adjacent.cells)
		new.coords <- unique(new.coords)
		colnames(new.coords) <- colnames(coords)
	#add the original coordinates at the end
		new.coords <- rbind(new.coords, coords)
		print(new.coords)
		
if(writeCSV==TRUE)	{write.csv(new.coords,row.names=FALSE, file=as.character(writepath))}		
		
return(new.coords)									
																					}}
																				
																					

brach.list <- c("brevidactylus", "cobos", "dalawangdaliri", "elerae", "isangdaliri", "ligtas", "minimus", "pathfinderi", "suluensis", "tiboliorum", "vermis", "vindumi", "wrighti")
brach.newcoords.list <- as.list(rep(NA, (length(brach.list))))
names(brach.newcoords.list) <- brach.list

#for (a in 1:length(brach.list))	{
  for (a in 5)	{ 
	brach.holder <- brach_loc[brach_loc[,1]==brach.list[a],2:3]
	brach.newcoords.list[[a]] <- adjcoords(env, brach.holder, adjacent.type=8, snumber=(8*nrow(brach.holder)), writeCSV=FALSE)
	brach.newcoords.list[[a]] <- cbind(rep(brach.list[a], times=nrow(brach.newcoords.list[[a]])), brach.newcoords.list[[a]])
	colnames(brach.newcoords.list[[a]]) <- c("species", "longitude", "latitude")
	
	#print(brach.newcoords.list[[a]])
	
	write.csv(brach.newcoords.list[[a]], paste(getwd(), paste("micro",brach.list[a],"csv", sep="."), sep="/"), row.names=FALSE)
	
	#return(brach.newcoords.list)
								}
	
#notes on Microendemic Pseudoreplication Presence Method
#- pull in raster stack of interest
#- pull in coordinates of interest
#- identify the pixels that correspond to the coordinates in the raster layer(s)
#- identify all adjacent pixels and determine their coordinates
#- obtain coordinates in output format
#- optional additional step: rather than take all adjacent pixels, sample a needed number of those pixels, specified "x"
#- use the raster package function called xyFromCell
																					