# script to check the number of raster cells in objects for all training areas!

library(raster)

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/MAXENT_2016/tif.worldclim.ersi.clipped.MCPs/")

train.files <- list.files(path = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/MAXENT_2016/tif.worldclim.ersi.clipped.MCPs", pattern = "\\.tif$", full.names = TRUE)
train <- as.list(rep(NA, length(train.files)))
cel <- {}

for (a in 1:length(train.files))	{
	
	train [[a]] <- raster(train.files[a])
	crs(train[[a]]) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
	hold <- train[[a]]
	cel.hold <- length(hold[values(hold)!="NA"])
	cel <- c(cel, cel.hold)
	
	
								}

cel.rd <- round(cel*.95)
