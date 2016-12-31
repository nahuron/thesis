#####Prep script for community sorter

library(gdistance)

#mac
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data")

#linux
setwd("~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data")

	#setClass(Class="CommunitiesStack", representation(sim.com.matrix="matrix", extent.stack="list", coms="list", com.extent="data.frame"))	#use @ notation to get at individual pieces

#mac
read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/Brachymeles.unique.locality.fall2016.csv", header=T)->brachloc

#linux
read.csv("~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/Brachymeles.unique.locality.fall2016.csv", header=T)->brachloc

#####run script for community sorter

#mac
sorted2 <- com.sorter(brachloc, gridsize=0.05, writeCSV=T, writepath=paste("/Users/nicholashuron/Desktop/niche_overlap_tests/test2","sim.com.matrix.csv",sep="/"))

#linux
sorted2 <- com.sorter(brachloc, gridsize=0.05, writeCSV=T, writepath=paste("~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community","sim_com_matrix_0.95.csv",sep="/"))

sorted2 <- com.sorter(brachloc, gridsize=0.05, writeCSV=F)

write.csv(sorted2$sim.com.matrix, file="~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/sim_com_matrix_0.05.csv")

#####Prep script for bias map editor v2.0

library(raster)

setwd("/Users/nicholashuron/Desktop/niche_overlap_tests/")

#linux
setwd("~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/")

{
#create environmental layer stack

env.files <- list.files(path = "~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/climate_data/", pattern = "\\.asc$", full.names = TRUE)
env <- stack(env.files)
crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"

bath <- raster("~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data/climate_data.tif")
crs(bath) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"

#clip the bathymetric data to the extent of interest

bath.clipped <- intersect(bath, env)
crs(bath.clipped) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"

bath.clipped[bath.clipped>=-120 & bath.clipped<=0]
value<-0;for(i in -120:0){freq(bath.clipped, value=i)->x; value=value+x}
value
length(bath.clipped[bath.clipped>=-120 & bath.clipped<=0])

bath.clipped[bath.clipped>=-120 & bath.clipped<=0] <- 0

rasterToPoints(bath.clipped, fun=function(x){0>=x & -120<=x})->bath.clipped2

#recode all cells with values between -120 and 0 with -130
rec <- matrix(c(-9820,0,-130))
rec <- t(rec)
bath.clippedv2 <- reclassify(bath.clipped, rec)
plot(bath.clippedv2)


#isolate the bicol-samar channel and bring into a separate item
bicsam.extent <- extent(123.8959, xmax=124.780, ymin=12.3278,ymax=12.877)
bicsam <- intersect(bath, bicsam.extent)
plot(bicsam)

#isolate the bohol-cebu channel

bohceb.extent <- extent(123.99035, xmax=124.1344, ymin=10.15742, ymax=10.2893)
bohceb <- intersect(bath, bohceb.extent)
plot(bohceb)

#recode those values to -130 m

bicsam[bicsam>=-120 & bicsam<=0] <- -130
bohceb[bohceb>=-120 & bohceb<=0] <- -130

bath.clippedv2 <- merge(bicsam, bath.clipped)
bath.clippedv2 <- merge(bohceb, bath.clippedv2)

rec2 <- matrix(c(-9821,-121.00001,NA, -121.99999,2685,1),ncol=2)
rec2 <- t(rec2)
bath.clippedv3 <- reclassify(bath.clippedv2, rec2)
crs(bath.clippedv3) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
dev.new()
plot(bath.clippedv3, col=c("red","blue"), breaks=c(NA,0,1))


#save file!
writeRaster(bath.clippedv3, filename="PAICS.asc", format="ascii", overwrite=T)

#obtain cell numbers for the desired extent and the main raster
#bicsam.cells <- cellsFromExtent(bath.clipped, extent=bicsam)

#bicsam.focalcells <- cellsFromExtent(bath.clipped, extent=bicsam[bicsam>=-120 & bicsam<=0])


tr1 <- transition(bath.clippedv3, transitionFunction=mean, 16)
dev.new()
image(raster(tr1))

tr2 <- geoCorrection(tr1, type="c", scl=TRUE)

}


#####run script for bias map editor v2.0

connect<-comchecker(sorted2$sim.com.matrix, sorted2$com.extent, bath.clippedv3, brachloc, writeCSV=F)



#####run full sorter and editor script loop

 setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2")
 
 #linux
 setwd("~/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2")

#for (aa in seq(from=1.00, to=0.05, by=-0.05))	{
for(aa in 0.05)	{

	print(paste0("Replicate: ", "revised_com_", as.character(aa), ".csv"))
		
	draftcom <- com.sorter(brachloc, gridsize=aa, writeCSV=T)
	revisedcom <- comchecker(draftcom$sim.com.matrix, draftcom$com.extent, bath.clippedv3, brachloc, writeCSV=F)
	
	#assign(paste0("revised_com_", as.character(aa)), revisedcom)
	
	write.csv(revisedcom, file=paste(getwd(), paste0("revised_com_", as.character(aa), ".csv"), sep="/"))
	
	draftcom <- {}
	revisedcom <- {}
	rm(list=as.character(paste0("revised_com_", as.character(aa))))
	
												}

#####run code to test phylogenetic relationships

library (ape)
library (geiger)
library(picante)

#setwd("/Volumes/NAHURONJUMP/brachymeles.genetic/")		#work in jump drive folders

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/")	#work in dropbox folders


brachtree <- read.tree("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/Brachymeles Chronogram/Newick.tre")

#brachtree <- read.nexus("brachnexus.nex")

#see what taxa are present in the tree
brachtree$tip.label
#visualize the tree
plot.phylo(brachtree, no.margin=T)

#branches to drop due to no need for outgroup
branches.outgroups <- c("'Lygosoma.bowringii.LSUHC6998'", "'Lygosoma.bowringii.LSUHC6970'", "'Lygosoma.quadrupes.LSUHC8002'", "'Lygosoma.quadrupes.LSUHC8403'", "'Lygosoma.LSUHC9321'", "'Tachydromus.sexilineatus.KUFS217'")
#drop outgroups and put resulting phylogeny in a new object
culled.brachtree <- drop.tip(brachtree, branches.outgroups)
#branches to drop due to multiple specimens per species
branches.dropped <- c("'sp.bon3.PalayPalay.ACD1368'", "'sp.bon1.Zambales.RMB4624'", "'miriamae.Thailand.DSM1293'", "'miriamae.Thailand.DSM1363'", "'talinis.Negros.RMB3305'", "'kadwa.Luzon.CDS4057'", "'tungaoi.Masbate.CDS5123'", "'boulengeri.Polillo.CDS1118'", "'taylori.Cebu.CDS3602'", "'mindorensis.Mindoro.ELR0435'", "'boholensis.Bohol.CDS4473'", "'makusog.Catanduanes.CDS2369'", "'schadenbergi.Zambnga.RMB10525'", "'orientalis.Bohol.CDS4692'", "'muntingkamay.Palali.ELR1278'", "'elerae.Kalinga.CDS3481'", "'tiboliorum.SMindanao.ACD5747'", "'tiboliorum.SMindanao.ACD5727'", "'tiboliorum.Mindanao.H1641'", "'samad.Leyte.CDS3390'", "'samad.Leyte.CDS3138'", "'samad.Leyte.CDS3143'", "'samad.Leyte.CDS3400'", "'samad.Leyte.CDS3151'", "'samad.Samar.CDS2552'", "'samad.Samar.CDS2543'", "'samad.Samar.CDS2762'", "'samad.Samar.CDS2706'", "'pathfinderi.Glan.CDS5192'", "'pathfinderi.Glan.CDS5193'", "'pathfinderi.Glan.CDS5199'", "'pathfinderi.Glan.CDS5242'", "'pathfinderi.Glan.CDS5241'", "'pathfinderi.Glan.CDS5201'", "'pathfinderi.Glan.CDS5198'", "'gracilis.Toril.H1175'", "'gracilis.Mindanao.Davao.H1176'", "'gracilis.Baracatan.H1169'", "'gracilis.Cotabato.ACD5262'", "'gracilis.Mindanao.H1428'", "'gracilis.Mindanao.H1332'", "'gracilis.Mindanao.H1429'", "'bicolor.Luzon.RMB11986'", "'bonitae.Luzon.RMB3681'", "'bonitae.Luzon.ACD6099'", "'bonitae.Polillo.CDS1116'", "'bonitae.Polillo.RMB8902'", "'burksi.Marinduque.CDS3595'", "'burksi.Marinduque.CDS3596'", "'burksi.Mindoro.RMB4896'", "'mapalanggaon.Masbate.CDS5101'", "'tridactylus.Sibalom.GVAG287'", "'tridactylus.Negros.CDS3619'", "'tridactylus.Negros.CDS3618'", "'tridactylus.Negros.CDS3617'", "'tridactylus.Alojipan.CDS1647'", "'tridactylus.Pandan.CDS1529'", "'tridactylus.Pandan.CDS1528'", "'isangdaliri.Aurora.RMB12650'", "'samarensis.Samar.CDS3034'", "'minimus.Catanduanes.CDS2371'", "'minimus.Catanduanes.CDS2372'", "'samarensis.Samar.CDS3036'", "'lukbani.Labo.RMB9681'", "'lukbani.Labo.RMB9629'", "'lukbani.Labo.RMB9661'", "'lukbani.Labo.RMB9909'", "'lukbani.Labo.RMB9908'", "'brevidactylus.Bicol.RMB4039'", "'cobos.Catanduanes.CDS5248'", "'cobos.Catanduanes.CDS5254'", "'bicolandia.Luzon.CDS4050'", "'bicolandia.Luzon.CDS4051'", "'bicolandia.Luzon.CDS4089'", "'cebuensis.Cebu.CDS3597'", "'cebuensis.Cebu.CDS3599'", "'libayani.Lapinig.CDS3645'", "'libayani.Lapinig.CDS3677'", "'libayani.Lapinig.CDS3679'", "'libayani.Lapinig.CDS3692'", "'libayani.Lapinig.CDS3700'", "'libayani.Lapinig.CDS3669'", "'paeforum.Leyte.CDS3417'", "'paeforum.Leyte.CDS3419'", "'paeforum.Leyte.CDS3416'", "'sp.bon2.Lubang.CDS3881'", "'sp.bon2.Lubang.CDS3928'", "'sp.bon4.Kalinga.CDS3484'", "'sp.bon4.Gonzaga.CDS3730'", "'sp.bon4.Gonzaga.CDS3731'", "'sp.bon4.NLuzon.RMB11924'", "'sp.bon4.NLuzon.ACD6175'", "'sp.bon4.NLuzon.ACD2002'", "'sp.bon4.CamigNort.RMB7324'", "'sp.bon4.CamigNort.RMB7286'", "'sp.bon3.PalayPalay.ACD1236'", "'bonitae.Polillo.RMB8901'", "'sp.apus1.Kalimantan.RMBR2040'", "'apus.Sabah.SP06915'", "'hilong.MtHilong.CDS5431'", "'hilong.MtHilong.CDS5432'", "'hilong.MtHilong.CDS5429'")
#drop branches and put resulting phylogeny in a new object
culled.brachtree <- drop.tip(culled.brachtree, branches.dropped)

#visualize the new tree
plot.phylo(culled.brachtree, no.margin=T)
culled.brachtree$tip.label

#save the tree to file
#write.tree(culled.brachtree, file="culledbrachtree_final.tre")
#write.nexus(culled.brachtree, file="culledbrachtree_final.nex")

library(picante)

#read in tree file
brach <- read.nexus("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/Brachymeles Chronogram/culledbrachtree_final.nex")

#revise tips to reflect species names
brach$tip.label <- gsub("\\..*","",brach$tip.label)

#test statistics
#creates phylo distance matrix
phydist_brach <- cophenetic(brach)

#calculate whole tree PSV value
allofem <- rep(1, times=length(brach$tip.label))
names(allofem) <- brach$tip.label
psv(samp=allofem, tree=brach)#0.854



#####run code to obtain ENM niche overlaps for empirical values

library(picante)

#stores filepath for each ascii file output of Maxent
asc.tocompare <- list.files(path= "/Users/nicholashuron/Desktop/niche_overlap_tests/test2/ascii", pattern= ".asc$", full.names=T)

#store shortened filenames to rename cols and rows
species.names <- list.files(path= "/Users/nicholashuron/Desktop/niche_overlap_tests/test2/ascii", pattern= ".asc$", full.names=F)
species.names <-  gsub("\\.asc$","",species.names)

#compare all of them
revised.niche.overlaps <- niche.overlap(asc.tocompare)

#rename results cols/rows for ease of reading
colnames(revised.niche.overlaps) <- species.names
rownames(revised.niche.overlaps) <- species.names

write.csv(revised.niche.overlaps, file="/Users/nicholashuron/Desktop/niche_overlap_tests/communities/revised.niche.overlaps.all.csv", col.names=T, row.names=T)



#EVALUATE BY FR
#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)

#remove as yet undescribed species
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="species2",]
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="species3",]
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="species4",]

#genetic removals
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="wrighti",]
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="suluensis",]
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="dalawangdaliri",]
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="vermis",]
brach_fr_key <- brach_fr_key[brach_fr_key[,1]!="vindumi",]


#isolating results to just luzon species
#pare down to Luzon
brach_fr_key_l <- brach_fr_key[brach_fr_key$FR.code=="L",]

#remove species not from PHL
brach_fr_key_l <- na.omit(brach_fr_key_l)
brach_fr_key_l <- droplevels.data.frame(brach_fr_key_l)

#isolating results to just mindanao species
#pare down to Luzon
brach_fr_key_m <- brach_fr_key[brach_fr_key$FR.code=="M",]

#remove species not from PHL
brach_fr_key_m <- na.omit(brach_fr_key_m)
brach_fr_key_m <- droplevels.data.frame(brach_fr_key_m)


#calculate Luzon tree PSV value
allofem_l <- rep(1, times=length(brach_fr_key_l$Species))
names(allofem_l) <- brach_fr_key_l$Species
psv(samp=allofem_l, tree=brach)  #0.851

#calculate Luzon tree PSV value
allofem_m <- rep(1, times=length(brach_fr_key_m$Species))
names(allofem_m) <- brach_fr_key_m$Species
psv(samp=allofem_m, tree=brach)  #0.739
