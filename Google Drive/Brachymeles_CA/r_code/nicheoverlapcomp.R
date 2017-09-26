#script for doing the background similarity test

library(phyloclim)

#must obtain datapoints as a dataframe

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/ENMS")
#compare ENMs for niche conservation

brach.niche <- c("bicolandia_0_Philippines.asc", "bonitae_Philippines_avg.asc", "boulengeri_Philippines_avg.asc", "burksi_Philippines_avg.asc", "cebuensis_Philippines.asc", "cobos_Philippines.asc", "gracilis_Philippines_avg.asc", "hilong_Philippines_avg.asc", "kadwa_Philippines_avg.asc", "libayani_Philippines.asc", "luzoni_Philippines_avg.asc", "mapalanggaon_Philippines.asc", "mindorensis_Philippines_avg.asc", "muntingkamay_Philippines.asc", "orientalis_Philippines_avg.asc", "paeforum_Philippines.asc", "samad_Philippines_avg.asc", "samarensis_Philippines.asc", "schadenbergi_Philippines_avg.asc", "talinis_Philippines_avg.asc", "taylori_Philippines_avg.asc", "tridactylus_Philippines_avg.asc", "tungaoi_Philippines.asc", "vulcani_Philippines_avg.asc", "boholensis_Philippines_avg.asc", "bicolor_Philippines.asc")

#test niche overlap function
niche.overlap(brach.niche) -> brach.niche.overlap

write.csv(brach.niche.overlap, "niche.overlap.matrix.csv")

#isolate the two stats
brach.niche.overlap[lower.tri(brach.niche.overlap)] -> brach.i
brach.niche.overlap[upper.tri(brach.niche.overlap)] -> brach.d

#isolate two communities of interest: com 31

#com 31

com31 <- c("boholensis", "orientalis", "samad")

com31.i <- c(brach.niche.overlap["boholensis_Philippines_avg.asc", "orientalis_Philippines_avg.asc"], brach.niche.overlap["boholensis_Philippines_avg.asc", "samad_Philippines_avg.asc"], brach.niche.overlap["orientalis_Philippines_avg.asc","samad_Philippines_avg.asc"])

# path to MAXENT
# --------------
maxent.exe <- "/Users/nicholashuron/Desktop/niche_overlap_tests/maxent/maxent.jar"
# a data frame of coordinates where two species
# have been detected ('presence points') and
# a raster stack of environmental covariables
# --------------------------------------

library(raster)
library(dismo)
library(maps)
library(rgdal)
library(maptools)

data(wrld_simpl)
plot(wrld_simpl, xlim=c(118,120), ylim=c(0,20))
points(sites.tal, pch=17, col="blue")
points(sites.tri, pch=17, col="red")

species <- c("talinis", "tridactylus")
sites <-read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/locality_data/Brachymeles.unique.locality.May2016.csv", header=TRUE, sep=",")

sites[sites[,1]=="tridactylus",]->sites.tri
sites.tri<- sites.tri[,-1]
sites.tri <- unique(sites.tri)

fold.tri <- kfold(sites.tri, k=5)
test.tri <- sites.tri[fold.tri==1,]
train.tri <- sites.tri[fold.tri!=1,]

sites[sites[,1]=="talinis",]->sites.tal
sites.tal<- sites.tal[,-1]
sites.tal <- unique(sites.tal)

#ascii files of bioclim variables to project results on
projected.bioclimvars <- list.files(path= "/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/", pattern= ".asc$", full.names=T)

stack(projected.bioclimvars) -> predictors

read.asciigrid(projected.bioclimvars[1]) ->bio1.pro
read.asciigrid(projected.bioclimvars[2]) ->bio12.pro
read.asciigrid(projected.bioclimvars[3]) ->bio14.pro
read.asciigrid(projected.bioclimvars[4]) ->bio2.pro
read.asciigrid(projected.bioclimvars[5]) ->bio3.pro

#ascii files of bioclim variables to train from
tri.bioclimvars <- list.files(path= "/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/training/tridactylus/", pattern= ".asc$", full.names=T)
tal.bioclimvars <- list.files(path= "/Users/nicholashuron/Desktop/niche_overlap_tests/climate_data/training/talinis/", pattern= ".asc$", full.names=T)

#tridactylus
read.asciigrid(tri.bioclimvars[1]) ->bio1.tri
read.asciigrid(tri.bioclimvars[2]) ->bio12.tri
read.asciigrid(tri.bioclimvars[3]) ->bio14.tri
read.asciigrid(tri.bioclimvars[4]) ->bio2.tri
read.asciigrid(tri.bioclimvars[5]) ->bio3.tri

stack(tri.bioclimvars) -> predictors.tri

#talinis
read.asciigrid(tal.bioclimvars[1]) ->bio1.tal
read.asciigrid(tal.bioclimvars[2]) ->bio12.tal
read.asciigrid(tal.bioclimvars[3]) ->bio14.tal
read.asciigrid(tal.bioclimvars[4]) ->bio2.tal
read.asciigrid(tal.bioclimvars[5]) ->bio3.tal


#directory to store results

outdir <- "/Users/nicholashuron/Desktop/niche_overlap_tests/outputtester"


setwd(outdir)
mx.tri <- maxent(x=predictors, p=train.tri, path=outdir, args=c("outputformat=raw", "replicates=5"), removeDuplicates=T)
, "projectionlayers=c(bio1.pro, bio2.pro, bio3.pro, bio12.pro, bio14.pro)"

bgst <-bg.similarity.test(p=sites, env=predictors, n=99, conf.level=0.95, app=maxent.exe, dir=outdir)


preds <- list.files(path = data.path, pattern = "[.]asc")





net <- niche.equivalency.test(samples, preds, reps, maxent.exe, dir)