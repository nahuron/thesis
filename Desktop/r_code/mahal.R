#-----------------------------------------------------------------------------------------
library(FactoMineR)
library(pcaMethods)
library(corrplot)
library(Hmisc)
library(maptools)
library(HDMD)

#create a new class
setClass(Class="com_distances", representation(sim.com.distances="list", emp.com.distances="list", sim.com.mdistances="vector", emp.com.mdistances="vector"))	#use @ notation to get at individual pieces

#set working space mac
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology")


#Old Version of reading in data
#-----------------------------------------------------------------------------------------
#load in dataset 1 (here it is the morphometric data provided via CDS per individual with species identifier): 17
read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_2016_brachymeles_morph_v6_rawmeans.txt", header=T, fill=T, nrows=600) -> brach_morph

#load in no TL version of dataset: 16
read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_2016_brachymeles_morph_v6_rawmeans_noTL.txt", header=T, fill=T, nrows=600) -> brach_morph

#load in absolutely reduced version of dataset: 10
read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_2016_brachymeles_morph_v5_rawmeans_reduced.txt", header=T, fill=T, nrows=600) -> brach_morph

read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_2016_brachymeles_morph_v5_rawmeans_reduced_noTL.txt", header=T, fill=T, nrows=600) -> brach_morph

#reduced noTL
read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_2016_brachymeles_morph_v5_rawmeans_reduced.txt", header=T, fill=T, nrows=600) -> brach_morph
brach_morph <- brach_morph[!colnames(brach_morph) %in% "TL.HL"]
#-----------------------------------------------------------------------------------------

#New Method for reading in data
#read in morphology object
read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_brachymeles_morph_raw_means.csv", header=T, fill=T, nrows=600) -> brach_morph

#drop species that aren't considered in the CA structure portion
brach_morph <- brach_morph[!brach_morph[,1]=="species2",]
brach_morph <- brach_morph[!brach_morph[,1]=="species3",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f._bonitae",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f..bonitae",]
brach_morph <- brach_morph[!brach_morph[,1]=="apus",]
brach_morph <- brach_morph[!brach_morph[,1]=="miriamae",]

brach_morph$Species <- droplevels(brach_morph$Species)

#drop characters not of interest
brach_morph <- brach_morph[!colnames(brach_morph) %in% "AGD"]
brach_morph <- brach_morph[!colnames(brach_morph) %in% "NL"]


#-----------------------------------------------------------------------------------------
#store nominal characters elsewhere
brach_limbs <- brach_morph[colnames(brach_morph) %in% c("Species", "Limbstate", "Fldig", "Hldig")]
#remove the nominal characters
brach_morph <- brach_morph[!colnames(brach_morph) %in% c("Limbstate", "Fldig", "Hldig")]


#-----------------------------------------------------------------------------------------
#REDUCE THE DATASET

#Exclude any entries that are missing data
brach_morph <- na.omit(brach_morph)

#drop species that aren't considered in the CA structure portion
brach_morph <- brach_morph[!brach_morph[,1]=="species2",]
brach_morph <- brach_morph[!brach_morph[,1]=="species3",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f._bonitae",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f..bonitae",]
#brach_morph <- brach_morph[!brach_morph[,1]=="vermis",]
#brach_morph <- brach_morph[!brach_morph[,1]=="vindumi",]
#brach_morph <- brach_morph[!brach_morph[,1]=="wrighti",]

#check number of entries that remain
length(brach_morph$Species)	#518 #here it is 515, why?

unique(brach_morph$Species) #45 levels, should be 42 since dalawangdaliri and suluensis are omitted; libayani dropped
brach_morph$Species <- droplevels(brach_morph$Species)

#now check to ensure levels were dropped
unique(brach_morph$Species) #45 levels, should be 43 since dalawangdaliri and suluensis are omitted; libayani dropped

#-----------------------------------------------------------------------------------------
#alternative to PC loading methods that uses Mahalanobis D^2 distances, which are adjusted for sample sizes and variances

#complete dataset

#noTL dataset

#reduced morphology dataset
d2.full <- pairwise.mahalanobis(x=brach_morph[,2:10], grouping=brach_morph$Species)

rownames(d2.full$distance) <- unique(brach_morph$Species);colnames(d2.full$distance) <- unique(brach_morph$Species)

#-----------------------------------------------------------------------------------------
#Obtain Community Means and test against nulls


com.morph.disp <- function(empirical, simulated, mahalmatrix){
  
  #create list objects of the intracommunity distances and mean intracommunity distances
  sim.com.distances <- as.list(rep(NA, (nrow(simulated))))
  names(sim.com.distances) <- paste("SIM", rownames(simulated), sep=".")
  sim.com.mdistances <- as.list(rep(NA, (nrow(simulated))))
  names(sim.com.mdistances) <- paste("SIM", rownames(simulated), sep=".")
  
  emp.com.distances <- as.list(rep(NA, (nrow(empirical)+1)))
  #emp.com.distances <- as.list(rep(NA, (nrow(empirical))))
  names(emp.com.distances) <- paste("EMP", rownames(empirical), sep=".")
  emp.com.mdistances <- as.list(rep(NA, (nrow(empirical))))
  names(emp.com.mdistances) <- paste("EMP", rownames(empirical), sep=".")
  
  #populate distances for simulated coms
  for (a in 1: nrow(simulated)){
    singlecom.holder <- names(simulated[a,simulated[a,]==1]) #obtain species in community by name
    #singlecom.distance.holder <- rep(NA, length(ncol(combn(singlecom.holder,2))))  #create list to store all distances
    singlecom.distance.holder <- {}
    #print(singlecom.holder)
      
      for (b in 1:(length(singlecom.holder)-1)){
        for (c in 2: (length(singlecom.holder))){
          if(b<c){
            #print(mahalmatrix[rownames(mahalmatrix) %in% singlecom.holder[b],colnames(mahalmatrix) %in% singlecom.holder[c]])
            singlecom.distance.holder <- c(singlecom.distance.holder, mahalmatrix[rownames(mahalmatrix) %in% singlecom.holder[b],colnames(mahalmatrix) %in% singlecom.holder[c]])
          }
        }
      }
    
    #store the distances for the individual simulated com
    sim.com.distances[[a]] <- singlecom.distance.holder
    sim.com.mdistances[[a]] <- mean(singlecom.distance.holder)
    
  }
  
  #populate distances for empirical coms
  for (d in 1: nrow(empirical)){
    singlecom.holder <- colnames(empirical[d,empirical[d,]==1]) #obtain species in community by name
    #singlecom.distance.holder <- rep(NA, length(ncol(combn(singlecom.holder,2))))  #create list to store all distances
    singlecom.distance.holder <- {}
    
    for (e in 1:(length(singlecom.holder)-1)){
      for (f in 2: (length(singlecom.holder))){
        if(e<f){
          #print(mahalmatrix[rownames(mahalmatrix) %in% singlecom.holder[e],colnames(mahalmatrix) %in% singlecom.holder[f]])
          singlecom.distance.holder <- c(singlecom.distance.holder, mahalmatrix[rownames(mahalmatrix) %in% singlecom.holder[e],colnames(mahalmatrix) %in% singlecom.holder[f]])
        }
      }
    }
    
    #store the distances for the individual empirical com
    emp.com.distances[[d]] <- singlecom.distance.holder
    emp.com.mdistances[[d]] <- mean(singlecom.distance.holder)
    
  }
  
  invisible(new("com_distances", sim.com.distances=sim.com.distances, emp.com.distances=emp.com.distances, sim.com.mdistances=as.vector(sim.com.mdistances, mode="numeric"), emp.com.mdistances=as.vector(emp.com.mdistances, mode="numeric")))		
  
}


#-----------------------------------------------------------------------------------------
#Loop to Run!

emp.com.files <- list.files(path = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern = "_fr.csv$", full.names = TRUE)
emp.com.files.short <- gsub("\\.csv$","",list.files(path = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern = "_fr.csv$", full.names = FALSE))

for(a in 1:length(emp.com.files)){
  #for (a in 1){
  #read csv file for emp_com
  emp.comm <- read.csv(file=paste0(emp.com.files[a]), header=T, row.names=1)
  #emp.com.short <- gsub("\\.csv$","",emp.comm)
  
  #edit the community matrix
  print(emp.com.files.short[a])
  print(paste0("There are ", nrow(emp.comm), " original empirical communities"))
  emp.comm <- emp.comm[lapply(emp.comm, class)!="factor"]
  emp.comm <- emp.comm[colnames(emp.comm)!="species2"]
  emp.comm <- emp.comm[colnames(emp.comm)!="species3",]
  emp.comm <- emp.comm[colnames(emp.comm)!="dalawangdaliri"]
  emp.comm <- emp.comm[colnames(emp.comm)!="c.f..bonitae"]
  emp.comm <- emp.comm[colnames(emp.comm)!="c.f._bonitae"]
  #emp.comm <- emp.comm[colnames(emp.comm)!="vermis"]
  #emp.comm <- emp.comm[colnames(emp.comm)!="vindumi"]
  #emp.comm <- emp.comm[colnames(emp.comm)!="wrighti"]
  emp.comm <- emp.comm[colnames(emp.comm)!="suluensis"]
  #emp.comm <- emp.comm[colnames(emp.comm)!="libayani"]
  emp.comm <- emp.comm[rowSums(emp.comm)>=2,]
  emp.comm <- na.omit(emp.comm)
  print(paste0("There are ", nrow(emp.comm), " empirical communities"))
  
  #generate null coms
  com.simulator(2,5,999,sort(unique(brach_loc[brach_loc$Species %in% unique(brach_morph$Species),"Species"])), writeCSV=F) -> null.coms
  
  #com morph function
  com.morph.holder <- com.morph.disp(emp.comm, null.coms,d2.full$distance)
  
  #store the mean emp.com morphological distances
  emp.com.mean.holder <- com.morph.holder@emp.com.mdistances
  
  emp.com.grandmean.holder <- mean(emp.com.mean.holder)
  
  #create object for p-values per each emp.com
  emp.com.mean.pvalues <- rep(NA, (length(emp.com.mean.holder)+1))
  
  for (b in 1:nrow(emp.comm)){
    #(length(answer@sim.com.mdistances[answer@sim.com.mdistances<=answer@emp.com.mdistances[1]])/length(answer@sim.com.mdistances))
    emp.com.mean.pvalues[b] <- 1-mean(com.morph.holder@sim.com.mdistances > emp.com.mean.holder[b])
  }
  
  #grand mean values
  emp.com.mean.pvalues[(length(emp.com.mean.holder)+1)] <- 1-mean(com.morph.holder@sim.com.mdistances > emp.com.grandmean.holder)
  emp.com.mean.holder <- c(emp.com.mean.holder, emp.com.grandmean.holder)
  
  print(com.morph.holder@sim.com.mdistances)
  print(mean(com.morph.holder@sim.com.mdistances))
  
  #add NA's to the com matrix for final output
  
  emp.comm.holder <- rbind(emp.comm, rep(NA, times=ncol(emp.comm)))
  
  
  emp.com.newdata <- cbind(emp.com.mean.holder, emp.com.mean.pvalues)
  emp.com.newdata <- cbind(emp.comm.holder, emp.com.newdata)
  
  #write to new csv in the new morphology folder
  
  #write.csv(emp.com.newdata, paste0("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/", emp.com.files.short[a], "_morph.csv"))
  
  #noTL
  write.csv(emp.com.newdata, paste0("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/", gsub("_fr","_rednoTL_fr",emp.com.files.short[a]), "_morph.csv"))
  
  #rm(list=c("emp.com.newdata", "emp.com.mean.pvalues", "emp.com.mean.holder", "com.morph.holder", "emp.comm"))
  
}

#comment the rm line to get a null distro
com.morph.disp(emp.comm, null.coms, d2.full$distance)->answer
