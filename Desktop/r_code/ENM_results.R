#prepare a table of test type, omission error, and test metric

#prepare list of directories to review
enm.dirs <- list.dirs("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/ENMS_Fall_2016", recursive = F)
enm.dirs.short <- list.dirs("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/ENMS_Fall_2016",full.names = F, recursive = F)

#rm the ascii dir
enm.dirs <- enm.dirs[-grep("ascii",enm.dirs)]
enm.dirs.short <- enm.dirs.short[-grep("ascii",enm.dirs.short)]
#rm the libayani micro dir
enm.dirs <- enm.dirs[-grep("libayani_micro",enm.dirs)]
enm.dirs.short <- enm.dirs.short[-grep("libayani_micro",enm.dirs.short)]

#read in the summary file for jackknifing method
enm.jackknife <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/ENMS_Fall_2016/all_jackknife_brachymeles_results.csv", header=T)

#aggregator (Species, Method, Omission, AUC, D, BestTest, BestTest Omission, BestTest AUC)
ENMrater <- data.frame(matrix(data=NA, ncol=8, nrow=length(enm.dirs.short)))
colnames(ENMrater) <- c("Species", "Method", "Omission", "AUC", "D", "Best", "Best_Omission", "Best_AUC")

for(a in 1: length(enm.dirs)){
  print(enm.dirs.short[a])
  #add in the species name in ENMrater
  ENMrater$Species[a] <- enm.dirs.short[a]
  #read in result matrix
  maxentresults <- read.csv(paste0(enm.dirs[a], "/maxentResults.csv"), header=T)
  
  #determine model type
  print(nrow(enm.jackknife[enm.jackknife$Species %in% ENMrater$Species[a],]))
  if(nrow(enm.jackknife[enm.jackknife$Species %in% ENMrater$Species[a],]) > 0){
    ENMrater$Method[a] <- "jackknife"
  }
  else {
    ENMrater$Method[a] <- "standard"
  }
  
  #add the Test AUC
  ENMrater$AUC[a] <- maxentresults$Test.AUC[nrow(maxentresults)]
  
  #add the Omission
  #Method: Jacknife
  
  
}