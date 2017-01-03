#PCA by species means (transform, then calculate means)


read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/brachymeles.morphology/huron_brachymeles_morph_raw_means.csv", header=T, fill=T, nrows=600) -> brach_morph

#drop species that aren't considered in the CA structure portion
brach_morph <- brach_morph[!brach_morph[,1]=="species2",]
brach_morph <- brach_morph[!brach_morph[,1]=="species3",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f._bonitae",]
brach_morph <- brach_morph[!brach_morph[,1]=="c.f..bonitae",]
brach_morph <- brach_morph[!brach_morph[,1]=="apus",]
brach_morph <- brach_morph[!brach_morph[,1]=="miriamae",]

brach_morph$Species <- droplevels(brach_morph$Species)

brach_morph <- brach_morph[!colnames(brach_morph) %in% "AGD"]

#Exclude any entries that are missing data
brach_morph <- na.omit(brach_morph)

#create means object
brach_morph_means <- as.data.frame(matrix(NA,nrow=length(unique(brach_morph$Species)), ncol=ncol(brach_morph)))
brach_morph_means[,1] <- unique(brach_morph$Species)
colnames(brach_morph_means)[1] <- colnames(brach_morph)[1]


#script to store means in means object
for ( a in 2:ncol(brach_morph)) {
  means.holder <- tapply(brach_morph[,a], INDEX=brach_morph$Species, FUN=mean)
  brach_morph_means[,a] <- means.holder
  colnames(brach_morph_means)[a] <- colnames(brach_morph)[a]
}
  
brach_morph_means <- na.omit(brach_morph_means)

##standardize ruling:
#divide by HL
#add 1
#ln() transform

