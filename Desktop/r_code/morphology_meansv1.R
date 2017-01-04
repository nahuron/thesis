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
brach_morph <- brach_morph[!colnames(brach_morph) %in% "NL"]


#Exclude any entries that are missing data
brach_morph <- na.omit(brach_morph)
brach_morph$Species <- droplevels(brach_morph$Species)


###standardize ruling:
##divide by HL
#new brach object for transformed values
brach_morph_final <- brach_morph
#store all HL
brach_HL <- brach_morph_final$HL
#rm HL from original object
brach_morph_final <- brach_morph_final[!colnames(brach_morph_final) %in% "HL"]
#divide certain columns by HL (SVL, TL, ForeL, HindL, MBW, MBD, TW, TD, HW, HD, ED, END, SNL, IND)
for(b in 2:ncol(brach_morph_final)){
  if(is.integer(brach_morph_final[,b])==FALSE && is.numeric(brach_morph_final[,b])==TRUE){
    #ln(continuous value then divide by HL then add 1)
    print(paste0("Continuous: ", colnames(brach_morph_final)[b]))
    brach_morph_final[,b] <- log(((brach_morph_final[,b]/brach_HL)+1))
  }
  #ln(discrete value then add 1)
  else if(is.integer(brach_morph_final[,b])==TRUE && is.numeric(brach_morph_final[,b])==TRUE){
    print(paste0("Discrete: ", colnames(brach_morph_final)[b]))
    brach_morph_final[,b] <- log(brach_morph_final[,b]+1)
  }
}

#create means object
brach_morph_means <- as.data.frame(matrix(NA,nrow=length(unique(brach_morph_final$Species)), ncol=ncol(brach_morph_final)))
brach_morph_means[,1] <- unique(brach_morph_final$Species)
colnames(brach_morph_means)[1] <- colnames(brach_morph_final)[1]


#script to store means in means object
for ( a in 2:ncol(brach_morph_final)) {
  means.holder <- tapply(brach_morph_final[,a], INDEX=brach_morph_final$Species, FUN=mean)
  brach_morph_means[,a] <- means.holder
  colnames(brach_morph_means)[a] <- colnames(brach_morph)[a]
}

brach_morph_means <- na.omit(brach_morph_means)

#obtain all loadings regardless of community to establish approximate min and max for figure
pca.all <- prcomp(brach_morph_means[,2:ncol(brach_morph_means)], scale=T)
rownames(pca.all$x) <- sort(unique(brach_morph_means$Species))

summary(pca.all)->pca.varholder
plot(pca.all$x[,1],pca.all$x[,2])

#-----------------------------------------------------------------------------------------

