#PCA by species means (transform, then calculate means)
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
#alternative PCA groupings:
#isolate characters

#noTL dataset (FullNT)
brach_morph <- brach_morph[colnames(brach_morph) %in% c("Species", "SVL","HL", "ForeL", "HindL", "MBW", "MBD", "TW", "TD", "HW", "HD", "ED", "END", "SNL", "IND", "Limbstate", "Fldig", "Hldig")]

#reduced dataset (Body)
brach_morph <- brach_morph[colnames(brach_morph) %in% c("Species", "SVL","HL", "TL", "ForeL", "HindL", "MBW", "TW", "HW", "Limbstate", "Fldig", "Hldig")]

#reduced and noTL dataset (BodyNT)
brach_morph <- brach_morph[colnames(brach_morph) %in% c("Species", "SVL","HL", "ForeL", "HindL", "MBW", "TW", "HW", "Limbstate", "Fldig", "Hldig")]

#no discrete characters (warning: labels will not work)
#store nominal characters elsewhere
brach_limbs <- brach_morph[colnames(brach_morph) %in% c("Species", "Limbstate", "Fldig", "Hldig")]
brach_morph <- brach_morph[!colnames(brach_morph) %in% c("Limbstate", "Fldig", "Hldig")]


#Exclude any entries that are missing data
brach_morph <- na.omit(brach_morph)
brach_morph$Species <- droplevels(brach_morph$Species)
#-----------------------------------------------------------------------------------------

###standardize ruling:
##divide by HL
#new brach object for transformed values
#brach_morph_final <- brach_morph
#store all HL
brach_HL <- brach_morph$HL
#rm HL from original object
brach_morph_final <- brach_morph[!colnames(brach_morph) %in% "HL"]
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

#-----------------------------------------------------------------------------------------
#create means object
brach_morph_means <- as.data.frame(matrix(NA,nrow=length(unique(brach_morph_final$Species)), ncol=ncol(brach_morph_final)))
brach_morph_means[,1] <- unique(brach_morph_final$Species)
colnames(brach_morph_means)[1] <- colnames(brach_morph_final)[1]

#store nominal characters elsewhere
if(any(colnames(brach_morph) %in% "Limbstate")==TRUE){
  brach_limbs <- brach_morph[colnames(brach_morph) %in% c("Species", "Limbstate", "Fldig", "Hldig")]
}
  
#script to store means in means object
for ( a in 2:ncol(brach_morph_final)) {
  means.holder <- tapply(brach_morph_final[,a], INDEX=brach_morph_final$Species, FUN=mean)
  brach_morph_means[,a] <- means.holder
  colnames(brach_morph_means)[a] <- colnames(brach_morph_final)[a]
}

brach_morph_means <- na.omit(brach_morph_means)

#-----------------------------------------------------------------------------------------
#obtain all loadings regardless of community to establish approximate min and max for figure
pca.all <- prcomp(brach_morph_means[,2:ncol(brach_morph_means)], scale=T)
rownames(pca.all$x) <- sort(unique(brach_morph_means$Species))

summary(pca.all)->pca.varholder
#-----------------------------------------------------------------------------------------

library(maptools)
library(scales)
library(TeachingDemos)
library(RColorBrewer)
n <- 40
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
sort(unique(col_vector)) -> col_vector
cols40 <- sort(c(1, 2, 3, 5, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 24, 28, 29, 31, 36, 41, 43, 45, 46, 51, 52, 53, 54, 57, 59, 63, 64, 6, 8, 21, 35, 34, 4, 67))

#create palette
sp.palette <- (col_vector[cols40])[1:length(unique(brach_morph$Species))]

#store colors and species
species.legend <- cbind(unique(as.character(as.factor(brach_morph[,"Species"]))), sp.palette)
#-----------------------------------------------------------------------------------------

#store digit counts for labels
brach_fldig <- tapply(brach_limbs$Fldig, INDEX=brach_limbs$Species, FUN=mean)
brach_fldig <- round(brach_fldig, digits=2)
brach_hldig <- tapply(brach_limbs$Hldig, INDEX=brach_limbs$Species, FUN=mean)
brach_hldig <- round(brach_hldig, digits=2)
brach_limb <- tapply(brach_limbs$Limbstate, INDEX=brach_limbs$Species, FUN=mean)
#-----------------------------------------------------------------------------------------

palette(sp.palette)
plot(pca.all$x[,1],pca.all$x[,2], 
     pch=23, cex=2, col="black", bg=as.factor(rownames(pca.all$x)), 
     xlim=c(min(pca.all$x[,1]*1.25),max(pca.all$x[,1]*1.25)), 
     ylim=c(min(pca.all$x[,2]*1.25),max(pca.all$x[,2]*1.25)), 
     main=expression(Morphometric~PCA~of~Species~of~italic(Brachymeles)), 
     xlab=paste0("PC 1 - ", (pca.varholder$importance[2,1]*100), " %"),
     ylab=paste0("PC 2 - ", (pca.varholder$importance[2,2]*100), " %"))
pointLabel(pca.all$x[,1],pca.all$x[,2],paste0(rownames(pca.all$x), " ",brach_limb," (", brach_fldig, ",", brach_hldig,")"), allowSmallOverlap=F, col="black", cex=0.7, doPlot=FALSE) -> xypts
shadowtext(xypts$x, xypts$y, labels=paste0(rownames(pca.all$x), " ",brach_limb," (", brach_fldig, ",", brach_hldig,")"), col=palette(sp.palette), bg="black", cex=0.7, r=0.05)

