#necessary packages
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
library(viridis)

#read in empirical overlap data
read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/revised.niche.overlaps.all.csv", header=T, row.names = 1) -> overlaps

#d statistic only
overlaps_d <- overlaps
overlaps_d[lower.tri(overlaps_d)==TRUE]<-NA

#i statistic only
overlaps_i <- overlaps
overlaps_i[upper.tri(overlaps_i)==TRUE]<-NA

rnames <- rownames(overlaps)
cnames <- colnames(overlaps)

my_palette <- viridis_pal(1,0,1,"D")
heatmap.2(as.matrix(overlaps), col = my_palette, Rowv=F, Colv = F, trace="none", density.info = "histogram")


#read in significance BST_D
read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/BST_D.csv", header=T, row.names = 1) -> overlaps_sigd

rnames <- rownames(overlaps_sigd)
cnames <- colnames(overlaps_sigd)

my_palette <- viridis_pal(1, 0,1, "D")
heatmap.2(as.matrix(overlaps_sigd), col=my_palette, Rowv = F, Colv = F, trace = "none", density.info = "histogram")

#read in significance BST_I
read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/BST_I.csv", header=T, row.names = 1) -> overlaps_sigi

rnames <- rownames(overlaps_sigi)
cnames <- colnames(overlaps_sigi)

my_palette <- viridis_pal(1, 0,1, "D")
heatmap.2(as.matrix(overlaps_sigi), col=my_palette, Rowv = F, Colv = F, trace = "none", density.info = "histogram")


#read in significance NET_ID
read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/NET_ID.csv", header=T, row.names = 1) -> overlaps_sigid

rnames <- rownames(overlaps_sigid)
cnames <- colnames(overlaps_sigid)

my_palette <- viridis_pal(1, 0,1, "D")
heatmap.2(as.matrix(overlaps_sigid), col=my_palette, Rowv = F, Colv = F, trace = "none", density.info = "histogram")


