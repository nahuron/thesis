#separate out coordinate data into new files by species



#function form of this script

sp.coordsorter <- function(coords, sp=NULL, writepath=getwd())	{
	
	colnames(coords) <- c("species", colnames(coords)[2], colnames(coords)[3])
	sp.names <- {}
	
	for (i in 1:length(unique(coords[,1])))	{
	
	sp.uni <- unique(coords[,1])
	sp.names <- coords[coords[,1]==sp.uni[i],]
	rownames(sp.names)<-{}
	write.csv(sp.names, file=paste(writepath, sp.uni[i], ".csv", sep=""), row.names=F)
	print(sp.names)
							}
							
																}
																
																															
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/locality_data")
read.csv("Brachymeles.unique.locality.winter2016.csv", header=T)->brach.coords
colnames(brach.coords) <- c("species", "longitude", "latitude")
brach.coords <- brach.coords[,1:3]

sp.coordsorter(coords=brach.coords, writepath= "/Users/nicholashuron/Desktop/niche_overlap_tests/garbage/")


							
#read all files that are background for a species and sample 10,000 points without replacement from those data

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/MAXENT_2016/Maxent_experiments_2016_06_07/")

as.character(list.files())-> foldernames
paste(getwd(),foldernames,sep="/") -> folderpaths
myfiles_comb <- {}				

#enclose in a for loop starting here

for (j in 1:length(folderpaths))	{

list.files(path=folderpaths[j], pattern="backgroundPredictions[.]csv")->temp
temp <- paste(folderpaths[j], temp, sep="/")
myfiles <- lapply(temp, read.csv)
do.call(rbind,myfiles) -> myfiles_comb
myfiles_comb <- myfiles_comb[,1:2]
myfiles_comb <- cbind(rep(as.character(foldernames[j]), times=length(nrow(myfiles_comb))), myfiles_comb)
colnames(myfiles_comb) <- c("species", "longitude", "latitude")
myfiles_comb[,1] <- gsub("_jackknife", "", myfiles_comb[,1])

if(nrow(myfiles_comb)>9999)	{

picklist <- sample(rownames(myfiles_comb), 10000)
myfiles_sampled <- myfiles_comb[picklist,]
							}	else	{
									picklist <- sample(rownames(myfiles_comb), 10000, replace=T)
									myfiles_sampled <- myfiles_comb[picklist,]
										}

rownames(myfiles_sampled) <- {}

write.csv(myfiles_sampled, file=paste("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/MAXENT_2016/niche_overlap_tests/input/", "background_", foldernames[j], ".csv", sep=""), row.names=F)


myfiles_comb <- {}
									}