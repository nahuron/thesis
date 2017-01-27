#script to add FR's to communities

#mac morphological
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph/")

#noTL
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.noTL/")

#D^2
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morphd2.noTL/")

#mac genetic
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

#linux genetic
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

#mac raw coms
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv3")


#read in example dataset
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)

#linux
brach_fr_key <- read.csv("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)


for (c in 1:length(mydata.filepath))  {
  #call data itself
  mydata <- read.csv(mydata.filepath[c], header=T, row.names=1)
  #create object for communities to have FR code
  com.fr <- rep(NA, times=nrow(mydata))
  
#Iterate Matching and coding in com.fr object
for (a in 1:(nrow(mydata)-1)) { #isolate row of interest
    com.match.hold <- 0
  for (b in which(apply(X=mydata, MARGIN=2, FUN=function(x) all(x[!is.na(x)==TRUE]%%1==0)), useNames = FALSE)) { #isolate cell within row of interest
  if (mydata[a,b]==1){
    m <- match(as.character(colnames(mydata[a,][b])), as.character(brach_fr_key[,1]))
    com.match.hold <- c(com.match.hold, as.character(brach_fr_key$FR.code[m]))
    com.match.hold <- com.match.hold[com.match.hold!=0]
    com.match.hold <- com.match.hold[!is.na(com.match.hold)]
    print(com.match.hold)
  }
  }
  com.match.hold <- unique(com.match.hold)
  print(com.match.hold)
  if (length(com.match.hold) > 1 ) {
    print(paste0("Community ", rownames(mydata)[a], " spans more than 1 Faunal Regions!"))
    print(com.match.hold)
    com.fr[a] <- as.character(length(com.match.hold))
    }
  else if (length(com.match.hold == 1)) {com.fr[a] <- as.character(brach_fr_key$FR.code[m])}
  
  mydata.revised <- cbind(mydata, com.fr)
  
  print(mydata)
  
  rm(com.match.hold)
}
  
  #write data results
  write.csv(mydata.revised, paste0(getwd(),"/", mydata.filepathshort[c], "_fr.csv"))
  
rm(mydata)
rm(mydata.revised)
rm(com.fr)
}

