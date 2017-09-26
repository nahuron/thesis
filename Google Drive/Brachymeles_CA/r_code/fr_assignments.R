#script to add FR's to communities

#noTL
setwd("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Morphological/com.morph.noTL/")

#D^2
setwd("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Morphological/com.morphd2.noTL/")

#mac genetic
setwd("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Genetic/com.genetic.3/")

#linux genetic
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic")

#mac raw coms
setwd("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Community/communitiesv3")


#read in example dataset
mydata.filepath <- list.files(path=getwd(), pattern = "\\.csv$", full.names = TRUE)
mydata.filepathshort <- gsub("\\.csv$","",list.files(path=getwd(), pattern = ".csv$", full.names = FALSE))

#read in FR key
brach_fr_key <- read.csv("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Community/brach_FR_key.csv", header=T)

#linux
brach_fr_key <- read.csv("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/brach_FR_key.csv", header=T)


for (c in 1:length(mydata.filepath))  {
  #call data itself
  mydata <- read.csv(mydata.filepath[c], header=T, row.names=1)
  #create object for communities to have FR code
  com.fr <- rep(NA, times=nrow(mydata))
  
#Iterate Matching and coding in com.fr object
  if(any(is.na(mydata[nrow(mydata),]))){
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
    
  }
  else{
    for (a in 1:(nrow(mydata))) { #isolate row of interest
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
    
  }
  
  #write data results
  write.csv(mydata.revised, paste0(getwd(),"/", mydata.filepathshort[c], "_fr.csv"))
  
rm(mydata)
rm(mydata.revised)
rm(com.fr)
}


#function to add FR's to object
add.FR <- function(communities, FR.key){
  #create object for communities to have FR code
  com.fr <- rep(NA, times=nrow(communities))
  
  #Iterate Matching and coding in com.fr object
  if(any(is.na(communities[nrow(communities),]))){
    for (a in 1:(nrow(communities)-1)) { #isolate row of interest
      com.match.hold <- 0
      for (b in which(apply(X=communities, MARGIN=2, FUN=function(x) all(x[!is.na(x)==TRUE]%%1==0)), useNames = FALSE)) { #isolate cell within row of interest
        if (communities[a,b]==1){
          m <- match(as.character(colnames(communities[a,][b])), as.character(FR.key[,1]))
          com.match.hold <- c(com.match.hold, as.character(FR.key$FR.code[m]))
          com.match.hold <- com.match.hold[com.match.hold!=0]
          com.match.hold <- com.match.hold[!is.na(com.match.hold)]
          print(com.match.hold)
        }
      }
      com.match.hold <- unique(com.match.hold)
      print(com.match.hold)
      if (length(com.match.hold) > 1 ) {
        print(paste0("Community ", rownames(communities)[a], " spans more than 1 Faunal Regions!"))
        print(com.match.hold)
        com.fr[a] <- as.character(length(com.match.hold))
      }
      else if (length(com.match.hold == 1)) {com.fr[a] <- as.character(FR.key$FR.code[m])}
      
      mydata.revised <- cbind(communities, com.fr)
      
      return(mydata.revised)
      
      rm(com.match.hold)
    }
    
  }
  else{
    for (a in 1:(nrow(communities))) { #isolate row of interest
      com.match.hold <- 0
      for (b in which(apply(X=communities, MARGIN=2, FUN=function(x) all(x[!is.na(x)==TRUE]%%1==0)), useNames = FALSE)) { #isolate cell within row of interest
        if (communities[a,b]==1){
          m <- match(as.character(colnames(communities[a,][b])), as.character(FR.key[,1]))
          com.match.hold <- c(com.match.hold, as.character(FR.key$FR.code[m]))
          com.match.hold <- com.match.hold[com.match.hold!=0]
          com.match.hold <- com.match.hold[!is.na(com.match.hold)]
          print(com.match.hold)
        }
      }
      com.match.hold <- unique(com.match.hold)
      print(com.match.hold)
      if (length(com.match.hold) > 1 ) {
        print(paste0("Community ", rownames(communities)[a], " spans more than 1 Faunal Regions!"))
        print(com.match.hold)
        com.fr[a] <- as.character(length(com.match.hold))
      }
      else if (length(com.match.hold == 1)) {com.fr[a] <- as.character(FR.key$FR.code[m])}
      
      mydata.revised <- cbind(communities, com.fr)
      
      print(mydata.revised)
      
      rm(com.match.hold)
    }
    
  }
  return(mydata.revised)
  #write data results
}
