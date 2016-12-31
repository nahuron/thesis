# separate out results of each dataset by faunal regions Mindanao vs Luzon

#requires computing of means?

factor_partition <- function(inputdirectory, inputfilepattern, factorby, outputdirectory=getwd())
{
  #obtain both the full list of file paths and the shortened file names
  full.inputfiles <- list.files(path=inputdirectory, pattern=inputfilepattern, full.names=TRUE)
  short.inputfiles <- list.files(path=inputdirectory, pattern=inputfilepattern, full.names = FALSE)
  
  #remove csv extension from short
  short.inputfiles <- gsub(".csv","",short.inputfiles)
  
  for(a in 1:length(full.inputfiles)){
    #read in full dataset
    fulldata.holder <- read.csv(full.inputfiles[a], header=T, row.names = 1)
    
    for (b in 1:length(factorby)){
      
      #grab the factored column
      subsetdata.holder <- fulldata.holder[sapply(fulldata.holder,is.factor)]
      
      #find the rows that contain the value of b in factorby
      subsetdata.finalholder <- which(subsetdata.holder==factorby[b])
      
      #isolate entries for final partition
      subsetdata.final <- fulldata.holder[subsetdata.finalholder,]
      
      #obtain mean values for metrics of the subset
      
      
      #store final object according to outputdirectory path
      csvname <- paste0(outputdirectory, "/",short.inputfiles[a],"_",factorby[b],".csv")
      write.csv(subsetdata.final, csvname)
      
      csvname <- {}; subsetdata.final <- {}; subsetdata.finalholder <- {}
    }
    
    subsetdata.holder <- {}; fulldata.holder <- {}
    
  }
  
}

#genetic data
factor_partition("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic", inputfilepattern = "_fr.csv", factorby = c("L","M"), outputdirectory = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/com.genetic.fr")

#morphological data
factor_partition("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph", inputfilepattern = "_fr.csv", factorby = c("L","M"), outputdirectory = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets//Morphological/com.morph.fr")

#morphological data with noTL
factor_partition("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.noTL", inputfilepattern = "_fr.csv", factorby = c("L","M"), outputdirectory = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Morphological/com.morph.noTL.fr")

#raw communities
factor_partition("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", inputfilepattern = "_fr.csv", factorby = c("L","M"), outputdirectory = "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2.fr")
