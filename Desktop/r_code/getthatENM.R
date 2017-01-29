#determine if communities contain species with ENM's that are unique
#load in community matrix
#load in NIT sig matrix

#create function that creates a list that contains the results for each pairwise comparison in each com

getthatENM <- function(com.matrix, sig.matrix){
  #list of species by community
  species.list <- list()
  length(species.list) <- nrow(com.matrix)
  names(species.list) <- rownames(com.matrix)
  
  #list of community NIT values
  coms.list <- list()
  length(coms.list) <- nrow(com.matrix)
  names(coms.list) <- rownames(com.matrix)
  
  #loop to obtain the species names
  for (a in 1: length(species.list)){
    species.list[[a]] <- colnames(com.matrix[a,which(com.matrix[a,]==1)])
  }
  
  #find the comparisons and store them in coms.list
  for(b in 1:length(species.list)){   #starts by going to com
    coms.list[[b]] <- sig.matrix[rownames(sig.matrix) %in% species.list[[b]],colnames(sig.matrix) %in% species.list[[b]]]
      
    }
  
  return(coms.list)
  
}

getthatENM(emp.comm,NET_sig)


#read in community options
community.files <- list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv3", full.names=T)
community.files <- community.files[-grep("_fr.csv", community.files)]
community.files.short <- gsub("\\.csv$","",list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv3", pattern="\\.csv$", full.names=F))
community.files.short <- community.files.short[-grep("_fr", community.files.short)]

for(n in 1: length(community.files.short)){
  emp.comm <- read.csv(paste0(community.files[n]), header=T, row.names=1)
  #remove species not in ENMs
  emp.comm <- emp.comm[,!colnames(emp.comm) %in% c("c.f._bonitae", "species2","species3")]
  
  getthatENM(emp.comm, NET_sig)
}