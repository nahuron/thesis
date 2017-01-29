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

getthatENM(emp.comm,BST_sig_d)