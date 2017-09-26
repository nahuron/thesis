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
  
  #function to add PAIC names or identifiers to list element names
  if(grep("com.fr", colnames(com.matrix)) > 0){
    names(coms.list) <- paste0(rownames(com.matrix), "_", com.matrix$com.fr)
  }
  
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

#read in the ENM comparison table
NET_sig <- read.csv("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Geographic/niche_overlap/NET_ID.csv",header = T, row.names = 1)

#read in community options
community.files <- list.files(path= "/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Community/communitiesv3", full.names=T)
community.files <- community.files[grep("_fr.csv", community.files)]
community.files.short <- gsub("\\.csv$","",list.files(path= "/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Community/communitiesv3", pattern="\\.csv$", full.names=F))
community.files.short <- community.files.short[grep("_fr", community.files.short)]

results.net <- list()
length(results.net) <- length(community.files.short)
names(results.net) <- community.files.short

for(n in 1: length(community.files.short)){
  emp.comm <- read.csv(paste0(community.files[n]), header=T, row.names=1)
  #remove species not in ENMs
  emp.comm <- emp.comm[,!colnames(emp.comm) %in% c("c.f._bonitae", "species2","species3")]
  
  results.net[[n]] <- getthatENM(emp.comm, NET_sig)
  
  print(community.files.short[n])
  print(getthatENM(emp.comm, NET_sig))
  rm(emp.comm)
}


#dummy lines for analyzing
for(o in 1:length(results.net)){
  print(lapply(results.net[[o]], function(x) sum(x[lower.tri(x)])))
}
lapply(results.net$revised_com_0.05, function(x) x[upper.tri(x)]<-NA)
lapply(results.net$revised_com_0.05, function(x) x[lower.tri(x)])
lapply(results.net$revised_com_1.00_fr, function(x) sum(x[lower.tri(x)]))

lapply(results.net$revised_com_0.05, function(x) sum(x[upper.tri(x)]))


#loop to pull I values for each community at grid size 0.05
for (b in 1:length(results.net)){
  print(names(results.net)[b])
  NETS <- rep(NA,length(results.net[[b]]))
  for(a in 1: length(results.net[[b]])){
    #print(results.net[[b]][[a]][lower.tri(results.net[[b]][[a]])])
    tester <- results.net[[b]][[a]][lower.tri(results.net[[b]][[a]])]
    #print(paste0(names(results.net)[b],":", a, " is...", (100-(sum(tester[!is.na(tester)])/length(tester[!is.na(tester)])*100)), "% of ENMs similar"))
    NETS[a] <- (100-(sum(tester[!is.na(tester)])/length(tester[!is.na(tester)])*100))
    print(mean(NETS[!is.na(NETS)]))
  }
  #read in each grid size communities matrix and add as a new row
  commatrix <- read.csv(community.files[b], header=T, row.names = 1)
  final.commatrix <- cbind(commatrix[,-ncol(commatrix)], rowSums(commatrix[,-ncol(commatrix)]), NETS, commatrix[,ncol(commatrix)])
  colnames(final.commatrix)[ncol(final.commatrix)] <- "com.fr"
  colnames(final.commatrix)[(ncol(final.commatrix)-2)] <- "n.coms"
  final.commatrix <- rbind(final.commatrix, rep(NA, ncol(final.commatrix)))
  final.commatrix$n.coms[nrow(final.commatrix)] <- mean(final.commatrix$n.coms[!is.na(final.commatrix$n.coms)])
  final.commatrix$NETS[nrow(final.commatrix)] <- mean(final.commatrix$NETS[!is.na(final.commatrix$NETS)])
  
  write.csv(final.commatrix, file=paste0("/Users/nicholashuron/Desktop/Huron, Nicholas/thesis/Datasets/Geographic/com.NET/",community.files.short[b],"_net_i.csv"))
}





###NOT USING!!!!!!
#D values
for (b in 1:length(results.net)){
  print(names(results.net)[b])
  NETS <- rep(NA,length(results.net[[b]]))
  for(a in 1: length(results.net[[b]])){
    #print(results.net[[b]][[a]][lower.tri(results.net[[b]][[a]])])
    tester <- results.net[[b]][[a]][upper.tri(results.net[[b]][[a]])]
    #print(paste0(names(results.net)[b],":", a, " is...", (100-(sum(tester[!is.na(tester)])/length(tester[!is.na(tester)])*100)), "% of ENMs similar"))
    NETS[a] <- (100-(sum(tester[!is.na(tester)])/length(tester[!is.na(tester)])*100))
    print(mean(NETS[!is.na(NETS)]))
  }
  #read in each grid size communities matrix and add as a new row
  commatrix <- read.csv(community.files[b], header=T, row.names = 1)
  final.commatrix <- cbind(commatrix[,-ncol(commatrix)], rowSums(commatrix[,-ncol(commatrix)]), NETS, commatrix[,ncol(commatrix)])
  colnames(final.commatrix)[ncol(final.commatrix)] <- "com.fr"
  colnames(final.commatrix)[(ncol(final.commatrix)-2)] <- "n.coms"
  write.csv(final.commatrix, file="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/com.NET/com_NET_d.csv")
}

