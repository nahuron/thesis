
library(ENMTools)

#compare niche.overlap for niche equivalence test

pvalue.overlap <- function(NET, emp.overlap) {   #emp overlap is a vector of comparisons from niolap in (d,i) format
  sim.overlap.d <- NET$reps.overlap[2:nrow(NET$reps.overlap),1]
  sim.overlap.i <- NET$reps.overlap[2:nrow(NET$reps.overlap),2]
  
  emp.overlap.d <- emp.overlap[1]
  emp.overlap.i <- emp.overlap[2]
  
  pvalues <- rep(NA, times=2)
  names(pvalues) <- c("D","I")
  
  pvalues[1] <- length(sim.overlap.d[sim.overlap.d<=emp.overlap.d])/length(sim.overlap.d)
  pvalues[2] <- length(sim.overlap.i[sim.overlap.i<=emp.overlap.i])/length(sim.overlap.i)
  
  return(pvalues)
  
}

###script to obtain all p-values across taxa included

#obtain list of possible .rData files to use in comparisons

NET.files <- list.files("/Volumes/NAHURON_THESIS/NET2/", pattern=".rda$", full.names = T)
NET.names <- list.files("/Volumes/NAHURON_THESIS/NET2/", pattern=".rda$", full.names = F)
  NET.names <-  gsub("\\.rda$","",NET.names)

#obtain the appropriate niolap file and split it
emp.overlaps <- read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/revised.niche.overlaps.all.csv",row.names=1,sep=",", header=T)


#get the species names from NET value
gsub("NET_","",NET.names)
strsplit(gsub("NET_","",NET.names[1]),split = "_vs_")


#create p-value object

pvalues.net <- as.data.frame(matrix(NA, ncol=4,nrow=length(NET.files)))
colnames(pvalues.net) <- c("Dpvalue","Ipvalue", "EmpD","EmpI")

##loop
#load in sequence each .rData file
#run check.overlap for i and d values and write to p-value object
#remove .rData file
#next iteration

for (a in 1:length(NET.files))  {
#for (a in 1)  {
  load(NET.files[a])
  print(paste0(as.character(NET.files[a])))
  NET.holder <- strsplit(gsub("NET_","",NET.names[a]),split = "_vs_")
  
  empir.i <- emp.overlaps[(rownames(emp.overlaps)==NET.holder[[1]][1]),colnames(emp.overlaps)==NET.holder[[1]][2]]
  empir.d <- emp.overlaps[(rownames(emp.overlaps)==NET.holder[[1]][2]),colnames(emp.overlaps)==NET.holder[[1]][1]]
  
  print(paste0("Empirical D is: ", empir.d," and Empirical I is:", empir.i))
  
  pvalues.net[a,1:2] <- pvalue.overlap(get(NET.names[a]), c(empir.d,empir.i))
  pvalues.net[a,3] <- empir.d
  pvalues.net[a,4] <- empir.i
  rownames(pvalues.net)[a] <- NET.names[a]
  rm(list=as.character(NET.names[a]))
  gc()
}

write.csv(pvalues.net, file="/Volumes/NAHURON_THESIS/pvalues.net.revised8.csv", row.names=T, col.names=T)


#Compare BST niche.overlap

###script to obtain all p-values across taxa included

#obtain list of possible .rData files to use in comparisons

BST.files <- list.files("/Volumes/Nicholas_Huron_Backup/BAST", pattern=".rda$", full.names = T)
BST.names <- list.files("/Volumes/Nicholas_Huron_Backup/BAST", pattern=".rda$", full.names = F)
  BST.names <-  gsub("\\.rda$","",BST.names)

#obtain the appropriate niolap file and split it

emp.overlaps <- read.table("/Users/nicholashuron/Desktop/niche_overlap_tests/communities/revised.niche.overlaps.csv",row.names=1,sep=",", header=T)

#emp.overlaps.d <- emp.overlaps[upper.tri(emp.overlaps, diag=FALSE)]
#emp.overlaps.i <- emp.overlaps[lower.tri(emp.overlaps, diag=FALSE)]

#get the species names from BST value
gsub("BAST_","",BST.names)
strsplit(gsub("BAST_","",BST.names),split = "_vs_")


#create p-value object

pvalues.bst <- as.data.frame(matrix(NA, ncol=4,nrow=length(BST.files)))
colnames(pvalues.bst) <- c("Dpvalue","Ipvalue", "EmpD","EmpI")

##loop
#load in sequence each .rData file
#run check.overlap for i and d values and write to p-value object
#remove .rData file
#next iteration

#for (a in 1:length(BST.files))  {
for (a in 1)  {
  load(BST.files[a])
  print(paste0(as.character(BST.files[a])))
  BST.holder <- strsplit(gsub("BAST_","",BST.names[a]),split = "_vs_")
  BST.holder[[1]][2] <- gsub("bg_", "", BST.holder[[1]][2])
  
  empir.d <- emp.overlaps[(rownames(emp.overlaps)==BST.holder[[1]][1]),colnames(emp.overlaps)==BST.holder[[1]][2]]
  empir.i <- emp.overlaps[(rownames(emp.overlaps)==BST.holder[[1]][2]),colnames(emp.overlaps)==BST.holder[[1]][1]]
  
  print(paste0("Empirical D is: ", empir.d," and Empirical I is:", empir.i))
  
  pvalues.bst[a,1:2] <- pvalue.overlap(get(BST.names[a]), c(empir.d,empir.i))
  pvalues.bst[a,3] <- empir.d
  pvalues.bst[a,4] <- empir.i
  rownames(pvalues.bst)[a] <- BST.names[a]
  rm(list=as.character(BST.names[a]))
  gc()
}

write.csv(pvalues.bst, file="/Volumes/Nicholas_Huron_Backup/BAST/pvalues.bst.csv", row.names=T, col.names=T)