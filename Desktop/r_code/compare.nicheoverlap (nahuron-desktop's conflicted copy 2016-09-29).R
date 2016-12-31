#compare niche.overlap

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
#mac
NET.files <- list.files("/Volumes/Nicholas_Huron_Backup/NIT", pattern=".rda$", full.names = T)
NET.names <- list.files("/Volumes/Nicholas_Huron_Backup/NIT", pattern=".rda$", full.names = F)
  NET.names <-  gsub("\\.rda$","",NET.names)

#linux
  NET.files <- list.files("/media/nicholas/Nicholas_Huron_Backup/NIT", pattern=".rda$", full.names = T)
  NET.names <- list.files("/media/nicholas/Nicholas_Huron_Backup/NIT", pattern=".rda$", full.names = F)
  NET.names <-  gsub("\\.rda$","",NET.names)
  
#obtain the appropriate niolap file and split it

#mac
emp.overlaps <- read.table("/Users/nicholashuron/Desktop/niche_overlap_tests/communities/revised.niche.overlaps.csv",row.names=1,sep=",", header=T)

#backup
emp.overlaps <- read.table("/media/nicholas/Nicholas_Huron_Backup/niche_overlap_tests/communities/revised.niche.overlaps.csv",row.names=1,sep=",", header=T)

#emp.overlaps.d <- emp.overlaps[upper.tri(emp.overlaps, diag=FALSE)]
#emp.overlaps.i <- emp.overlaps[lower.tri(emp.overlaps, diag=FALSE)]

library(ENMTools)

#get the species names from NET value
gsub("NET_","",NET.names)
strsplit(gsub("NET_","",NET.names[1]),split = "_vs_")


#create p-value object

pvalues <- as.data.frame(matrix(NA, ncol=2,nrow=length(NET.files)))

##loop
#load in sequence each .rData file
#run check.overlap for i and d values and write to p-value object
#remove .rData file
#next iteration

for (a in 1:length(NET.files))  {
#for (a in 1:5)  {
  load(NET.files[a])
  NET.holder <- strsplit(gsub("NET_","",NET.names[a]),split = "_vs_")
  
  empir.d <- emp.overlaps[(rownames(emp.overlaps)==NET.holder[[1]][1]),colnames(emp.overlaps)==NET.holder[[1]][2]]
  empir.i <- emp.overlaps[(rownames(emp.overlaps)==NET.holder[[1]][2]),colnames(emp.overlaps)==NET.holder[[1]][1]]
  
  pvalues[a,1:2] <- pvalue.overlap(get(NET.names[a]), c(empir.d,empir.i))
  rm(list=as.character(NET.names[a]))
  gc()
}

