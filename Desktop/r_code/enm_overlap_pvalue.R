#mac
setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic")

#linux
setwd("/home/nicholas/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic")

#BST (row is sp a, col is background of sp b)
pvalues.bst <- read.csv("pvalues.bst.csv", header=T)

BST_names <- pvalues.bst[,1]

sp_names <- unique(unlist(strsplit(gsub("BAST_","",BST_names),split = "_vs_")))
sp_names <- unique(gsub("bg_", "", sp_names))

BST_sig_d <- matrix(NA, nrow=length(sp_names), ncol=length(sp_names))
BST_sig_i <- matrix(NA, nrow=length(sp_names), ncol=length(sp_names))

colnames(BST_sig_d) <- sp_names
rownames(BST_sig_d) <- sp_names
colnames(BST_sig_i) <- sp_names
rownames(BST_sig_i) <- sp_names

#loop that determines the significance and codes the correct cell as -1, 0, or 1
for (a in 1:length(BST_names))  {
#  for (a in 1)  {
  
  #isolate the pair of interest
  sp_name_holder <- unlist(strsplit(gsub("BAST_","",BST_names[a]),split = "_vs_"))
  print(sp_name_holder)
  #remove the bg from the second species
  sp_name_holder <- unique(gsub("bg_", "", sp_name_holder))
  
  results_holder <- pvalues.bst[pvalues.bst[,1]==BST_names[a],]
  print(results_holder)
  
  #D metric
  if(results_holder$Dpvalue >= 0.975)	{
    BST_sig_d[rownames(BST_sig_d)==sp_name_holder[1], (colnames(BST_sig_d)==sp_name_holder[2])] <- 1
  }
  else if(results_holder$Dpvalue <= 0.025)	{
    BST_sig_d[rownames(BST_sig_d)==sp_name_holder[1], (colnames(BST_sig_d)==sp_name_holder[2])] <- -1
  }
  else if(results_holder$Dpvalue >= 0.025 & results_holder$Dpvalue <= 0.975)	{
    BST_sig_d[rownames(BST_sig_d)==sp_name_holder[1], (colnames(BST_sig_d)==sp_name_holder[2])] <- 0
  }
  
  write.csv(BST_sig_d, file="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/BST_D.csv")

  #I metric
  if(results_holder$Ipvalue >= 0.975)	{
    BST_sig_i[rownames(BST_sig_i)==sp_name_holder[1], (colnames(BST_sig_i)==sp_name_holder[2])] <- 1
  }
  else if(results_holder$Ipvalue <= 0.025)	{
    BST_sig_i[rownames(BST_sig_i)==sp_name_holder[1], (colnames(BST_sig_i)==sp_name_holder[2])] <- -1
  }
  else if(results_holder$Ipvalue >= 0.025 & results_holder$Ipvalue <= 0.975)	{
    BST_sig_i[rownames(BST_sig_i)==sp_name_holder[1], (colnames(BST_sig_i)==sp_name_holder[2])] <- 0
  }

write.csv(BST_sig_i, file="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/BST_I.csv")

    
}

#plots for BST
#first loop to read in first dataset
{  pvaluesorted <- unlist(lapply(pvalues.bst[-as.numeric(nrow(pvalues.bst)),"Dpvalue"], function (x) x <= 0.025 | x >= 0.975))
  
  sig.hold <- NA
  nonsig.hold <- NA
  
  for(b in 1:length(pvaluesorted)){
    if (pvaluesorted[b]==TRUE){
      sig.hold <- c(sig.hold, b)
      sig.hold <- sig.hold[!is.na(sig.hold)]
    } 
    else  {
      nonsig.hold <- c(nonsig.hold, b)
      nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
    }
  }

  sigs <- density(pvalues.bst$EmpD[sig.hold], kernel="triangular", from=0, to=1)
  nonsigs <- density(pvalues.bst$EmpD[nonsig.hold], kernel="triangular", from=0, to=1)
  print(pvalues.bst$EmpD[sig.hold])
  print(pvalues.bst$EmpD[nonsig.hold])
  
  hist(pvalues.bst$EmpD[sig.hold], col=rgb(0,1,1,0.8), freq=FALSE, ylim=c(0,2.5), main="Empirical Overlap of Schoener's D", xlab="Overlap")
  hist(pvalues.bst$EmpD[nonsig.hold], col=rgb(1,0,0,0.75), add=TRUE, freq=FALSE)
  lines(nonsigs, lty=5)
  lines(sigs)
  
  plot(sigs, ylim=c(0,2.5))
  lines(nonsigs, lty=5)
  

  pvaluesorted <- unlist(lapply(pvalues.bst[-as.numeric(nrow(pvalues.bst)),"Ipvalue"], function (x) x <= 0.025 | x >= 0.975))
  
  sig.hold <- NA
  nonsig.hold <- NA
  
  for(b in 1:length(pvaluesorted)){
    if (pvaluesorted[b]==TRUE){
      sig.hold <- c(sig.hold, b)
      sig.hold <- sig.hold[!is.na(sig.hold)]
    } 
    else  {
      nonsig.hold <- c(nonsig.hold, b)
      nonsig.hold <- nonsig.hold[!is.na(nonsig.hold)]
    }
  }
  
  sigs <- density(pvalues.bst$EmpI[sig.hold], kernel="triangular", from=0, to=1)
  nonsigs <- density(pvalues.bst$EmpI[nonsig.hold], kernel="triangular", from=0, to=1)
  print(pvalues.bst$EmpI[sig.hold])
  print(pvalues.bst$EmpI[nonsig.hold])
  
  plot(sigs, ylim=c(0,2.5))
  lines(nonsigs, lty=5)

hist(pvalues.bst$EmpI[sig.hold], col=rgb(0,1,1,0.8), freq=FALSE, ylim=c(0,2.5), main="Empirical Overlap of Hellinger-Based I", xlab="Overlap")
hist(pvalues.bst$EmpI[nonsig.hold], col=rgb(1,0,0,0.75), add=TRUE, freq=FALSE)
lines(nonsigs, lty=5)
lines(sigs)

}

  
#NST
  pvalues.nst <- read.csv("/Volumes/NAHURON_THESIS/pvalues.net.revised1.csv", header=T)
  
  NET_names <- sort(pvalues.nst[,1])
  
  sp_names <- sort(unique(unlist(strsplit(gsub("NET_","",NET_names),split = "_vs_"))))
  
  NET_sig <- matrix(NA, nrow=length(sp_names), ncol=length(sp_names))
  colnames(NET_sig) <- sp_names
  rownames(NET_sig) <- sp_names
  
  #loop!
  for (a in 1:length(NET_names))  {
    #for (a in 15)  {
    
    #isolate the pair of interest
    sp_name_holder <- unlist(strsplit(gsub("NET_","",NET_names[a]),split = "_vs_"))
    print(sp_name_holder)
    
    results_holder <- pvalues.nst[pvalues.nst[,1]==NET_names[a],]
    print(results_holder)
    
    #D metric
    if(results_holder$Dpvalue <= 0.025)	{   #ENMs are less similar than the null
      NET_sig[(rownames(NET_sig)==sp_name_holder[2]), (colnames(NET_sig)==sp_name_holder[1])] <- 1
    }
    else if(results_holder$Dpvalue >= 0.025)	{ #ENMs are equivalent to one another due to high overlap
      NET_sig[(rownames(NET_sig)==sp_name_holder[2]), (colnames(NET_sig)==sp_name_holder[1])] <- 0
    }
    
    #I metric
    if(results_holder$Ipvalue <= 0.025){
      NET_sig[(rownames(NET_sig)==sp_name_holder[1]), (colnames(NET_sig)==sp_name_holder[2])] <- 1
    }
    else if(results_holder$Ipvalue >= 0.025)	{
     NET_sig[(rownames(NET_sig)==sp_name_holder[1]), (colnames(NET_sig)==sp_name_holder[2])] <- 0
    }
  }

  write.csv(NET_sig, file="/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Geographic/NET_ID.csv")
  
  
#NIT   
gsub("_vs_","_", pvalues.nst[,1])

NET.files <- list.files("/Volumes/NAHURON_THESIS/NET", pattern=".rda$", full.names = T)
NET.names <- list.files("/Volumes/NAHURON_THESIS/NET", pattern=".rda$", full.names = F)
  NET.names <-  gsub("\\.rda$","",NET.names)
  
  #get the species names from NET value
gsub("NET_","",NET.names)
strsplit(gsub("NET_","",NET.names[1]),split = "_vs_")
