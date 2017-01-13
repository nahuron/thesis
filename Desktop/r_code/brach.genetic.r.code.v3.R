#megaloop or lapply

#read in community options
community.files <- list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", full.names=T)
community.files <- community.files[-grep("_fr.csv", community.files)]
community.files.short <- gsub("\\.csv$","",list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2", pattern="\\.csv$", full.names=F))
community.files.short <- community.files.short[-grep("_fr", community.files.short)]

#Luzon only
community.files <- list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2.fr", pattern="L.csv", full.names=T)
community.files.short <- gsub("\\.csv$","",list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2.fr", pattern="L.csv", full.names=F))

#Mindanao only
community.files <- list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2.fr", pattern="M.csv", full.names=T)
community.files.short <- gsub("\\.csv$","",list.files(path= "/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Community/communitiesv2.fr", pattern="M.csv", full.names=F))


#create results object
#tests to run: ses.mpd, ses.mntd, psv all with new null community method
#object should have individual empirical community metric value and p-value (n rows by 6 columns)

	results.df <- data.frame(matrix(data=NA, ncol=7, nrow=50))
	colnames(results.df) <- c("ses.mpd", "ses.mpd.p", "ses.mntd", "ses.mntd.p", "psv", "psv.p", "n.coms")
	results.list <- list()
	for(b in 1:length(community.files.short)){results.list[[b]] <- results.df}
	rm(results.df)

for (a in 1:length(community.files))	{
#for (a in 1:1)	{
	#read in community replicate
  emp.comm <- read.csv(paste0(community.files[a]), header=T, row.names=1)
  emp.comm <- emp.comm[lapply(emp.comm, class)!="factor"]
	print(paste0("There are ", nrow(emp.comm), " original empirical communities"))
	
	#remove species without phy data
	emp.comm <- emp.comm[colnames(emp.comm)!="species2"]; emp.comm <- emp.comm[colnames(emp.comm)!="dalawangdaliri"]; emp.comm <- emp.comm[colnames(emp.comm)!="c.f..bonitae"];emp.comm <- emp.comm[colnames(emp.comm)!="vermis"];emp.comm <- emp.comm[colnames(emp.comm)!="vindumi"]; emp.comm <- emp.comm[colnames(emp.comm)!="wrighti"];emp.comm <- emp.comm[colnames(emp.comm)!="suluensis"];#emp.comm <- emp.comm[colnames(emp.comm)!="libayani"]
	emp.comm <- emp.comm[rowSums(emp.comm)>=2,]
	print(paste0("There are ", nrow(emp.comm), " empirical communities"))
	print(emp.comm)
	
	#add the number of species per community to each row in output list
	results.list[[a]]$n.coms[1:nrow(emp.comm)] <- rowSums(emp.comm)
	print(rowSums(emp.comm))
	
	#conduct SES.MPD with modified null
	if(length(grep("_L", community.files.short[a]))>=1) {
	  ses.mpd.hold <- ses.mpd2(samp=emp.comm, dis=phydist_brach, subsetted=as.character(brach_fr_key_l$Species), null.model = "other", abundance.weighted = FALSE, runs = 1000, comminmax=c(2,5), alpha=0.05)
	}
	else if(length(grep("_M", community.files.short[a]))>=1)  {
	  ses.mpd.hold <- ses.mpd2(emp.comm, phydist_brach, subsetted=as.character(brach_fr_key_m$Species), null.model = "other", abundance.weighted = FALSE, runs = 1000, comminmax=c(2,5), alpha=0.05)
	}
	else{
	  ses.mpd.hold <- ses.mpd2(emp.comm, phydist_brach, null.model = "other", abundance.weighted = FALSE, runs = 1000, comminmax=c(2,5), alpha=0.05)
	}
	
	#populate output list with SES.MPD and p-value for each community
	results.list[[a]]$ses.mpd[1:length(ses.mpd.hold[[4]]$mpd.obs)] <- ses.mpd.hold[[4]]$mpd.obs.z
	results.list[[a]]$ses.mpd.p[1:length(ses.mpd.hold[[4]]$mpd.obs.p)] <- ses.mpd.hold[[4]]$mpd.obs.p
	
	#conduct SES.MNTD with modified null
	if(length(grep("_L", community.files.short[a]))>=1) {
	  ses.mntd.hold <- ses.mntd2(emp.comm, phydist_brach,subsetted=as.character(brach_fr_key_l$Species), null.model = "other", abundance.weighted = FALSE, runs = 1000, comminmax=c(2,5), alpha=0.05)
	}
	else if(length(grep("_M", community.files.short[a]))>=1)  {
	  ses.mntd.hold <- ses.mntd2(emp.comm, phydist_brach,subsetted=as.character(brach_fr_key_m$Species), null.model = "other", abundance.weighted = FALSE, runs = 1000, comminmax=c(2,5), alpha=0.05)
	}
	else{
	  ses.mntd.hold <- ses.mntd2(emp.comm, phydist_brach, null.model = "other", abundance.weighted = FALSE, runs = 1000, comminmax=c(2,5), alpha=0.05)
	}
	
	#populate output list with SES.MNTD and p-value for each community
	results.list[[a]]$ses.mntd[1:length(ses.mntd.hold[[4]]$mntd.obs)] <- ses.mntd.hold[[4]]$mntd.obs.z
	results.list[[a]]$ses.mntd.p[1:length(ses.mntd.hold[[4]]$mntd.obs.p)] <- ses.mntd.hold[[4]]$mntd.obs.p
	
	#calculate PSV for empirical community and populate output list
	psv.bycom.hold <- psv(emp.comm, brach)
	results.list[[a]]$psv[1:length(psv.bycom.hold$PSVs)] <- psv.bycom.hold$PSVs
	
	#obtain null distrubtion to compare to empirical PSVs
	if(length(grep("_L", community.files.short[a]))>=1) {
	  psv.hold <- phylostruct.rev(samp=emp.comm, tree=brach, comminmax=c(2,5), env=NULL, sp=as.character(brach_fr_key_l$Species), metric="psv", null.model="other", runs=1000, alpha= 0.05)
	}
	else if(length(grep("_M", community.files.short[a]))>=1)  {
	  psv.hold <- phylostruct.rev(emp.comm, brach, c(2,5), env=NULL, sp=as.character(brach_fr_key_m$Species), metric="psv", null.model="other", runs=1000, alpha= 0.05)
	}
	else{
	  psv.hold <- phylostruct.rev(emp.comm, brach, c(2,5), env=NULL, metric="psv", null.model="other", runs=1000, alpha= 0.05)
	}
	
	  psv.pvalues.hold <- rep(NA, length(psv.bycom.hold$PSVs))
	  
	#loop to generate p-values from null distrubtion and empirical values
		for (b in 1: length(psv.pvalues.hold))	{ psv.pvalues.hold[b] <- 1-mean(psv.hold$null.means > psv.bycom.hold$PSVs[b])}
	results.list[[a]]$psv.p[1:length(psv.pvalues.hold)] <- psv.pvalues.hold
	
	#store mean values of metrics at bottom of element for output list
	rownames(results.list[[a]][1+nrow(emp.comm),]) <- "Means"
	
	#ses.mpd mean
	results.list[[a]]$ses.mpd[1+nrow(emp.comm)] <- ses.mpd.hold$mean.mpd.obs.z
	
	#ses.mntd mean
	results.list[[a]]$ses.mntd[1+nrow(emp.comm)] <- ses.mntd.hold$mean.mntd.obs.z

	#psv mean
	results.list[[a]]$psv[1+nrow(emp.comm)] <- psv.hold$mean.obs
	
	#species richness mean
	results.list[[a]]$n.coms[1+nrow(emp.comm)] <- mean(rowSums(emp.comm))
	
	#ses.mpd.p for ses.mpd mean
	results.list[[a]]$ses.mpd.p[1+nrow(emp.comm)] <- 1 - mean(ses.mpd.hold$mpd.rand.z > ses.mpd.hold$mean.mpd.obs.z)
	
	#ses.mntd.p for ses.mntd mean
	results.list[[a]]$ses.mntd.p[1+nrow(emp.comm)] <- 1 - mean(ses.mntd.hold$mntd.rand.z > ses.mntd.hold$mean.mntd.obs.z)
	
	#psv.p for psv mean
	results.list[[a]]$psv.p[1+nrow(emp.comm)] <- 1-mean(psv.hold$null.means > psv.hold$mean.obs)
	
	#wipe leftover elements that contain "NA"
	results.list[[a]] <- results.list[[a]][!rowSums(is.na(results.list[[a]])) > 0,]
	print(results.list[[a]][!rowSums(is.na(results.list[[a]])) > 0,])
	
	rm(list=c("emp.comm", "ses.mpd.hold", "ses.mntd.hold", "psv.hold", "psv.bycom.hold", "psv.pvalues.hold"))
										}

for (c in 1:length(results.list))	{
#for (c in 1)	{
  emp.comm <- read.csv(paste0(community.files[c]), header=T, row.names=1)
  emp.comm <- emp.comm[lapply(emp.comm, class)!="factor"]
  #remove species without phy data
  emp.comm <- emp.comm[colnames(emp.comm)!="species2"]; emp.comm <- emp.comm[colnames(emp.comm)!="dalawangdaliri"]; emp.comm <- emp.comm[colnames(emp.comm)!="c.f..bonitae"];emp.comm <- emp.comm[colnames(emp.comm)!="vermis"];emp.comm <- emp.comm[colnames(emp.comm)!="vindumi"]; emp.comm <- emp.comm[colnames(emp.comm)!="wrighti"];emp.comm <- emp.comm[colnames(emp.comm)!="suluensis"];emp.comm <- emp.comm[colnames(emp.comm)!="libayani"]
  emp.comm <- emp.comm[rowSums(emp.comm)>=2,]
  
  results.list[[c]] <- cbind(rbind(emp.comm, rep(NA, times=ncol(emp.comm))), results.list[[c]])
  
  csvname <- paste0("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Genetic/",community.files.short[c], "_phylo.csv")
	write.csv(results.list[[c]], csvname)
									}

#PSV (closer to 0 indicates relatedness and clustering)
psv.result <- psv(comm, brach, compute.var=TRUE)

