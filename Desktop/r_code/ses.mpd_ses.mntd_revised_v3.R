ses.mpd2 <- function (samp, dis, subsetted=colnames(dis), null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap", "other"), abundance.weighted = FALSE, runs = 999, iterations = 1000, comminmax=NULL, alpha=0.05) 
{
    dis <- as.matrix(dis)
    mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted)
    mean.mpd.obs <- mean(mpd.obs)
    null.model <- match.arg(null.model)
    mpd.rand <- switch(null.model, 
    	taxa.labels = t(replicate(runs, mpd(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted))), 
        richness = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))), 
        frequency = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted))), 
        sample.pool = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))), 
        phylogeny.pool = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "richness"), taxaShuffle(dis), abundance.weighted))), 
        independentswap = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted))), trialswap = t(replicate(runs, mpd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted))),
        other = t(replicate(runs, mpd(as.matrix(com.simulator(comminmax[1], comminmax[2], nrow(samp), subsetted)), dis, abundance.weighted))))
        #print(mpd.rand)
    mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
    mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2, FUN = rank)[1, ]
    mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
    quantiles.null <- quantile(rowMeans(mpd.rand), probs = c(alpha/2, 1 - (alpha/2)))
    
    mpd.rand.z <- (rowMeans(mpd.rand)- mean(mpd.rand.mean))/mean(mpd.rand.sd)
    
    print(rowMeans(mpd.rand))
    print(mpd.rand.mean)
    print(mpd.rand.sd)
    
    return(list(null.means=rowMeans(mpd.rand), mean.mpd.obs= mean.mpd.obs, quantiles.null = quantiles.null, data.frame(ntaxa = specnumber(samp), mpd.obs, mpd.rand.mean, mpd.rand.sd, mpd.obs.rank, mpd.obs.z, mpd.obs.p = mpd.obs.rank/(runs + 1), runs = runs, row.names = row.names(samp)), mean.mpd.obs.z = mean(mpd.obs.z), mpd.rand.z = mpd.rand.z))
}

ses.mntd2 <- function (samp, dis, subsetted=colnames(dis), null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap", "other"), abundance.weighted = FALSE, runs = 999, iterations = 1000, comminmax=NULL, alpha=0.05) 
{
    dis <- as.matrix(dis)
    mntd.obs <- mntd(samp, dis, abundance.weighted)
    mean.mntd.obs <- mean(mntd.obs)
    null.model <- match.arg(null.model)
    mntd.rand <- switch(null.model, 
    taxa.labels = t(replicate(runs, mntd(samp, taxaShuffle(dis), abundance.weighted))), 
    richness = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))), 
    frequency = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted))), 
    sample.pool = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))), 
    phylogeny.pool = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "richness"), taxaShuffle(dis), abundance.weighted))), 
    independentswap = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted))), 
    trialswap = t(replicate(runs, mntd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted))),
    other = t(replicate(runs, mntd(as.matrix(com.simulator(comminmax[1], comminmax[2], nrow(samp), subsetted)), dis, abundance.weighted)))	)
    mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
    mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
    mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
    mntd.obs.rank <- apply(X = rbind(mntd.obs, mntd.rand), MARGIN = 2, FUN = rank)[1, ]
    mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)
    quantiles.null <- quantile(rowMeans(mntd.rand), probs = c(alpha/2, 1 - (alpha/2)))
    
    mntd.rand.z <- (rowMeans(mntd.rand)- mean(mntd.rand.mean))/mean(mntd.rand.sd)
    
    return(list(null.means=rowMeans(mntd.rand), mean.mntd.obs= mean.mntd.obs, quantiles.null = quantiles.null, data.frame(ntaxa = specnumber(samp), mntd.obs, mntd.rand.mean, mntd.rand.sd, mntd.obs.rank, mntd.obs.z, mntd.obs.p = mntd.obs.rank/(runs + 1), runs = runs, row.names = row.names(samp)), mean.mntd.obs.z = mean(mntd.obs.z), mntd.rand.z = mntd.rand.z))
}