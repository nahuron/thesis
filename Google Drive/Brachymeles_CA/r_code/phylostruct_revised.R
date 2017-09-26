phylostruct.rev <- function (samp, tree, comminmax=NULL, env = NULL, metric = c("psv", "psr", "pse", "psc", "sppregs"), null.model = c("frequency", "richness", "independentswap", "trialswap", "other"), runs = 100, it = 1000, alpha = 0.05, fam = "binomial") 
{
	metric <- match.arg(metric)
    null.model <- match.arg(null.model)
    if (metric == "sppregs") {
        nulls <- t(replicate(runs, sppregs(randomizeMatrix(samp, null.model = null.model, iterations = it), env, tree, fam = fam)$correlations))
        obs <- sppregs(samp, env, tree, fam = fam)$correlations
        mean.null <- apply(nulls, 2, mean)
        quantiles.null <- t(apply(nulls, 2, quantile, probs = c(alpha/2, 1 - (alpha/2))))
        if ((null.model != "independentswap") && (null.model !=  "trialswap")) {
            it = NA
        }
        return(list(metric = metric, null.model = null.model, runs = runs, it = it, obs = obs, mean.null = mean.null, quantiles.null = quantiles.null, phylo.structure = NULL, nulls = nulls))
    }
    else 	{
    			if (missing(comminmax))	{
    										nulls <- switch(metric, psv = replicate(runs, mean(psv(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree, compute.var = FALSE)[, 1], na.rm = TRUE)), psr = replicate(runs, mean(psr(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree, compute.var = FALSE)[, 1], na.rm = TRUE)), pse = replicate(runs, mean(pse(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree)[, 1], na.rm = TRUE)), psc = replicate(runs, mean(psc(as.matrix(randomizeMatrix(samp, null.model = null.model, iterations = it)), tree)[, 1], na.rm = TRUE)))
        									quantiles.null <- quantile(nulls, probs = c(alpha/2, 1 - (alpha/2)))
        									mean.null <- mean(nulls)
        									mean.obs <- switch(metric, psv = mean(psv(samp, tree, compute.var = FALSE)[, 1], na.rm = TRUE), psr = mean(psr(samp, tree, compute.var = FALSE)[, 1], na.rm = TRUE), pse = mean(pse(samp, tree)[, 1], na.rm = TRUE), psc = mean(psc(samp, tree)[, 1], na.rm = TRUE))
        									if (mean.obs <= quantiles.null[1]) {
            								phylo.structure = "underdispersed"
        									}
        									else {
            									if (mean.obs >= quantiles.null[2]) {
                									phylo.structure = "overdispersed"
            									}
            								else {
                								phylo.structure = "random"
            									}
        									}
        									if ((null.model != "independentswap") && (null.model != "trialswap")) {
            									it = NA
        									}
        								return(list(metric = metric, null.model = null.model, runs = runs, it = it, mean.obs = mean.obs, mean.null = mean.null, quantiles.null = quantiles.null, phylo.structure = phylo.structure, null.means = nulls))		
    									}
    			else if(null.model == "other")	{
    					Species <- unique(tree$tip.label)
    					#null.coms <- com.simulator(comminmax[1], comminmax[2], nrow(samp), Species)
    					nulls <- switch(metric, psv = replicate(runs, mean(psv(as.matrix(com.simulator(comminmax[1], comminmax[2], nrow(samp), Species)), tree, compute.var = FALSE)[, 1], na.rm = TRUE)), psr = replicate(runs, mean(psr(as.matrix(com.simulator(comminmax[1], comminmax[2], nrow(samp), Species)), tree, compute.var = FALSE)[, 1], na.rm = TRUE)), pse = replicate(runs, mean(pse(as.matrix(com.simulator(comminmax[1], comminmax[2], nrow(samp), Species)), tree)[, 1], na.rm = TRUE)), psc = replicate(runs, mean(psc(as.matrix(com.simulator(comminmax[1], comminmax[2], nrow(samp), Species)), tree)[, 1], na.rm = TRUE)))
    					quantiles.null <- quantile(nulls, probs = c(alpha/2, 1 - (alpha/2)))
    					mean.null <- mean(nulls)
    					mean.obs <- switch(metric, psv = mean(psv(samp, tree, compute.var = FALSE)[, 1], na.rm = TRUE), psr = mean(psr(samp, tree, compute.var = FALSE)[, 1], na.rm = TRUE), pse = mean(pse(samp, tree)[, 1], na.rm = TRUE), psc = mean(psc(samp, tree)[, 1], na.rm = TRUE))
        									if (mean.obs <= quantiles.null[1]) {
            								phylo.structure = "underdispersed"
        									}
        									else {
            									if (mean.obs >= quantiles.null[2]) {
                									phylo.structure = "overdispersed"
            									}
            								else {
                								phylo.structure = "random"
            									}
        									}
        									if ((null.model != "independentswap") && (null.model != "trialswap")) {
            									it = NA
        									}
        									return(list(metric = metric, null.model = null.model, runs = runs, it = it, mean.obs = mean.obs, mean.null = mean.null, quantiles.null = quantiles.null, phylo.structure = phylo.structure, null.means = nulls))
    					}						
    									
    		}
}
    		
phylostruct.rev(comm, brach, c(2,5), env=NULL, metric="psv", null.model="other", runs=1000, alpha= 0.05)
phylostruct.rev(comm, brach, env=NULL, metric="psv", null.model="richness", runs=100, alpha= 0.05)