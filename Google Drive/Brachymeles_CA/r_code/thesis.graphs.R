#graphing the results

setwd("/Users/nicholashuron/Dropbox/Huron, Nick Masters/")

dev.new()
phylograph <- read.csv("phylogenetic.comparisons.csv", header=1)

jpeg(filename="Figure1R.jpg", height= 600, width=600, units="px", quality=200)

table(phylograph$dispersion,phylograph$metric)->metrics.tograph
names <- c("MNTD\n0.055822478", "MPD\n0.094342662", "PSV\n0.86684")
par(mar = c(5,5,2,2) + 0.1)
barplot(metrics.tograph,beside=T, col=c(0,1,16), ylim=c(0,30), names.arg=names, cex.axis=1.2, cex.lab=1.75, xlab="(Mean Values Below Metrics)", ylab="Observations")

legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16), cex=2)


dev.off()

#luzon
jpeg(filename="Figure2Rluz.jpg", height= 600, width=600, units="px", quality=200)

phylograph.luzon <- phylograph[phylograph[,3]=="luzon",]
luzon.tograph <- table(phylograph.luzon$dispersion, phylograph.luzon$metric)
luzon.mean <- sapply(split(as.numeric(phylograph.luzon$statistic), as.factor(phylograph.luzon$metric)), mean)
luzon.names <- paste(names(luzon.mean), as.character(round(luzon.mean, digits=8)), sep="\n")
par(mar = c(5,5,2,2) + 0.1)
barplot(luzon.tograph,beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=luzon.names, cex.axis=1.2, cex.lab=1.75, xlab="(Mean Values Below Metrics)", ylab="Observations")
legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16), cex=2)
dev.off()


#mindanao
jpeg(filename="Figure2Rmina.jpg", height= 600, width=600, units="px", quality=200)

phylograph.mina <- phylograph[phylograph[,3]=="mindanao",]
mina.tograph <- table(phylograph.mina$dispersion, phylograph.mina$metric)
mina.mean <- sapply(split(as.numeric(phylograph.mina$statistic), as.factor(phylograph.mina$metric)), mean)
mina.names <- paste(names(mina.mean), as.character(round(mina.mean, digits=8)), sep="\n")
par(mar = c(5,5,2,2) + 0.1)
barplot(mina.tograph,beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=mina.names, cex.axis=1.2, cex.lab=1.75, xlab="(Mean Values Below Metrics)", ylab="Observations")
legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16), cex=2)
dev.off()

#mindoro

jpeg(filename="Figure2Rmino.jpg", height= 600, width=600, units="px", quality=200)

phylograph.mino <- phylograph[phylograph[,3]=="mindoro",]
mino.tograph <- table(phylograph.mino$dispersion, phylograph.mino$metric)
mino.mean <- sapply(split(as.numeric(phylograph.mino$statistic), as.factor(phylograph.mino$metric)), mean)
mino.names <- paste(names(mino.mean), as.character(round(mino.mean, digits=8)), sep="\n")
barplot(mino.tograph,beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=mino.names, main="Mindoro FR Phylogenetic Dispersion")
legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16))

dev.off()

#lubang

jpeg(filename="Figure2Rlub.jpg", height= 600, width=600, units="px", quality=200)

phylograph.lub <- phylograph[phylograph[,3]=="lubang",]
lub.tograph <- table(phylograph.lub$dispersion, phylograph.lub$metric)
lub.mean <- sapply(split(as.numeric(phylograph.lub$statistic), as.factor(phylograph.lub$metric)), mean)
lub.names <- paste(names(lub.mean), as.character(round(lub.mean, digits=8)), sep="\n")
barplot(lub.tograph,beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=lub.names, main="Lubang Phylogenetic Dispersion")
legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16))

dev.off()

#visayan

jpeg(filename="Figure2Rvis.jpg", height= 600, width=600, units="px", quality=200)

phylograph.vis <- phylograph[phylograph[,3]=="visayan",]
vis.tograph <- table(phylograph.vis$dispersion, phylograph.vis$metric)
vis.mean <- sapply(split(as.numeric(phylograph.vis$statistic), as.factor(phylograph.vis$metric)), mean)
vis.names <- paste(names(vis.mean), as.character(round(vis.mean, digits=8)), sep="\n")
barplot(vis.tograph,beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=vis.names, main="Visayan FR Phylogenetic Dispersion")
legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16))

dev.off()

#Romblon IG

jpeg(filename="Figure2Rrom.jpg", height= 600, width=600, units="px", quality=200)


phylograph.rom <- phylograph[phylograph[,3]=="romblon",]
rom.tograph <- table(phylograph.rom$dispersion, phylograph.rom$metric)
rom.mean <- sapply(split(as.numeric(phylograph.rom$statistic), as.factor(phylograph.rom$metric)), mean)
rom.names <- paste(names(rom.mean), as.character(round(rom.mean, digits=8)), sep="\n")
barplot(rom.tograph,beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=rom.names, main="Romblon IG Phylogenetic Dispersion")
legend("topright", legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16))

dev.off()

#Plot them all in one figure

phylocom.byfr <- list(luzon.tograph, mina.tograph, mino.tograph, vis.tograph, rom.tograph, lub.tograph)
names(phylocom.byfr) <- c("Luzon Coms", "Mindanao Coms", "Mindoro Coms", "Visayan Coms", "Romblon Coms", "Lubang Coms")
phylocom.byfr.means <- c(luzon.mean, mina.mean, mino.mean, vis.mean, rom.mean, lub.mean)
phylocom.byfr.names <- list(luzon.names, mina.names, mino.names, vis.names, rom.names, lub.names)

dev.new()

jpeg(filename="Figure2RALL.jpg", height= 600, width=800, units="px", quality=300)
par(mfrow=c(2,3)) 
for (i in 1:(length(phylocom.byfr)))	{
	par(mar = c(5,5,3,2) + 0.1)
	barplot(phylocom.byfr[[i]],cex.axis=1.3, cex.lab=1.5, beside=T, col=c(0,1,16), ylim=c(0,15), names.arg=phylocom.byfr.names[[i]], main=names(phylocom.byfr)[[i]], xlab="(Mean Values Below Metrics)", ylab="Observations", cex.main=2)
	#abline(v = mean(phylocom.byfr.means[[i]]), col = "blue", lwd = 2, lty=2)
											}
legend("topright", cex=1.4, legend=c("Overdispersed","Underdispersed", "Random"), fill=c(0,1,16))
dev.off()

#hellinger's i

dev.new()

jpeg(filename="Figure3R.jpg", height= 600, width=600, units="px", quality=200)

brach.i <- read.csv("hellingersI.csv", header=1)
library(moments)
skewness(brach.i[,3]) -> brachi.skew
brachi.skew <- paste("Skewness", round(brachi.skew, digits=5), sep=": ") 
mline<-mean(brach.i[,3])

par(mar = c(6,5,5,5) + 0.1)
hist(brach.i[,3], breaks=20, cex.axis=1.2, cex.lab=1.5, xlim=c(0,1.0), ylim=c(0,40), xaxt='n', main="Niche Overlap as Measured by\nHellinger's I", cex.main=2, xlab=c("I Statistic",as.character(brachi.skew)))
axis(side=1, cex.axis=1.2, at=seq(0, 1, 0.1))
abline(v = mline, col = "red", lwd = 2, lty=2)
text(mline+0.075, 35, cex=1.4, labels="Mean", col="red")

dev.off()

noverlap <- read.csv("hellingersI.csv", header=1, na.strings="")


#group by Same or Diferent Faunal Region
SFR <- noverlap[noverlap[,4]==1,]
SFR.mean <- mean(SFR[,3])

DFR <- noverlap[noverlap[,4]==0,]
DFR.mean <- mean(DFR[,3])


#group by Same or Different Community
SC <- SFR[SFR[,5]=="y",]
SC.mean <- mean(SC[,3])

DC <- noverlap[noverlap[,5]=="n",]
DC.mean <- mean(DC[,3])

#partition DC by same FR or different
DCSFR <- DC[DC[,4]==1,]
DCSFR.mean <- mean(DCSFR[,3])

DCDFR <- DC[DC[,4]==0,]
DCDFR.mean <- mean(DCDFR[,3])

#plot the mean I by FR for Communities

FRI <- SC$Faunal.Region.ID

i.luz <- subset(SC, FRI=='luzon')
i.mina <-subset(SC, FRI=='mindanao')
i.mino <-subset(SC, FRI=='mindoro')
i.vis <- subset(SC, FRI=='visayan')
i.rom <- subset(SC, FRI=='romblon')

i.luz.mean <- mean(i.luz[,3])
i.mina.mean <- mean(i.mina[,3])
i.mino.mean <- mean(i.mino[,3])
i.vis.mean <- mean(i.vis[,3])
i.rom.mean <- mean(i.rom[,3])

i.luz.mean
i.mina.mean
i.mino.mean
i.vis.mean
i.rom.mean

niche.overlap.means <- c(SFR.mean, DFR.mean, SC.mean, DC.mean, DCSFR.mean, DCDFR.mean, i.luz.mean, i.mina.mean, i.mino.mean, i.vis.mean, i.rom.mean)
names(niche.overlap.means) <- c("Same FR", "Different FR", "Same Com", "Different Com", "SFR DC", "DFR DC", "Luzon Coms", "Mindanao Coms", "Mindoro Coms", "Visayan Coms", "Romblon Coms")
niche.overlap.dists <- list(SFR[,3], DFR[,3], SC[,3], DC[,3], DCSFR[,3], DCDFR[,3], i.luz[,3], i.mina[,3], i.mino[,3], i.vis[,3], i.rom[,3])
names(niche.overlap.dists) <- c("Same FR", "Different FR", "Same Com", "Different Com", "SFR DC", "DFR DC", "Luzon Coms", "Mindanao Coms", "Mindoro Coms", "Visayan Coms", "Romblon Coms")

dev.new()
par(mfrow=c(4,3)) 
for (i in 1:(length(niche.overlap.dists)))	{
	#dev.new()
	hist(niche.overlap.dists[[i]], main=names(niche.overlap.dists)[[i]], cex.axis=0.8, xlim=c(-0.1,1.1), ylim=c(0,50))
	abline(v = mean(niche.overlap.means[i]), col = "blue", lwd = 2, lty=2)
											}

#split up the graphs by topics

noverlap.dists.large <- list(SFR[,3], DFR[,3], SC[,3], DC[,3], DCSFR[,3], DCDFR[,3])
names(noverlap.dists.large) <- c("Same FR", "Different FR", "Same Com", "Different Com", "SFR DC", "DFR DC")
noverlap.large.means <- c(SFR.mean, DFR.mean, SC.mean, DC.mean, DCSFR.mean, DCDFR.mean)

dev.new()

jpeg(filename="Figure3ALLCOMPARE.jpg", height= 1000, width=1400, units="px", quality=200)

par(mfrow=c(2,3)) 
for (i in 1:(length(noverlap.dists.large)))	{
	par(mar = c(6,5,5,5) + 0.1)
	hist(noverlap.dists.large[[i]], xaxt='n', xlab="I Statistic", main=names(noverlap.dists.large)[[i]], cex.main=2, cex.lab=1.3, cex.axis=1.3, xlim=c(0,1), ylim=c(0,50))
	axis(side=1, cex.axis=1.3, at=seq(0, 1, 0.1))
	abline(v = noverlap.large.means[i], col = "red", lwd = 2, lty=2)
	text(noverlap.large.means[i]+0.075, 45, cex=1.4, labels="Mean", col="red")
											}
											
dev.off()



noverlap.dists.byfr <- list(i.luz[,3], i.mina[,3], i.mino[,3], i.vis[,3], i.rom[,3])
names(noverlap.dists.byfr) <- c("Luzon Coms", "Mindanao Coms", "Mindoro Coms", "Visayan Coms", "Romblon Coms")
noverlap.byfr.means <- c(i.luz.mean, i.mina.mean, i.mino.mean, i.vis.mean, i.rom.mean)

dev.new()

jpeg(filename="Figure3ALLFR.jpg", height= 800, width=1000, units="px", quality=200)

par(mfrow=c(2,3)) 
for (i in 1:(length(noverlap.dists.byfr)))	{
	par(mar = c(6,5,5,5) + 0.1)
	hist(noverlap.dists.byfr[[i]], xaxt='n', xlab="I Statistic", main=names(noverlap.dists.byfr)[[i]], cex.main=1.5, cex.axis=0.8, cex.lab=1.25, xlim=c(0,1), ylim=c(0,6))
	axis(side=1, cex.axis=1.2, at=seq(0, 1, 0.1))
	abline(v = noverlap.byfr.means[i], col = "red", lwd = 2, lty=2)
	text(noverlap.byfr.means[i]+0.1, 5.5, cex=1.4, labels="Mean", col="red")
											}

dev.off()

#luzon only
dev.new()

jpeg(filename="Figure3luz.jpg", height= 600, width=600, units="px", quality=200)

par(mar = c(6,5,5,5) + 0.1)
hist(i.luz[,3], cex.axis=1.2, cex.lab=1.5, xlim=c(0,1.0), ylim=c(0,6), xaxt='n', main="Niche Overlap as Measured by\nHellinger's I", cex.main=2, xlab="I Statistic")
axis(side=1, cex.axis=1.2, at=seq(0, 1, 0.1))
abline(v = i.luz.mean, col = "red", lwd = 2, lty=2)
text(i.luz.mean+0.075, 5.5, cex=1.4, labels="Mean", col="red")

dev.off()

#mindanao only
dev.new()

jpeg(filename="Figure3mina.jpg", height= 600, width=600, units="px", quality=200)

par(mar = c(6,5,5,5) + 0.1)
hist(i.mina[,3], cex.axis=1.2, cex.lab=1.5, xlim=c(0,1.0), ylim=c(0,6), xaxt='n', main="Niche Overlap as Measured by\nHellinger's I", cex.main=2, xlab="I Statistic")
axis(side=1, cex.axis=1.2, at=seq(0, 1, 0.1))
abline(v = i.mina.mean, col = "red", lwd = 2, lty=2)
text(i.mina.mean+0.075, 5.5, cex=1.4, labels="Mean", col="red")

dev.off()

#single histogram
dev.new()

#FR Comparison
hist(SFR[,3], ylim=c(0,50), xlim=c(-0.1,1.1),col=rgb(0,0,0,.5))
hist(DFR[,3], add=T, col=rgb(1,1,1,0.5))
text(mean(SFR[,3])+0.1, 40, labels="Same FR\nMean", col="blue")
abline(v = mean(SFR[,3]), col = "blue", lwd = 2, lty=2)
abline(v=mean(DFR[,3]),lty=2, lwd=2, col="red")
text(mean(DFR[,3])-0.12, 48, labels="Different FR\nMean", col="red")
legend("topright", legend=c("Same Faunal Region","Different Faunal Region"), fill=c(rgb(0,0,0,.5),rgb(1,1,1,0.5)))
legend("topright", legend=c("Same Faunal Region","Different Faunal Region"), fill=c(rgb(1,1,1,.5),rgb(1,1,1,0.5)))
