library (ape)
library (geiger)

#setwd("/Volumes/NAHURONJUMP/brachymeles.genetic/")		#work in jump drive folders

#setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/")	#work in dropbox folders

setwd("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/Datasets/Brachymeles Chronogram/")

#brachtree <- read.tree("Newick.tre")

brachtree <- read.nexus("Brachymeles_BEAST_Combtree.tre")

#brachtree <- read.nexus("brachnexus.nex")

#see what taxa are present in the tree
brachtree$tip.label

#visualize the tree
plot.phylo(brachtree, no.margin=T)

#branches to drop due to no need for outgroup
branches.outgroups <- c("'Lygosoma.bowringii.LSUHC6998'", "'Lygosoma.bowringii.LSUHC6970'", "'Lygosoma.quadrupes.LSUHC8002'", "'Lygosoma.quadrupes.LSUHC8403'", "'Lygosoma.LSUHC9321'", "'Tachydromus.sexilineatus.KUFS217'")

#drop outgroups and put resulting phylogeny in a new object
culled.brachtree <- drop.tip(brachtree, branches.outgroups)

#Revised branches dropped because of multiple individuals per species
branches.dropped <- c("'apus.Kalimantan.RMBR2040'", "'apus.Sabah.SP06915'", "'miriamae.Thailand.DSM1293'", "'miriamae.Thailand.DSM1363'", "'talinis.Negros.RMB3305'","'kadwa.Luzon.CDS4057'","'tungaoi.Masbate.CDS5123'", "'boulengeri.Polillo.CDS1118'","'taylori.Cebu.CDS3602'", "'mindorensis.Mindoro.ELR0435'", "'muntingkamay.Palali.ELR1278'","'elerae.Kalinga.CDS3481'","'boholensis.Bohol.CDS4473'","'makusog.Catanduanes.CDS2369'","'schadenbergi.Zambnga.RMB10525'", "'orientalis.Bohol.CDS4692'", "'bicolor.Luzon.RMB11986'", "'tiboliorum.SMindanao.ACD5747'", "'tiboliorum.SMindanao.ACD5727'", "'tiboliorum.Mindanao.H1641'","'samad.Leyte.CDS3390'", "'samad.Leyte.CDS3138'", "'samad.Leyte.CDS3143'", "'samad.Leyte.CDS3400'", "'samad.Leyte.CDS3151'", "'samad.Samar.CDS2552'", "'samad.Samar.CDS2543'", "'samad.Samar.CDS2762'", "'samad.Samar.CDS2706'","'hilong.MtHilong.CDS5431'", "'hilong.MtHilong.CDS5432'","'pathfinderi.Glan.CDS5192'", "'pathfinderi.Glan.CDS5193'", "'pathfinderi.Glan.CDS5199'", "'pathfinderi.Glan.CDS5242'", "'pathfinderi.Glan.CDS5241'", "'pathfinderi.Glan.CDS5201'", "'pathfinderi.Glan.CDS5198'","'gracilis.Toril.H1175'", "'gracilis.Mindanao.Davao.H1176'","'gracilis.Baracatan.H1169'", "'gracilis.Cotabato.ACD5262'", "'gracilis.Mindanao.H1428'", "'gracilis.Mindanao.H1332'", "'gracilis.Mindanao.H1429'", "'cebuensis.Cebu.CDS3597'", "'cebuensis.Cebu.CDS3599'","'paeforum.Leyte.CDS3417'", "'paeforum.Leyte.CDS3419'", "'paeforum.Leyte.CDS3416'","'libayani.Lapinig.CDS3645'", "'libayani.Lapinig.CDS3677'", "'libayani.Lapinig.CDS3679'", "'libayani.Lapinig.CDS3692'", "'libayani.Lapinig.CDS3700'", "'libayani.Lapinig.CDS3669'", "'samarensis.Samar.CDS3034'", "'samarensis.Samar.CDS3036'", "'minimus.Catanduanes.CDS2371'", "'minimus.Catanduanes.CDS2372'","'lukbani.Labo.RMB9681'", "'lukbani.Labo.RMB9629'", "'lukbani.Labo.RMB9661'", "'lukbani.Labo.RMB9909'", "'lukbani.Labo.RMB9908'","'brevidactylus.Bicol.RMB4039'","'cobos.Catanduanes.CDS5248'", "'cobos.Catanduanes.CDS5254'","'bicolandia.Luzon.CDS4050'", "'bicolandia.Luzon.CDS4051'", "'bicolandia.Luzon.CDS4089'", "'isangdaliri.Luzon.RMB12650'","'mapalanggaon.Masbate.CDS5101'", "'tridactylus.Sibalom.GVAG287'", "'tridactylus.Negros.CDS3619'", "'tridactylus.Negros.CDS3618'", "'tridactylus.Negros.CDS3617'", "'tridactylus.Alojipan.CDS1647'", "'tridactylus.Pandan.CDS1529'", "'tridactylus.Pandan.CDS1528'","'bonitae.Luzon.ACD6099'","'bonitae.Luzon.RMB3680'","'bonitae.Luzon.RMB3681'","'bonitae.Polillo.RMB8902'","'bonitae.Polillo.RMB8901'","'burksi.Marinduque.CDS3595'","'burksi.Marinduque.CDS3596'", "'burksi.Mindoro.RMB4896'","'ligtas.Lubang.CDS3881'", "'ligtas.Lubang.CDS3928'","'species2.PalayPalay.ACD1236'", "'ilocandia.NLuzon.ACD6175'", "'ilocandia.NLuzon.ACD2002'","'ilocandia.CamigNort.RMB7324'", "'ilocandia.CamigNort.RMB7286'","'ilocandia.Kalinga.CDS3483'","'ilocandia.Gonzaga.CDS3731'","'ilocandia.NLuzon.RMB11924'","'hilong.MtHilong.CDS5429'", "'species1.Zambales.RMB4624'")

#drop branches and put resulting phylogeny in a new object
culled.brachtree <- drop.tip(culled.brachtree, branches.dropped)

#branches to drop due to multiple specimens per species
branches.dropped <- c("'sp.bon3.PalayPalay.ACD1368'", "'sp.bon1.Zambales.RMB4624'", "'miriamae.Thailand.DSM1293'", "'miriamae.Thailand.DSM1363'", "'talinis.Negros.RMB3305'", "'kadwa.Luzon.CDS4057'", "'tungaoi.Masbate.CDS5123'", "'boulengeri.Polillo.CDS1118'", "'taylori.Cebu.CDS3602'", "'mindorensis.Mindoro.ELR0435'", "'boholensis.Bohol.CDS4473'", "'makusog.Catanduanes.CDS2369'", "'schadenbergi.Zambnga.RMB10525'", "'orientalis.Bohol.CDS4692'", "'muntingkamay.Palali.ELR1278'", "'elerae.Kalinga.CDS3481'", "'tiboliorum.SMindanao.ACD5747'", "'tiboliorum.SMindanao.ACD5727'", "'tiboliorum.Mindanao.H1641'", "'samad.Leyte.CDS3390'", "'samad.Leyte.CDS3138'", "'samad.Leyte.CDS3143'", "'samad.Leyte.CDS3400'", "'samad.Leyte.CDS3151'", "'samad.Samar.CDS2552'", "'samad.Samar.CDS2543'", "'samad.Samar.CDS2762'", "'samad.Samar.CDS2706'", "'pathfinderi.Glan.CDS5192'", "'pathfinderi.Glan.CDS5193'", "'pathfinderi.Glan.CDS5199'", "'pathfinderi.Glan.CDS5242'", "'pathfinderi.Glan.CDS5241'", "'pathfinderi.Glan.CDS5201'", "'pathfinderi.Glan.CDS5198'", "'gracilis.Toril.H1175'", "'gracilis.Mindanao.Davao.H1176'", "'gracilis.Baracatan.H1169'", "'gracilis.Cotabato.ACD5262'", "'gracilis.Mindanao.H1428'", "'gracilis.Mindanao.H1332'", "'gracilis.Mindanao.H1429'", "'bicolor.Luzon.RMB11986'", "'bonitae.Luzon.RMB3681'", "'bonitae.Luzon.ACD6099'", "'bonitae.Polillo.RMB8902'", "'burksi.Marinduque.CDS3595'", "'burksi.Marinduque.CDS3596'", "'burksi.Mindoro.RMB4896'", "'mapalanggaon.Masbate.CDS5101'", "'tridactylus.Sibalom.GVAG287'", "'tridactylus.Negros.CDS3619'", "'tridactylus.Negros.CDS3618'", "'tridactylus.Negros.CDS3617'", "'tridactylus.Alojipan.CDS1647'", "'tridactylus.Pandan.CDS1529'", "'tridactylus.Pandan.CDS1528'", "'isangdaliri.Aurora.RMB12650'", "'samarensis.Samar.CDS3034'", "'minimus.Catanduanes.CDS2371'", "'minimus.Catanduanes.CDS2372'", "'samarensis.Samar.CDS3036'", "'lukbani.Labo.RMB9681'", "'lukbani.Labo.RMB9629'", "'lukbani.Labo.RMB9661'", "'lukbani.Labo.RMB9909'", "'lukbani.Labo.RMB9908'", "'brevidactylus.Bicol.RMB4039'", "'cobos.Catanduanes.CDS5248'", "'cobos.Catanduanes.CDS5254'", "'bicolandia.Luzon.CDS4050'", "'bicolandia.Luzon.CDS4051'", "'bicolandia.Luzon.CDS4089'", "'cebuensis.Cebu.CDS3597'", "'cebuensis.Cebu.CDS3599'", "'libayani.Lapinig.CDS3645'", "'libayani.Lapinig.CDS3677'", "'libayani.Lapinig.CDS3679'", "'libayani.Lapinig.CDS3692'", "'libayani.Lapinig.CDS3700'", "'libayani.Lapinig.CDS3669'", "'paeforum.Leyte.CDS3417'", "'paeforum.Leyte.CDS3419'", "'paeforum.Leyte.CDS3416'", "'sp.bon2.Lubang.CDS3881'", "'sp.bon2.Lubang.CDS3928'", "'sp.bon4.Kalinga.CDS3484'", "'sp.bon4.Gonzaga.CDS3730'", "'sp.bon4.Gonzaga.CDS3731'", "'sp.bon4.NLuzon.RMB11924'", "'sp.bon4.NLuzon.ACD6175'", "'sp.bon4.NLuzon.ACD2002'", "'sp.bon4.CamigNort.RMB7324'", "'sp.bon4.CamigNort.RMB7286'", "'sp.bon3.PalayPalay.ACD1236'", "'bonitae.Polillo.RMB8901'", "'sp.apus1.Kalimantan.RMBR2040'", "'apus.Sabah.SP06915'", "'hilong.MtHilong.CDS5431'", "'hilong.MtHilong.CDS5432'", "'hilong.MtHilong.CDS5429'","'ilocandia.CamigNort.RMB7286'", "'ilocandia.CamigNort.RMB7324'", "'ilocandia.Kalinga.CDS3483'", "'ilocandia.Kalinga.CDS3484'", "'ilocandia.NLuzon.ACD2002'", "'ilocandia.NLuzon.ACD6175'", "'ilocandia.NLuzon.RMB11924'", "'isangdaliri.Luzon.RMB12650'", "'ligtas.Lubang.CDS3928'", "'species2.PalayPalay.ACD1368'", "species1.Zambales.RMB4624")

#drop branches and put resulting phylogeny in a new object
culled.brachtree <- drop.tip(culled.brachtree, branches.dropped)

#visualize the new tree
plot.phylo(culled.brachtree, no.margin=T)

culled.brachtree$tip.label

#save the tree to file
write.tree(culled.brachtree, file="culledbrachtree_final.tre")
write.nexus(culled.brachtree, file="culledbrachtree_final.nex")


#analyze data
library(picante)

#read community table in (NEEDS TO BE REVISED IN NUMBERS and ALTERED)
comm <- read.table("/Users/nicholashuron/Dropbox/STUDENT FOLDERS/Huron, Nick/Huron_Nick_Masters/brachymeles.communitymatrix.txt", header=T, row.names=1)
commmat <- as.matrix(comm)

#read in tree file
brach <- read.nexus("culledbrachtree_final.nex")


#test statistics
#creates phylo distance matrix
phydist_brach <- cophenetic(brach)

brach$tip.label <- gsub("\\..*", "", brach$tip.label)

colnames(phydist_brach) <- gsub("\\..*", "", colnames(phydist_brach))
rownames(phydist_brach) <- gsub("\\..*", "", rownames(phydist_brach))

#null model info guide
#frequency == takes the same number of total instances for a species but varies them across communities (allowing for as many as 0 species in a community or as many as 6)
#richness == takes the same number of total species per community but varies the number of times species can exist in communities (can have a higher number occurrences across communities)
#taxa.labels == switches labels in the matrix (NOT PSV)
#sample.pool == randomize community data matrix by drawing species from those that occur in at least 1 community (NOT PSV)
#phylogeny.pool == randomize community data matrix by drawing species from a pool of a posible species in the phylogeny (mpd matrix) with equal probability (NOT PSV)

#standardize the effect sizes should use nulls either taxa labels, sample pool, phylogeny pool
#MPD == NRI*-1 (negative values indicate clustering) tree wide sensitivity
ses.mpd.taxa.labels <- ses.mpd(comm, phydist_brach, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

ses.mpd.frequency <- ses.mpd(comm, phydist_brach, null.model = "frequency",abundance.weighted = FALSE, runs = 1000)

ses.mpd.richness <- ses.mpd(comm, phydist_brach, null.model = "richness",abundance.weighted = FALSE, runs = 1000)

ses.mpd.sample.pool <- ses.mpd(comm, phydist_brach, null.model = "sample.pool",abundance.weighted = FALSE, runs = 1000)

ses.mpd.phylogeny.pool <- ses.mpd(comm, phydist_brach, null.model = "phylogeny.pool",abundance.weighted = FALSE, runs = 1000)



#MNTD == NTI sensitivity to clustering/evenness at the tips of a tree
ses.mntd.taxa.labels <- ses.mntd(comm, phydist_brach, null.model = "taxa.labels",abundance.weighted = FALSE, runs = 1000)

ses.mntd.frequency <- ses.mntd(comm, phydist_brach, null.model = "frequency",abundance.weighted = FALSE, runs = 1000)

ses.mntd.richness <- ses.mntd(comm, phydist_brach, null.model = "richness",abundance.weighted = FALSE, runs = 1000)

ses.mntd.sample.pool <- ses.mntd(comm, phydist_brach, null.model = "sample.pool",abundance.weighted = FALSE, runs = 1000)

ses.mntd.phylogeny.pool <- ses.mntd(comm, phydist_brach, null.model = "phylogeny.pool",abundance.weighted = FALSE, runs = 1000)



#PSV (closer to 0 indicates relatedness and clustering)
psv.result <- psv(comm, brach, compute.var=TRUE)

psv.result

#signed proportion of total deviation from mean PSV that each species contributes
psv.spp(comm, brach)


save.xlsx("community_results.xlsx", psv.result)

phylostruct(comm, brach, env=NULL, metric="psv", null.model="frequency", runs=1000, alpha= 0.05)

phylostruct(comm, brach, env=NULL, metric="psv", null.model="richness", runs=750, alpha=0.05)


hist(tester$null.means, xlim=c(0,1), ylim=c(0,300))
abline(v=psv.result$PSVs, col="blue"); text(x=(psv.result$PSVs+.2),y=310,labels= rownames(psv.result))


psv.result[psv.result[,1]>0.9,]





#create for loop to use as a means to obtain p values for individual community PSV's

#obtain null distribution
psv(randomizeMatrix(comm, null.model="frequency"), brach)



