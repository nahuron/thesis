#random bits of code from Brachymeles CA


for (a in 1:length(results.gen)){
  #code all values of MPD
  results.hold <- sign(results.gen[[a]]$ses.mpd)
  #code all values of MNTD
  results.hold <- rbind(results.hold, sign(results.gen[[a]]$ses.mntd))
  #code all values of PSV
  psv.hold <- results.gen[[a]]$psv
  for(b in 1:length(psv.hold)){
    if(psv.hold[b]<0.854){
      psv.hold[b] <- -1
    }
    else{
      psv.hold[b] <- 1
    }}
  results.hold <- rbind(results.hold, psv.hold)
  print(a);print(abs(colSums(results.hold)));print(results.hold);print(rbind(results.gen[[a]]$ses.mpd,results.gen[[a]]$ses.mntd, results.gen[[a]]$psv))
  print(length(which(abs(colSums(results.hold))==1)))
}



par(mfrow=c(2,2))
#Mean Test Omission
hist(enms[enms$Method %in% c("Standard", "Jackknife"),2], col=rgb(1,0,0,0.5), xlim=c(0,1), xlab = "Test Omission", main="Mean Model")
hist(enms[enms$Method %in% c("PS Standard", "PS Jackknife"),2], add=T, col=rgb(0,0,1,0.5))
#Best Model Test Omission
hist(enms[enms$Method %in% c("Standard", "Jackknife"),4], col=rgb(1,0,0,0.5), xlim=c(0,1), xlab = "Test Omission", main="Best Model")
hist(enms[enms$Method %in% c("PS Standard", "PS Jackknife"),4], add=T, col=rgb(0,0,1,0.5))
#Mean Test AUC
hist(enms[enms$Method %in% c("Standard", "Jackknife"),3], col=rgb(1,0,0,0.5), xlim=c(0,1), xlab = "Test AUC", main="Mean Model")
hist(enms[enms$Method %in% c("PS Standard", "PS Jackknife"),3], add=T, col=rgb(0,0,1,0.5))
#Best Model Test AUC
hist(enms[enms$Method %in% c("Standard", "Jackknife"),5], col=rgb(1,0,0,0.5), xlim=c(0,1), xlab = "Test AUC", main="Best Model")
hist(enms[enms$Method %in% c("PS Standard", "PS Jackknife"),5], add=T, col=rgb(0,0,1,0.5))