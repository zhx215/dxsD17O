rm(list=ls())
load("C:\\Users\\zhyxi\\Desktop\\17O project\\era_resampling_monthly_v2.RData")
load("C:\\Users\\zhyxi\\Desktop\\17O project\\lmdz4_resampling.RData")

par(mfrow=c(3,2),mar=c(3,3,2,2),xpd=T)
# 1
plot(NA,NA,xlim=c(12,110),ylim=c(0,10),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(20,110,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,10,2),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("oceanic RH (%)",1,padj=1.8,cex=0.8)
mtext("Probability density (%)",2,padj=-1.8,cex=0.8)
lines(rhs_unit,colSums(rhs_freq_slice)/sum(colSums(rhs_freq_slice))*100)
legend("topleft","(a)",x.intersp=0.05,bty = "n",cex=1.2)

# 2
plot(NA,NA,xlim=c(12,110),ylim=c(-27,28),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(20,110,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,20,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("oceanic RH (%)",1,padj=1.8,cex=0.8)
mtext(expression(paste("Oceanic T"[d]," (\u00B0C)")),2,padj=-1.2,cex=0.8)
for (i in 1:length(rhs_unit)){
  if (sum(colSums(td_freq_per_rhs_slice[,i,]))>0){
    smp <- sample(d2m_ocean_unit,1e5,replace=T,prob=colSums(td_freq_per_rhs_slice[,i,]))
    result <- quantile(smp,c(0.025,0.16,0.5,0.84,0.975))
    lines(rep(rhs_unit[i],2),c(result[1],result[5]),lwd=1,col="gray")
    lines(rep(rhs_unit[i],2),c(result[2],result[4]),lwd=1)
    points(rhs_unit[i],result[3],pch=21,bg="white",cex=0.7)
  }
}
legend("bottomright",c("median","63%","95%"),lwd=c(NA,1,1),
       pch=c(21,NA,NA),pt.cex=c(0.7,NA,NA),col=c("black","black","gray"),
       bty="n")
legend("topleft","(b)",x.intersp=0.05,bty = "n",cex=1.2)

# 3
plot(NA,NA,xlim=c(-27,28),ylim=c(0,14),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-25,25,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,14,2),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Oceanic T"[d]," (\u00B0C)")),1,padj=1.6,cex=0.8)
mtext("Probability density (%)",2,padj=-1.8,cex=0.8)
lines(d2m_ocean_unit,
      colSums(td_freq_per_rhs_slice[,39,])/sum(colSums(td_freq_per_rhs_slice[,39,]))*100,col="deeppink")
lines(d2m_ocean_unit,
      colSums(td_freq_per_rhs_slice[,59,])/sum(colSums(td_freq_per_rhs_slice[,59,]))*100)
lines(d2m_ocean_unit,
      colSums(td_freq_per_rhs_slice[,79,])/sum(colSums(td_freq_per_rhs_slice[,79,]))*100,col="deepskyblue")
box(lwd=1)
legend(-26.9,12,c("oceanic RH = 50%","oceanic RH = 70%","oceanic RH = 90%"),
       lwd=c(1,1,1),col=c("deeppink","black","deepskyblue"),box.lwd=NA)
legend("topleft","(c)",x.intersp=0.05,bty = "n",cex=1.2)

# 4
plot(NA,NA,xlim=c(-27,28),ylim=c(-40,-9),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-25,25,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,-10,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Oceanic T"[d]," (\u00B0C)")),1,padj=1.6,cex=0.8)
mtext(expression(paste("Oceanic vapor ",delta^18,"O (","\u2030",")")),2,padj=-0.9,cex=0.8)
for (i in 20:length(td_lmdz_unit)){
  smp <- sample(d18O_lmdz_unit,1e5,replace=T,prob=d18O_freq_per_td[i,])
  result <- quantile(smp,c(0.025,0.16,0.5,0.84,0.975))
  lines(rep(td_lmdz_unit[i],2),c(result[1],result[5]),lwd=1,col="gray")
  lines(rep(td_lmdz_unit[i],2),c(result[2],result[4]),lwd=1)
  points(td_lmdz_unit[i],result[3],pch=21,bg="white",cex=0.7)
}
legend("topleft","(d)",x.intersp=0.05,bty = "n",cex=1.2)

# 5
plot(NA,NA,xlim=c(-67,27),ylim=c(0,14),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-60,20,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,14,2),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Continental T"[d]," (\u00B0C)")),1,padj=1.6,cex=0.8)
mtext("Probability density (%)",2,padj=-1.8,cex=0.8)
lines(d2m_land_unit,colSums(d2m_land_freq_slice)/sum(colSums(d2m_land_freq_slice))*100)
legend("topleft","(e)",x.intersp=0.05,bty = "n",cex=1.2)

# 6
get_rh <- function(td,t){
  return(exp(17.67*td/(td+243.5)-17.67*t/(t+243.5))*100)
}

plot(NA,NA,xlim=c(0,27),ylim=c(0,40),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,25,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,40,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Continental T"[d]," (\u00B0C)")),1,padj=1.6,cex=0.8)
mtext("Surface temperature (%)",2,padj=-1.7,cex=0.8)

t40 <- t60 <- t80 <- vector()
t <- seq(0,50,0.01)
td <- seq(0,27,0.01)
for (i in 1:length(td)){
  i40 <- which.min(abs(get_rh(td[i],t)-40))
  i60 <- which.min(abs(get_rh(td[i],t)-60))
  i80 <- which.min(abs(get_rh(td[i],t)-80))
  t40[i] <- t[i40]
  t60[i] <- t[i60]
  t80[i] <- t[i80]
}
lines(td,t80,lty=5,lwd=1,col="aquamarine3")
lines(td,t60,lty=5,lwd=1)
lines(td[1:2387],t40[1:2387],lty=5,lwd=1,col="coral3")

for (i in 68:length(d2m_land_unit)){
  smp <- sample(t2m_land_unit,1e5,replace=T,prob=colSums(t2m_freq_per_d2m_land_slice[,i,]))
  result <- quantile(smp,c(0.025,0.16,0.5,0.84,0.975))
  lines(rep(d2m_land_unit[i],2),c(result[1],result[5]),lwd=1,col="gray")
  lines(rep(d2m_land_unit[i],2),c(result[2],result[4]),lwd=1)
  points(d2m_land_unit[i],result[3],pch=21,bg="white",cex=0.7)
}

legend("bottomright",c("surface RH = 40%","surface RH = 60%","surface RH = 80%"),
       lwd=1,lty=5,col=c("coral3","black","aquamarine3"),bty="n")
legend("topleft","(f)",x.intersp=0.05,bty = "n",cex=1.2)
