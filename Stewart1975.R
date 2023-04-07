rm(list=ls())

a18 <- function(temp){
  eq <- exp(1137/temp^2-0.4156/temp-0.00207)
  return(eq)
}
a17 <- function(temp){
  eq <- exp(1137/temp^2-0.4156/temp-0.00207)^0.529
  return(eq)
}
a2 <- function(temp){
  eq <- exp(24844/temp^2-76.248/temp+0.05261)
  return(eq)
}
diff16_18 <- 1.0285
diff16_17 <- 1.0285^0.518
diff1_2 <- 1.0251

get_D17O <- function(d18O,d17O){
  (log(d17O/1000+1)-0.528*log(d18O/1000+1))*10^6
}
get_d17O <- function(d18O,D17O){
  1000*(exp(D17O/10^6+0.528*log(d18O/1000+1))-1)
}

dxs_i <- 10
D17O_i <- 20
d18_i <- -8
d17_i <- get_d17O(d18_i,D17O_i)
d2_i <- dxs_i+8*d18_i

temp <- seq(0,30,0.1)+273.15
rh <- seq(0.4,1,0.001)
f <- c(0.95,0.8,0.65)
d18O_change <- dxs_change <- D17O_change <- array(NA,dim=c(length(temp),length(rh),3))

for (i in 1:length(temp)){
  for (j in 1:length(rh)){
    for (k in 1:length(f)){
      gamma <- a18(temp[i])*rh[j]/(1-a18(temp[i])*(diff16_18)^0.58*(1-rh[j]))
      beta <- (1-a18(temp[i])*diff16_18^0.58*(1-rh[j]))/(a18(temp[i])*diff16_18^0.58*(1-rh[j]))
      d18_f <- (f[k]^beta-1)*(d18_i+1000)*(1-gamma/a18(temp[i]))+d18_i

      gamma <- a17(temp[i])*rh[j]/(1-a17(temp[i])*(diff16_17)^0.58*(1-rh[j]))
      beta <- (1-a17(temp[i])*diff16_17^0.58*(1-rh[j]))/(a17(temp[i])*diff16_17^0.58*(1-rh[j]))
      d17_f <- (f[k]^beta-1)*(d17_i+1000)*(1-gamma/a17(temp[i]))+d17_i
      
      gamma <- a2(temp[i])*rh[j]/(1-a2(temp[i])*(diff1_2)^0.58*(1-rh[j]))
      beta <- (1-a2(temp[i])*diff1_2^0.58*(1-rh[j]))/(a2(temp[i])*diff1_2^0.58*(1-rh[j]))
      d2_f <- (f[k]^beta-1)*(d2_i+1000)*(1-gamma/a2(temp[i]))+d2_i
      
      d18O_change[i,j,k] <- d18_f-d18_i
      dxs_change[i,j,k] <- (d2_f-8*d18_f)-dxs_i
      D17O_change[i,j,k] <- get_D17O(d18_f,d17_f)-D17O_i
    }
  }
}

par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
# 1
plot(NA,NA,xlim=c(0,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp-273.15,rh,d18O_change[,,1],add=T,levels=c(0.6,1.2,1.3),
        col="blue",labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,d18O_change[,,2],add=T,levels=c(1,3,5),
        labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,d18O_change[,,3],add=T,levels=seq(2,8,2),
        col="brown1",labcex=0.8,vfont=NULL,method="edge")
text(-1,102.5/100,expression(paste("(a) ",Delta,"(",delta^18,"O)"," (\u2030)")),cex=0.9,pos=4)
legend("bottomleft",c("E = 0.05","E = 0.2","E = 0.35"),col=c("blue","black","brown1"),lwd=1,cex=0.9)

# 2
plot(NA,NA,xlim=c(0,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp-273.15,rh,dxs_change[,,1],add=T,levels=-c(1,3,5,6),
        col="blue",labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,dxs_change[,,2],add=T,levels=-c(5,10,15,23),
        labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,dxs_change[,,3],add=T,levels=-seq(10,40,10),
        col="brown1",labcex=0.8,vfont=NULL,method="edge")
text(-1,102.2/100,expression(paste("(b) ",Delta,"(d-excess) (","\u2030",")")),cex=0.9,pos=4)

# 3
plot(NA,NA,xlim=c(0,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp-273.15,rh,D17O_change[,,1],add=T,levels=-c(3,4),
        col="blue",labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,D17O_change[,,2],add=T,levels=-c(8,16,20),
        labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,D17O_change[,,3],add=T,levels=-c(10,25,35,39),
        col="brown1",labcex=0.8,vfont=NULL,method="edge")
text(-1,102.5/100,expression(paste("(c) ",Delta,"(",Delta*minute^17,"O)"," (per meg)")),cex=0.9,pos=4)

# 4
plot(NA,NA,xlim=c(0,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp-273.15,rh,D17O_change[,,1]/dxs_change[,,1],add=T,levels=c(0.6,0.8,1,2),
        col="blue",labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,D17O_change[,,2]/dxs_change[,,2],add=T,levels=c(0.8,1,1.5,2),
        labcex=0.8,vfont=NULL)
contour(temp-273.15,rh,D17O_change[,,3]/dxs_change[,,3],add=T,levels=c(0.8,1,1.5,2),
        col="brown1",labcex=0.8,vfont=NULL,method="edge")
text(-1,102.5/100,expression(paste("(d) ",Delta,"(",Delta*minute^17,"O)/",
                                 Delta,"(d-excess)")),cex=0.9,pos=4)