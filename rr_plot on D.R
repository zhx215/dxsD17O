rm(list=ls())
setwd("C:\\Users\\zhyxi\\Desktop\\17O project")
source("raindrop re-evaporation model.R")

#
rh <- seq(0.4,1,0.01)
temp <- seq(0,30,1)
d <- c(0.8,1.6,2.4)

### run rr matrix ###
# rr_result_matrix <- array(NA,dim=c(length(temp),length(rh),length(d),5))
# count <- 0
# for (i in 1:length(temp)){
#   for (j in 1:length(rh)){
#     for (k in 1:length(d)){
#       rr_result_matrix[i,j,k,] <- reevap(101325,temp[i]+273.15,rh[j],d[k],0,0,0)
#       count <- count + 1
#       print(count/(length(temp)*length(rh)*length(d))*100)
#     }
#   }
# }
# print(Sys.time())

load("C:/Users/zhyxi/Desktop/17O project/rr_D2.RData")
# 17O converter
get_D17O <- function(d18O,d17O){
  (log(d17O/1000+1)-0.528*log(d18O/1000+1))*10^6
}

get_d17O <- function(d18O,D17O){
  1000*(exp(D17O/10^6+0.528*log(d18O/1000+1))-1)
}

#
d18O_i <- -8
dxs_i <- 10
D17O_i <- 20
d2H_i <- d18O_i*8+dxs_i
d17O_i <- get_d17O(d18O_i,D17O_i)
d18O_f <- (d18O_i+1000)*rr_result_matrix[,,,2]-1000
d17O_f <- (d17O_i+1000)*rr_result_matrix[,,,3]-1000
d2H_f <- (d2H_i+1000)*rr_result_matrix[,,,4]-1000
d18O_change <- d18O_f-d18O_i
dxs_change <- d2H_f-8*d18O_f-dxs_i
D17O_change <- get_D17O(d18O_f,d17O_f)-D17O_i

#
for (i in 1:length(temp)){
  for (j in 1:length(rh)){
    for (k in 1:length(d)){
      if (rr_result_matrix[i,j,k,1]==100){
        d18O_change[i,j,k] <- NA
        dxs_change[i,j,k] <- NA
        D17O_change[i,j,k] <- NA
      }
    }
  }
}
d18O_change[,60:61,] <- NA;d18O_change[1,,] <- NA
dxs_change[,60:61,] <- NA;dxs_change[1,,] <- NA
D17O_change[,60:61,] <- NA;D17O_change[1,,] <- NA

#
lo_bd <- vector()
for (i in 2:length(temp)){
  lo_bd[i] <- min(which(!is.na(d18O_change[i,,1])))
}

par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
# 1
plot(NA,NA,xlim=c(1,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(5,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp,rh,d18O_change[,,3],add=T,levels=c(0.5,1,2),
        col="blue",labcex=0.8,vfont=NULL,method="edge")

contour(temp,rh,d18O_change[,,2],add=T,levels=c(0.5,1,2,4),
        labcex=0.8,vfont=NULL)

contour(temp,rh,d18O_change[,,1],add=T,levels=c(1,2,4,8),
        col="brown1",labcex=0.8,vfont=NULL)
lines(temp[14:31],rh[lo_bd][14:31],lty=2,col="brown1")

text(0,102.5/100,expression(paste("(a) ",Delta,"(",delta^18,"O)"," (\u2030)")),cex=0.9,pos=4)
legend("topleft",c("Dm = 0.8 mm","Dm = 1.6 mm","Dm = 2.4 mm","NA"),
       col=c("brown1","black","blue","brown1"),lwd=1,cex=0.9,lty=c(1,1,1,2))

points(20,0.71,bg="brown1",pch=21,cex=1.3)
points(20,0.47,bg="blue",pch=21,cex=1.3)

# 2
plot(NA,NA,xlim=c(1,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(5,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp,rh,dxs_change[,,3],add=T,levels=c(-2.5,-5,-10),
        col="blue",labcex=0.8,vfont=NULL,method="edge")

contour(temp,rh,dxs_change[,,2],add=T,levels=c(-5,-10,-20),
        labcex=0.8,vfont=NULL)

contour(temp,rh,dxs_change[,,1],add=T,levels=c(-5,-10,-20,-40),
        col="brown1",labcex=0.8,vfont=NULL)
lines(temp[14:31],rh[lo_bd][14:31],lty=2,col="brown1")

text(0,102.2/100,expression(paste("(b) ",Delta,"(d-excess) (","\u2030",")")),cex=0.9,pos=4)

points(20,0.71,bg="brown1",pch=21,cex=1.3)
points(20,0.47,bg="blue",pch=21,cex=1.3)

# 3
plot(NA,NA,xlim=c(1,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(5,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp,rh,D17O_change[,,3],add=T,levels=c(-20,-10,-5,-2.5),
        col="blue",labcex=0.8,vfont=NULL)

contour(temp,rh,D17O_change[,,2],add=T,levels=c(-40,-20,-10,-5),
        labcex=0.8,vfont=NULL)

contour(temp,rh,D17O_change[,,1],add=T,levels=c(-80,-40,-20,-10,-5),
        col="brown1",labcex=0.8,vfont=NULL)
lines(temp[14:31],rh[lo_bd][14:31],lty=2,col="brown1")

text(0,102.5/100,expression(paste("(c) ",Delta,"(",Delta*minute^17,"O)"," (per meg)")),cex=0.9,pos=4)

points(20,0.71,bg="brown1",pch=21,cex=1.3)
points(20,0.47,bg="blue",pch=21,cex=1.3)

# 4
plot(NA,NA,xlim=c(1,30),ylim=c(0.4,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(5,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0.4,1,0.1),cex.axis=0.9,padj=1.1,tck=-0.02,labels=seq(40,100,10))
mtext("Temperature (\u00B0C)",1,padj=2.2,cex=0.8)
mtext("RH (%)",2,padj=-2.5,cex=0.8)

contour(temp,rh,D17O_change[,,3]/dxs_change[,,3],add=T,levels=c(0.7,0.9,1.1,1.3),
        col="blue",labcex=0.8,vfont=NULL)

contour(temp,rh,D17O_change[,,2]/dxs_change[,,2],add=T,levels=c(0.8,1,1.2,1.4,1.6),
        labcex=0.8,vfont=NULL)

contour(temp,rh,D17O_change[,,1]/dxs_change[,,1],add=T,levels=c(1,1.5,2,2.5),
        col="brown1",labcex=0.8,vfont=NULL)
lines(temp[14:31],rh[lo_bd][14:31],lty=2,col="brown1")

text(0,102.5/100,expression(paste("(d) ",Delta,"(",Delta*minute^17,"O)/",
                                 Delta,"(d-excess)")),cex=0.9,pos=4)

points(20,0.71,bg="brown1",pch=21,cex=1.3)
text(20,0.68,1.71,col="brown1",cex=0.9)
points(20,0.47,bg="blue",pch=21,cex=1.3)
text(20,0.44,0.94,col="blue",cex=0.9)