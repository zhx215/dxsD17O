rm(list=ls())
library(openxlsx)
library(pracma)

get_d17O <- function(d18O,D17O){
  1000*(exp(D17O/10^6+0.528*log(d18O/1000+1))-1)
}
get_D17O <- function(d18O,d17O){
  (log(d17O/1000+1)-0.528*log(d18O/1000+1))*10^6
}

#
closure <- function(sst,h,o18=0,o17=0,o2=0,m=0.25,v18o=0,v17o=0,v2o=0){
  eq18 <- function(temp){1/exp(1137/temp^2-0.4156/temp-0.00207)}
  eq17 <- function(temp){1/exp(1137/temp^2-0.4156/temp-0.00207)^0.529}
  eq2 <- function(temp){1/exp(24844/temp^2-76.248/temp+0.05261)}
  k18 <- (1/1.0285)^m
  k17 <- (1/1.0285^0.518)^m
  k2 <- (1/1.0251)^m
  e18 <- (eq18(sst)*k18*(o18+1000)-v18o*h*k18)/(1-h*(1-k18))-1000
  e17 <- (eq17(sst)*k17*(o17+1000)-v17o*h*k17)/(1-h*(1-k17))-1000
  e2 <- (eq2(sst)*k2*(o2+1000)-v2o*h*k2)/(1-h*(1-k2))-1000
  return(c(e18,e17,e2))
}

sst_list <- seq(0,30,1)
rh_list <- seq(20,110,1)/100
e_array <- array(0,dim=c(3,length(sst_list),length(rh_list)))
for (i in 1:length(sst_list)){
  for (j in 1:length(rh_list)){
    e_array[,i,j] <- closure(sst_list[i]+273.15,rh_list[j])
  }
}
rh_array <- t(replicate(length(sst_list),rh_list))
sst_array <- replicate(length(rh_list),sst_list)

mlr <- lm(as.vector(e_array[3,,]-8*e_array[1,,])~
            as.vector(rh_array*100)+as.vector(sst_array))
summary(mlr)

mlr <- lm(as.vector(get_D17O(e_array[1,,],e_array[2,,]))~
            as.vector(rh_array*100)+as.vector(sst_array))
summary(mlr)

# d-excess = -0.50*RH+0.35*SST+42 # R2=0.9998
# D17O = -0.72*RH+0.09*SST+60 # R2=0.9996

# m effect
e_array2 <- e_array3 <- array(0,dim=c(3,2,length(rh_list)))
for (j in 1:length(rh_list)){
  e_array2[,1,j] <- closure(sst_list[16]+273.15,rh_list[j],0,0,0,0.15)
  e_array2[,2,j] <- closure(sst_list[16]+273.15,rh_list[j],0,0,0,0.35)
  e_array3[,1,j] <- closure(sst_list[16]+273.15,rh_list[j],0,0,0,0.25,0,0.014,10)
  e_array3[,2,j] <- closure(sst_list[16]+273.15,rh_list[j],0,0,0,0.25,0,-0.014,-10)
}

# data
data_D17O <- read.xlsx("C:\\Users\\zhyxi\\Desktop\\17O project\\ocean vapor data.xlsx",sheet=1)
data_dxs <- read.xlsx("C:\\Users\\zhyxi\\Desktop\\17O project\\ocean vapor data.xlsx",sheet=2)
data_offset <- read.xlsx("C:\\Users\\zhyxi\\Desktop\\17O project\\ocean vapor data.xlsx",sheet=3)

#
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
#1
plot(NA,NA,pch=21,bg="white",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",
     xlim=c(20,110),
     ylim=c(-20,50))
axis(1,seq(20,100,20),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,50,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("Oceanic RH (%)",1,padj=2.2,cex=0.8)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.8)

polygon(c(rh_list*100,rh_list[91:1]*100),
        c(e_array[3,1,]-8*e_array[1,1,],e_array[3,31,91:1]-8*e_array[1,31,91:1]),
        col="darkolivegreen1",border=NA)

lines(data_dxs[,1],data_dxs[,2],lwd=2)
lines(data_dxs[,1],data_dxs[,5])
lines(data_dxs[,1],data_dxs[,6])

lines(rh_list*100,e_array2[3,1,]-8*e_array2[1,1,],lty=2,col="blue")
lines(rh_list[6:91]*100,e_array2[3,2,6:91]-8*e_array2[1,2,6:91],lty=2,col="blue")
lines(rh_list*100,e_array3[3,1,]-8*e_array3[1,1,],lty=4,col="red")
lines(rh_list*100,e_array3[3,2,]-8*e_array3[1,2,],lty=4,col="red")

points(data_dxs[,7],data_dxs[,8],pch=21,bg="white")

box(lwd=1)
legend("bottomleft",c("regression line",
                      expression(paste("PI (2",sigma,")")),
                      "SST effect",
                      "aerodynamic effect",
                      "gradient effect"),
       lwd=c(2,1,NA,1,1),lty=c(1,1,NA,2,4),cex=0.9,col=c(rep("black",3),"blue","red"))
rect(23.5,-9.8,31,-7.8,col="darkolivegreen1",border=NA)
text(22,53,"(a)",cex=0.9)

text(100,47,"m = 0.35",cex=0.9,col="blue")
text(100,43,"m = 0.15",cex=0.9,col="blue")
text(99,39,expression(paste(""^2,"g = 10\u2030")),cex=0.9,col="red")
text(99,35,expression(paste(""^2,"g = -10\u2030")),cex=0.9,col="red")
text(89,47,"1",cex=0.9);points(89,47,cex=2)
text(89,43,"2",cex=0.9);points(89,43,cex=2)
text(86,39,"3",cex=0.9);points(86,39,cex=2)
text(86,35,"4",cex=0.9);points(86,35,cex=2)

points(29,47,cex=2,pch=21,bg="white");text(29,47,"1",cex=0.9)
points(24,20,cex=2,pch=21,bg="white");text(24,20,"2",cex=0.9)
points(107,6,cex=2,pch=21,bg="white");text(107,6,"4",cex=0.9)
points(107,-15,cex=2,pch=21,bg="white");text(107,-15,"3",cex=0.9)

#2
plot(NA,NA,pch=21,bg="white",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",
     xlim=c(20,110),
     ylim=c(-20,50))
axis(1,seq(20,100,20),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,50,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("Oceanic RH (%)",1,padj=2.2,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

polygon(c(rh_list*100,rh_list[91:1]*100),
        c(get_D17O(e_array[1,1,],e_array[2,1,]),get_D17O(e_array[1,31,91:1],e_array[2,31,91:1])),
        col="darkolivegreen1",border=NA)

lines(data_D17O[,1],data_D17O[,2],lwd=2)
lines(data_D17O[,1],data_D17O[,5])
lines(data_D17O[,1],data_D17O[,6])

climate_error <- sqrt(((data_D17O[,6]-data_D17O[,5])/4)^2-5^2)
lines(data_D17O[,1],data_D17O[,2]+2*climate_error,col="gray")
lines(data_D17O[,1],data_D17O[,2]-2*climate_error,col="gray")

lines(rh_list*100,get_D17O(e_array2[1,1,],e_array2[2,1,]),lty=2,col="blue")
lines(rh_list[23:89]*100,get_D17O(e_array2[1,2,23:89],e_array2[2,2,23:89]),lty=2,col="blue")
lines(rh_list[1:76]*100,get_D17O(e_array3[1,1,1:76],e_array3[2,1,1:76]),lty=4,col="red")
lines(rh_list*100,get_D17O(e_array3[1,2,],e_array3[2,2,]),lty=4,col="red")

points(data_D17O[,7],data_D17O[,8],pch=21,bg="white")

box(lwd=1)
legend("bottomleft",expression(paste("PI (2",sigma,") w/o analytical error")),
       lwd=1,lty=1,cex=0.9,col="gray",text.width=48)
text(22,53,"(b)",cex=0.9)

text(100,47,"m = 0.35",cex=0.9,col="blue")
text(100,43,"m = 0.15",cex=0.9,col="blue")
text(96,39,expression(paste(""^17,"g = 0.014\u2030")),cex=0.9,col="red")
text(96,35,expression(paste(""^17,"g = -0.014\u2030")),cex=0.9,col="red")
text(89,47,"1",cex=0.9);points(89,47,cex=2)
text(89,43,"2",cex=0.9);points(89,43,cex=2)
text(80,39,"3",cex=0.9);points(80,39,cex=2)
text(80,35,"4",cex=0.9);points(80,35,cex=2)

points(45,47,cex=2,pch=21,bg="white");text(45,47,"1",cex=0.9)
points(24,22,cex=2,pch=21,bg="white");text(24,22,"2",cex=0.9)
points(107,0,cex=2,pch=21,bg="white");text(107,0,"4",cex=0.9)
points(95,-18,cex=2,pch=21,bg="white");text(95,-18,"3",cex=0.9)

#3
plot(NA,NA,pch=21,bg="white",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",
     xlim=c(-10,15),
     ylim=c(-20,25))
axis(1,seq(-10,15,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,20,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("d-excess residual (","\u2030",")")),1,padj=1.9,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O residual"," (per meg)")),2,padj=-1.3,cex=0.8)

points(data_offset[,2],data_offset[,1],pch=21,bg="white")
text(-9.42,27,"(c)",cex=0.9)