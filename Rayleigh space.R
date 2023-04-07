rm(list=ls())
library(alphahull)

# fractionation
a18eff <- function(temp,mode="k",l=0.004) {
  if (temp >= 273.15) {
    eq <- exp(1137/temp^2-0.4156/temp-0.00207)
    value <- eq
  }
  else if (temp <= 250.15) {
    eq <- exp(11.839/temp-0.028224)
    Si <- 1-l*(temp-273.15)
    diff16_18 <- 1.0285
    k <- Si/(eq*diff16_18*(Si-1)+1)
    if (mode=="no k") {k <- 1}
    value <- k*eq
  }
  else {
    eq_273 <- exp(1137/273.15^2-0.4156/273.15-0.00207)
    eq_274 <- exp(1137/274.15^2-0.4156/274.15-0.00207)
    eq_250 <- exp(11.839/250.15-0.028224)
    eq_249 <- exp(11.839/249.15-0.028224)
    Si_250 <- 1-l*(250.15-273.15)
    Si_249 <- 1-l*(249.15-273.15)
    diff16_18 <- 1.0285
    k_250 <- Si_250/(eq_250*diff16_18*(Si_250-1)+1)
    k_249 <- Si_249/(eq_249*diff16_18*(Si_249-1)+1)
    if (mode=="no k") {k_250 <- k_249 <- 1}
    x1 <- 250.15
    x2 <- 273.15
    k1 <- eq_250*k_250-eq_249*k_249
    k2 <- eq_274-eq_273
    y1 <- eq_250*k_250
    y2 <- eq_273
    answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
    value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
  }
  return(value)
}

a17eff <- function(temp,mode="k",l=0.004) {
  if (temp >= 273.15) {
    eq <- exp(1137/temp^2-0.4156/temp-0.00207)^0.529
    value <- eq
  }
  else if (temp <= 250.15) {
    eq <- exp(11.839/temp-0.028224)^0.529
    Si <- 1-l*(temp-273.15)
    diff16_17 <- 1.0285^0.518
    k <- Si/(eq*diff16_17*(Si-1)+1)
    if (mode=="no k") {k <- 1}
    value <- k*eq
  }
  else {
    eq_273 <- exp(1137/273.15^2-0.4156/273.15-0.00207)^0.529
    eq_274 <- exp(1137/274.15^2-0.4156/274.15-0.00207)^0.529
    eq_250 <- exp(11.839/250.15-0.028224)^0.529
    eq_249 <- exp(11.839/249.15-0.028224)^0.529
    Si_250 <- 1-l*(250.15-273.15)
    Si_249 <- 1-l*(249.15-273.15)
    diff16_17 <- 1.0285^0.518
    k_250 <- Si_250/(eq_250*diff16_17*(Si_250-1)+1)
    k_249 <- Si_249/(eq_249*diff16_17*(Si_249-1)+1)
    if (mode=="no k") {k_250 <- k_249 <- 1}
    x1 <- 250.15
    x2 <- 273.15
    k1 <- eq_250*k_250-eq_249*k_249
    k2 <- eq_274-eq_273
    y1 <- eq_250*k_250
    y2 <- eq_273
    answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
    value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
  }
  return(value)
}

a2eff <- function(temp,mode="k",l=0.004) {
  if (temp >= 273.15) {
    eq <- exp(24844/temp^2-76.248/temp+0.05261)
    value <- eq
  }
  else if (temp <= 250.15) {
    eq <- exp(16289/temp^2-0.0945)
    Si <- 1-l*(temp-273.15)
    diff1_2 <- 1.0251
    k <- Si/(eq*diff1_2*(Si-1)+1)
    if (mode=="no k") {k <- 1}
    value <- k*eq
  }
  else {
    eq_273 <- exp(24844/273.15^2-76.248/273.15+0.05261)
    eq_274 <- exp(24844/274.15^2-76.248/274.15+0.05261)
    eq_250 <- exp(16289/250.15^2-0.0945)
    eq_249 <- exp(16289/249.15^2-0.0945)
    Si_250 <- 1-l*(250.15-273.15)
    Si_249 <- 1-l*(249.15-273.15)
    diff1_2 <- 1.0251
    k_250 <- Si_250/(eq_250*diff1_2*(Si_250-1)+1)
    k_249 <- Si_249/(eq_249*diff1_2*(Si_249-1)+1)
    if (mode=="no k") {k_250 <- k_249 <- 1}
    x1 <- 250.15
    x2 <- 273.15
    k1 <- eq_250*k_250-eq_249*k_249
    k2 <- eq_274-eq_273
    y1 <- eq_250*k_250
    y2 <- eq_273
    answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
    value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
  }
  return(value)
}

# 17O converter
get_D17O <- function(d18O,d17O){
  (log(d17O/1000+1)-0.528*log(d18O/1000+1))*10^6
}

get_d17O <- function(d18O,D17O){
  1000*(exp(D17O/10^6+0.528*log(d18O/1000+1))-1)
}

# Rayleigh
CC_relation <- function(temp){
  es <- 611.2*exp(17.67*temp/(temp+243.5))
  return(es)
}

Rayleigh_computer4b <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,step,advection=T,lambda){
  if (advection==T) {times <- 1} else (times <- 0.5)
  td_list <- seq(initial_td,final_td,length.out=step)
  f_list <- CC_relation(td_list)
  f_list <- f_list/f_list[1]
  
  d18Ov_list <- d18Ov
  d2Hv_list <- dxsv+8*d18Ov
  d17Ov_list <- get_d17O(d18Ov,D17Ov)
  
  d18Op_list <- (d18Ov_list+1000)*a18eff(initial_td+273.15,"k",lambda)-1000
  d2Hp_list <- (d2Hv_list+1000)*a2eff(initial_td+273.15,"k",lambda)-1000
  d17Op_list <- (d17Ov_list+1000)*a17eff(initial_td+273.15,"k",lambda)-1000
  
  for (i in 2:step){
    d18Ov_list[i] <- (d18Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a18eff(td_list[i]+273.15,"k",lambda)^times-1)-1000
    d17Ov_list[i] <- (d17Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a17eff(td_list[i]+273.15,"k",lambda)^times-1)-1000
    d2Hv_list[i] <- (d2Hv_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a2eff(td_list[i]+273.15,"k",lambda)^times-1)-1000
    
    d18Op_list[i] <- (d18Ov_list[i]+1000)*a18eff(td_list[i]+273.15,"k",lambda)-1000
    d17Op_list[i] <- (d17Ov_list[i]+1000)*a17eff(td_list[i]+273.15,"k",lambda)-1000
    d2Hp_list[i] <- (d2Hv_list[i]+1000)*a2eff(td_list[i]+273.15,"k",lambda)-1000
  }
  
  return(list(d18Op_list,d17Op_list,d2Hp_list,
              d18Ov_list,d17Ov_list,d2Hv_list,f_list,td_list))
}

### T > 0 result ###
td_initial <- seq(0,28,1)
d18Ov_initial <- seq(-25,-8,1)
dxsp_change <- D17Op_change <- array(NA,dim=c(3,length(td_initial),length(d18Ov_initial)))

for (i in 1:length(td_initial)){
  for (j in 1:length(d18Ov_initial)){
    result <- Rayleigh_computer4b(d18Ov_initial[j],10,10,td_initial[i],0,1000,T,0.004)
    f0.7 <- which(result[[7]]<0.7)[1]
    f0.4 <- which(result[[7]]<0.4)[1]
    dxsp_change[1,i,j] <- (result[[3]][1]-8*result[[1]][1])-10
    D17Op_change[1,i,j] <- get_D17O(result[[1]][1],result[[2]][1])-10
    dxsp_change[2,i,j] <- (result[[3]][f0.7]-8*result[[1]][f0.7])-10
    D17Op_change[2,i,j] <- get_D17O(result[[1]][f0.7],result[[2]][f0.7])-10
    dxsp_change[3,i,j] <- (result[[3]][f0.4]-8*result[[1]][f0.4])-10
    D17Op_change[3,i,j] <- get_D17O(result[[1]][f0.4],result[[2]][f0.4])-10
  }
}

# plot
colors <- c("brown2","burlywood2","chartreuse2","cadetblue2")
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)

# 1
plot(NA,NA,xlim=c(0,28),ylim=c(-25,-8),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,25,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-25,-10,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Initial T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste("Initial vapor ",delta^18,"O (","\u2030",")")),2,padj=-1.4,cex=0.8)

contour(td_initial,d18Ov_initial,dxsp_change[1,,],add=T,labcex=0.8,levels=seq(-8,12,4),
        vfont=NULL,col=colors[1],method="edge")
contour(td_initial,d18Ov_initial,dxsp_change[2,,],add=T,labcex=0.8,levels=seq(-6,4,2),
        vfont=NULL,col=colors[2],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change[3,,],add=T,labcex=0.8,levels=seq(-1.5,0,0.5),
        vfont=NULL,col=colors[3],method="flattest")
text(-0.9,-7.3,expression(paste("(a) ",Delta,"(d-excess) (","\u2030",")")),cex=0.9,pos=4)

# 2
plot(NA,NA,xlim=c(0,28),ylim=c(-25,-8),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,25,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-25,-10,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Initial T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste("Initial vapor ",delta^18,"O (","\u2030",")")),2,padj=-1.4,cex=0.8)

contour(td_initial,d18Ov_initial,D17Op_change[1,,],add=T,drawlabels=F,levels=seq(10,11,1),
        vfont=NULL,col=colors[1],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change[2,,],add=T,drawlabels=F,levels=seq(11,13,1),
        vfont=NULL,col=colors[2],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change[3,,],add=T,drawlabels=F,levels=seq(13,15,1),
        vfont=NULL,col=colors[3],method="edge")

text(5.9,-8.3,"11",cex=1,pos=3,col=colors[1])
text(17,-8.3,"10",cex=1,pos=3,col=colors[1])
text(7.4,-8.3,"13",cex=1,pos=3,col=colors[2])
text(13.5,-8.3,"12",cex=1,pos=3,col=colors[2])
text(20.8,-8.3,"11",cex=1,pos=3,col=colors[2])
text(15,-8.3,"15",cex=1,pos=3,col=colors[3])
text(19.2,-8.3,"14",cex=1,pos=3,col=colors[3])
text(24,-8.3,"13",cex=1,pos=3,col=colors[3])
text(-1,-6.7,expression(paste("(b) ",Delta,"(",Delta*minute^17,"O)"," (per meg)")),cex=0.9,pos=4)

# 3 
plot(NA,NA,xlim=c(-10,15),ylim=c(0,25),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-10,15,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,25,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(Delta,"(d-excess) (","\u2030",")")),1,padj=1.7,cex=0.8)
mtext(expression(paste(Delta,"(",Delta*minute^17,"O)"," (per meg)")),2,padj=-1.4,cex=0.8)

polygon(c(as.vector(dxsp_change[1,1,]),
          as.vector(dxsp_change[1,,18]),
          as.vector(dxsp_change[1,29,18:1]),
          as.vector(dxsp_change[1,29:1,1])),
        c(as.vector(D17Op_change[1,1,]),
          as.vector(D17Op_change[1,,18]),
          as.vector(D17Op_change[1,29,18:1]),
          as.vector(D17Op_change[1,29:1,1])),border=colors[1],col=NA)
polygon(c(as.vector(dxsp_change[2,1,]),
          as.vector(dxsp_change[2,,18]),
          as.vector(dxsp_change[2,29,18:1]),
          as.vector(dxsp_change[2,29:1,1])),
        c(as.vector(D17Op_change[2,1,]),
          as.vector(D17Op_change[2,,18]),
          as.vector(D17Op_change[2,29,18:1]),
          as.vector(D17Op_change[2,29:1,1])),border=colors[2],col=NA)
polygon(c(as.vector(dxsp_change[3,1,]),
          as.vector(dxsp_change[3,1:29,18]),
          as.vector(dxsp_change[3,29,18:1]),
          as.vector(dxsp_change[3,29:1,1])),
        c(as.vector(D17Op_change[3,1,]),
          as.vector(D17Op_change[3,1:29,18]),
          as.vector(D17Op_change[3,29,18:1]),
          as.vector(D17Op_change[3,29:1,1])),border=colors[3],col=NA)
text(-10.8,26.2,"(c)",cex=0.9,pos=4)

plot.new()
legend("topleft",c("f = 1",
                   expression(paste("f = 0.7 (initial T"[d]," > 5\u00B0C)")),
                   expression(paste("f = 0.4 (initial T"[d]," > 13\u00B0C)"))),
       col=colors[1:3],lwd=1,cex=0.9)






### T < 0 result ###
td_initial <- seq(-20,0,1)
d18Ov_initial <- seq(-40,-15,1)
dxsp_change <- D17Op_change <- array(NA,dim=c(4,length(td_initial),length(d18Ov_initial)))

for (i in 1:length(td_initial)){
  for (j in 1:length(d18Ov_initial)){
    result <- Rayleigh_computer4b(d18Ov_initial[j],10,10,td_initial[i],-40,1000,T,0.004)
    f0.7 <- which(result[[7]]<0.7)[1]
    f0.4 <- which(result[[7]]<0.4)[1]
    f0.2 <- which(result[[7]]<0.2)[1]
    dxsp_change[1,i,j] <- (result[[3]][1]-8*result[[1]][1])-10
    D17Op_change[1,i,j] <- get_D17O(result[[1]][1],result[[2]][1])-10
    dxsp_change[2,i,j] <- (result[[3]][f0.7]-8*result[[1]][f0.7])-10
    D17Op_change[2,i,j] <- get_D17O(result[[1]][f0.7],result[[2]][f0.7])-10
    dxsp_change[3,i,j] <- (result[[3]][f0.4]-8*result[[1]][f0.4])-10
    D17Op_change[3,i,j] <- get_D17O(result[[1]][f0.4],result[[2]][f0.4])-10
    dxsp_change[4,i,j] <- (result[[3]][f0.2]-8*result[[1]][f0.2])-10
    D17Op_change[4,i,j] <- get_D17O(result[[1]][f0.2],result[[2]][f0.2])-10
  }
}

# plot
colors <- c("brown2","burlywood2","chartreuse2","cadetblue2")
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)

# 1
plot(NA,NA,xlim=c(-20,0),ylim=c(-40,-15),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-20,0,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,-15,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Initial T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste("Initial vapor ",delta^18,"O (","\u2030",")")),2,padj=-1.4,cex=0.8)

contour(td_initial,d18Ov_initial,dxsp_change[1,,],add=T,labcex=0.8,levels=seq(-12,16,4),
        vfont=NULL,col=colors[1],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change[2,,],add=T,labcex=0.8,levels=seq(-8,8,4),
        vfont=NULL,col=colors[2],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change[3,,],add=T,labcex=0.8,levels=seq(-5,-3,1),
        vfont=NULL,col=colors[3],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change[4,,],add=T,labcex=0.8,levels=seq(-12,4,4),
        vfont=NULL,col=colors[4],method="flattest")
text(-20.7,-14,expression(paste("(a) ",Delta,"(d-excess) (","\u2030",")")),cex=0.9,pos=4)

# 2
plot(NA,NA,xlim=c(-20,0),ylim=c(-40,-15),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-20,0,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,-15,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Initial T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste("Initial vapor ",delta^18,"O (","\u2030",")")),2,padj=-1.4,cex=0.8)

contour(td_initial,d18Ov_initial,D17Op_change[1,,],add=T,drawlabels=F,levels=seq(15,55,10),
        vfont=NULL,col=colors[1],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change[2,,],add=T,drawlabels=F,levels=c(seq(20,40,10),49.8),
        vfont=NULL,col=colors[2],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change[3,,],add=T,drawlabels=F,levels=seq(35,40,5),
        vfont=NULL,col=colors[3],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change[4,,],add=T,drawlabels=F,levels=seq(20,40,10),
        vfont=NULL,col=colors[4],method="edge")

text(-19.4,-15.4,"55",cex=1,pos=3,col=colors[1])
text(-15.6,-15.4,"45",cex=1,pos=3,col=colors[1])
text(-12.8,-15.4,"35",cex=1,pos=3,col=colors[1])
text(-8.9,-15.4,"25",cex=1,pos=3,col=colors[1])
text(-4.5,-15.4,"15",cex=1,pos=3,col=colors[1])

text(-16.8,-15.4,"50",cex=1,pos=3,col=colors[2])
text(-11.6,-15.4,"40",cex=1,pos=3,col=colors[2])
text(-6.9,-15.4,"30",cex=1,pos=3,col=colors[2])
text(-2,-15.4,"20",cex=1,pos=3,col=colors[2])

text(-5.7,-15.4,"40",cex=1,pos=3,col=colors[3])
text(-0.8,-15.4,"35",cex=1,pos=3,col=colors[3])

text(-18,-15.4,"20",cex=1,pos=3,col=colors[4])
text(-10.4,-15.4,"30",cex=1,pos=3,col=colors[4])
text(-3.3,-15.4,"40",cex=1,pos=3,col=colors[4])

text(-20.7,-13.1,expression(paste("(b) ",Delta,"(",Delta*minute^17,"O)"," (per meg)")),cex=0.9,pos=4)

# 3 
plot(NA,NA,xlim=c(-25,25),ylim=c(10,60),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-20,20,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(10,60,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(Delta,"(d-excess) (","\u2030",")")),1,padj=1.9,cex=0.8)
mtext(expression(paste(Delta,"(",Delta*minute^17,"O)"," (per meg)")),2,padj=-1.4,cex=0.8)

polygon(c(as.vector(dxsp_change[1,1,]),
          as.vector(dxsp_change[1,,26]),
          as.vector(dxsp_change[1,21,26:1]),
          as.vector(dxsp_change[1,21:1,1])),
        c(as.vector(D17Op_change[1,1,]),
          as.vector(D17Op_change[1,,26]),
          as.vector(D17Op_change[1,21,26:1]),
          as.vector(D17Op_change[1,21:1,1])),border=colors[1],col=NA)
polygon(c(as.vector(dxsp_change[2,1,]),
          as.vector(dxsp_change[2,,26]),
          as.vector(dxsp_change[2,21,26:1]),
          as.vector(dxsp_change[2,21:1,1])),
        c(as.vector(D17Op_change[2,1,]),
          as.vector(D17Op_change[2,,26]),
          as.vector(D17Op_change[2,21,26:1]),
          as.vector(D17Op_change[2,21:1,1])),border=colors[2],col=NA)
plot(ahull(as.vector(dxsp_change[3,,]),as.vector(D17Op_change[3,,]),alpha=2.5),
     wpoints=F,add=T,lwd=1,col=colors[3])
polygon(c(as.vector(dxsp_change[4,1,]),
          as.vector(dxsp_change[4,,26]),
          as.vector(dxsp_change[4,21,26:1]),
          as.vector(dxsp_change[4,21:1,1])),
        c(as.vector(D17Op_change[4,1,]),
          as.vector(D17Op_change[4,,26]),
          as.vector(D17Op_change[4,21,26:1]),
          as.vector(D17Op_change[4,21:1,1])),border=colors[4],col=NA)
text(-26.7,62.4,"(c)",cex=0.9,pos=4)

plot.new()
legend("topleft",c("f = 1","f = 0.7","f = 0.4","f = 0.2"),
       col=colors[1:4],lwd=1,cex=0.9)