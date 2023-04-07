rm(list=ls())
library(alphahull)

# fractionation
a18eff <- function(temp,vmt,mit,l=0.004) {
  temp1 <- temp[temp<=mit]
  temp2 <- temp[(temp<vmt)&(temp>mit)]
  temp3 <- temp[temp>=vmt]
  
  eq <- exp(1137/temp3^2-0.4156/temp3-0.00207)
  value3 <- eq

  eq <- exp(11.839/temp1-0.028224)
  Si <- 1-l*(temp1-273.15)
  diff16_18 <- 1.0285
  k <- Si/(eq*diff16_18*(Si-1)+1)
  value1 <- k*eq

  eq_273 <- exp(1137/vmt^2-0.4156/vmt-0.00207)
  eq_274 <- exp(1137/(vmt+1)^2-0.4156/(vmt+1)-0.00207)
  eq_250 <- exp(11.839/mit-0.028224)
  eq_249 <- exp(11.839/(mit-1)-0.028224)
  Si_250 <- 1-l*(mit-273.15)
  Si_249 <- 1-l*((mit-1)-273.15)
  diff16_18 <- 1.0285
  k_250 <- Si_250/(eq_250*diff16_18*(Si_250-1)+1)
  k_249 <- Si_249/(eq_249*diff16_18*(Si_249-1)+1)
  x1 <- mit
  x2 <- vmt
  k1 <- eq_250*k_250-eq_249*k_249
  k2 <- eq_274-eq_273
  y1 <- eq_250*k_250
  y2 <- eq_273
  answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
  value2 <- answer[1]*temp2^3+answer[2]*temp2^2+answer[3]*temp2+answer[4]

  return(c(value3,value2,value1))
}

a17eff <- function(temp,vmt,mit,l=0.004) {
  temp1 <- temp[temp<=mit]
  temp2 <- temp[(temp<vmt)&(temp>mit)]
  temp3 <- temp[temp>=vmt]
  eq <- exp(1137/temp3^2-0.4156/temp3-0.00207)^0.529
  value3 <- eq

  eq <- exp(11.839/temp1-0.028224)^0.529
  Si <- 1-l*(temp1-273.15)
  diff16_17 <- 1.0285^0.518
  k <- Si/(eq*diff16_17*(Si-1)+1)
  value1 <- k*eq

  eq_273 <- exp(1137/vmt^2-0.4156/vmt-0.00207)^0.529
  eq_274 <- exp(1137/(vmt+1)^2-0.4156/(vmt+1)-0.00207)^0.529
  eq_250 <- exp(11.839/mit-0.028224)^0.529
  eq_249 <- exp(11.839/(mit-1)-0.028224)^0.529
  Si_250 <- 1-l*(mit-273.15)
  Si_249 <- 1-l*((mit-1)-273.15)
  diff16_17 <- 1.0285^0.518
  k_250 <- Si_250/(eq_250*diff16_17*(Si_250-1)+1)
  k_249 <- Si_249/(eq_249*diff16_17*(Si_249-1)+1)
  x1 <- mit
  x2 <- vmt
  k1 <- eq_250*k_250-eq_249*k_249
  k2 <- eq_274-eq_273
  y1 <- eq_250*k_250
  y2 <- eq_273
  answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
  value2 <- answer[1]*temp2^3+answer[2]*temp2^2+answer[3]*temp2+answer[4]

  return(c(value3,value2,value1))
}

a2eff <- function(temp,vmt,mit,l=0.004,mode="k") {
  temp1 <- temp[temp<=mit]
  temp2 <- temp[(temp<vmt)&(temp>mit)]
  temp3 <- temp[temp>=vmt]
  
  eq <- exp(24844/temp3^2-76.248/temp3+0.05261)
  value3 <- eq

  eq <- exp(16289/temp1^2-0.0945)
  Si <- 1-l*(temp1-273.15)
  diff1_2 <- 1.0251
  k <- Si/(eq*diff1_2*(Si-1)+1)
  value1 <- k*eq

  eq_273 <- exp(24844/vmt^2-76.248/vmt+0.05261)
  eq_274 <- exp(24844/(vmt+1)^2-76.248/(vmt+1)+0.05261)
  eq_250 <- exp(16289/mit^2-0.0945)
  eq_249 <- exp(16289/(mit-1)^2-0.0945)
  Si_250 <- 1-l*(mit-273.15)
  Si_249 <- 1-l*((mit-1)-273.15)
  diff1_2 <- 1.0251
  k_250 <- Si_250/(eq_250*diff1_2*(Si_250-1)+1)
  k_249 <- Si_249/(eq_249*diff1_2*(Si_249-1)+1)
  x1 <- mit
  x2 <- vmt
  k1 <- eq_250*k_250-eq_249*k_249
  k2 <- eq_274-eq_273
  y1 <- eq_250*k_250
  y2 <- eq_273
  answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
  value2 <- answer[1]*temp2^3+answer[2]*temp2^2+answer[3]*temp2+answer[4]

  return(c(value3,value2,value1))
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

#
Rayleigh_computer5 <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,step,advection=T,vmt,mit,lambda){
  if (advection==T) {times <- 1} else {times <- 0.5}
  td_list <- seq(initial_td,final_td,length.out=step)
  f_list <- CC_relation(td_list)
  f_list <- f_list/f_list[1]
  
  a18_list <- a18eff(td_list+273.15,vmt,mit,lambda)
  a17_list <- a17eff(td_list+273.15,vmt,mit,lambda)
  a2_list <-a2eff(td_list+273.15,vmt,mit,lambda)
  
  d18Ov_list <- d18Ov
  d2Hv_list <- dxsv+8*d18Ov
  d17Ov_list <- get_d17O(d18Ov,D17Ov)
  
  d18Op_list <- (d18Ov_list+1000)*a18_list[1]-1000
  d2Hp_list <- (d2Hv_list+1000)*a2_list[1]-1000
  d17Op_list <- (d17Ov_list+1000)*a17_list[1]-1000
  
  for (i in 2:step){
    d18Ov_list[i] <- (d18Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a18_list[i]^times-1)-1000
    d17Ov_list[i] <- (d17Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a17_list[i]^times-1)-1000
    d2Hv_list[i] <- (d2Hv_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a2_list[i]^times-1)-1000
    
    d18Op_list[i] <- (d18Ov_list[i]+1000)*a18_list[i]-1000
    d17Op_list[i] <- (d17Ov_list[i]+1000)*a17_list[i]-1000
    d2Hp_list[i] <- (d2Hv_list[i]+1000)*a2_list[i]-1000
  }
  
  return(list(d18Op_list,d17Op_list,d2Hp_list,
              d18Ov_list,d17Ov_list,d2Hv_list,f_list,td_list))
}

# T < 0 result
td_initial <- seq(-20,0,1)
d18Ov_initial <- seq(-40,-15,1)

vmt <- seq(263.15,273.15,1)
mit <- seq(233.15,258.15,1)
lambda <- seq(0.002,0.006,0.0002)

dxsp_change <- D17Op_change <- 
  array(NA,dim=c(length(vmt),length(mit),length(lambda),4,length(td_initial),length(d18Ov_initial)))

for (i in 1:length(vmt)){
  for (j in 1:length(mit)){
    for (k in 1:length(lambda)){
      for (l in 1:length(td_initial)){
        for (m in 1:length(d18Ov_initial)){
          result <- Rayleigh_computer5(d18Ov_initial[m],10,10,td_initial[l],-40,300,T,vmt[i],mit[j],lambda[k])
          f0.7 <- which(result[[7]]<0.7)[1]
          f0.4 <- which(result[[7]]<0.4)[1]
          f0.2 <- which(result[[7]]<0.2)[1]
          dxsp_change[i,j,k,1,l,m] <- (result[[3]][1]-8*result[[1]][1])-10
          D17Op_change[i,j,k,1,l,m] <- get_D17O(result[[1]][1],result[[2]][1])-10
          dxsp_change[i,j,k,2,l,m] <- (result[[3]][f0.7]-8*result[[1]][f0.7])-10
          D17Op_change[i,j,k,2,l,m] <- get_D17O(result[[1]][f0.7],result[[2]][f0.7])-10
          dxsp_change[i,j,k,3,l,m] <- (result[[3]][f0.4]-8*result[[1]][f0.4])-10
          D17Op_change[i,j,k,3,l,m] <- get_D17O(result[[1]][f0.4],result[[2]][f0.4])-10
          dxsp_change[i,j,k,4,l,m] <- (result[[3]][f0.2]-8*result[[1]][f0.2])-10
          D17Op_change[i,j,k,4,l,m] <- get_D17O(result[[1]][f0.2],result[[2]][f0.2])-10
        }
      }
    }
    print(j)
  }
} # around 1 hour
#save.image("C:/Users/zhyxi/Desktop/17O project/snow sensitivity result.RData")

load("C:/Users/zhyxi/Desktop/17O project/snow sensitivity result.RData")

dxsp_change_sd <- D17Op_change_sd <- 
  array(NA,dim=c(4,length(td_initial),length(d18Ov_initial)))
for (i in 1:4){
  for (j in 1:length(td_initial)){
    for (k in 1:length(d18Ov_initial)){
      dxsp_change_sd[i,j,k] <- sd(as.vector(dxsp_change[,,,i,j,k]))
      D17Op_change_sd[i,j,k] <- sd(as.vector(D17Op_change[,,,i,j,k]))
    }
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

contour(td_initial,d18Ov_initial,dxsp_change_sd[1,,],add=T,labcex=0.8,levels=seq(1,3,1),
        vfont=NULL,col=colors[1],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change_sd[2,,],add=T,labcex=0.8,levels=seq(1,3,1),
        vfont=NULL,col=colors[2],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change_sd[3,,],add=T,labcex=0.8,levels=seq(1.4,2,0.2),
        vfont=NULL,col=colors[3],method="flattest")
contour(td_initial,d18Ov_initial,dxsp_change_sd[4,,],add=T,labcex=0.8,levels=seq(1.2,2,0.2),
        vfont=NULL,col=colors[4],method="flattest")
text(-20.7,-14,expression(paste("(a) ",sigma,"(d-excess) (","\u2030",")")),cex=0.9,pos=4)

# 2
plot(NA,NA,xlim=c(-20,0),ylim=c(-40,-15),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-20,0,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,-15,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("Initial T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste("Initial vapor ",delta^18,"O (","\u2030",")")),2,padj=-1.4,cex=0.8)

contour(td_initial,d18Ov_initial,D17Op_change_sd[1,,],add=T,drawlabels=F,levels=seq(2,10,4),
        vfont=NULL,col=colors[1],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change_sd[2,,],add=T,drawlabels=F,levels=seq(2,6,2),
        vfont=NULL,col=colors[2],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change_sd[3,,],add=T,drawlabels=F,levels=seq(3,6,1),
        vfont=NULL,col=colors[3],method="edge")
contour(td_initial,d18Ov_initial,D17Op_change_sd[4,,],add=T,drawlabels=F,levels=seq(6,12,3),
        vfont=NULL,col=colors[4],method="edge")

text(-15.7,-15.4,"10",cex=1,pos=3,col=colors[1])
text(-10.6,-15.4,"6",cex=1,pos=3,col=colors[1])
text(-5.7,-15.4,"2",cex=1,pos=3,col=colors[1])

text(-8.4,-15.4,"6",cex=1,pos=3,col=colors[2])
text(-4.3,-15.4,"4",cex=1,pos=3,col=colors[2])
text(-1.2,-15.4,"2",cex=1,pos=3,col=colors[2])

text(-12.8,-15.4,"3",cex=1,pos=3,col=colors[3])
text(-9.9,-15.4,"4",cex=1,pos=3,col=colors[3])
text(-7.7,-15.4,"5",cex=1,pos=3,col=colors[3])
text(-5,-15.4,"6",cex=1,pos=3,col=colors[3])
text(-2,-15.4,"6",cex=1,pos=3,col=colors[3])

text(-17.4,-15.4,"12",cex=1,pos=3,col=colors[4])
text(-12,-15.4,"9",cex=1,pos=3,col=colors[4])
text(-7,-15.4,"6",cex=1,pos=3,col=colors[4])

text(-20.7,-13.1,expression(paste("(b) ",sigma,"(",Delta*minute^17,"O)"," (per meg)")),cex=0.9,pos=4)

# 3
plot(NA,NA,xlim=c(0,14),ylim=c(0,14),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,14,2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,14,2),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(sigma,"(d-excess) (","\u2030",")")),1,padj=1.9,cex=0.8)
mtext(expression(paste(sigma,"(",Delta*minute^17,"O)"," (per meg)")),2,padj=-1.4,cex=0.8)

polygon(c(as.vector(dxsp_change_sd[1,1,]),
          as.vector(dxsp_change_sd[1,,26]),
          as.vector(dxsp_change_sd[1,21,26:1]),
          as.vector(dxsp_change_sd[1,21:1,1])),
        c(as.vector(D17Op_change_sd[1,1,]),
          as.vector(D17Op_change_sd[1,,26]),
          as.vector(D17Op_change_sd[1,21,26:1]),
          as.vector(D17Op_change_sd[1,21:1,1])),border=colors[1],col=NA)
plot(ahull(as.vector(dxsp_change_sd[2,,]),as.vector(D17Op_change_sd[2,,]),alpha=2.5),
     wpoints=F,add=T,lwd=1,col=colors[2])
plot(ahull(as.vector(dxsp_change_sd[3,,]),as.vector(D17Op_change_sd[3,,]),alpha=2.5),
     wpoints=F,add=T,lwd=1,col=colors[3])
plot(ahull(as.vector(dxsp_change_sd[4,,]),as.vector(D17Op_change_sd[4,,]),alpha=2.5),
     wpoints=F,add=T,lwd=1,col=colors[4])
text(0.28,14.65,"(c)",cex=0.9)

plot.new()
legend("topleft",c("f = 1","f = 0.7","f = 0.4","f = 0.2"),
       col=colors[1:4],lwd=1,cex=0.9)