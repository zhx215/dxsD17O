rm(list=ls())
library(pracma)
library(raster)
library(scales)
load("C:\\Users\\zhyxi\\Desktop\\17O project\\era_resampling_monthly_v2.RData")
load("C:\\Users\\zhyxi\\Desktop\\17O project\\lmdz4_resampling.RData")
load("C:\\Users\\zhyxi\\Desktop\\17O project\\RR and rh.RData")
par(mfrow=c(3,3),mar=c(3,3.5,1.5,1),xpd=T)
colorbar <- colorRampPalette(c("blue","white","red"))
lmin <- -20
lmax <- 27

# fractionation
a18eff <- function(temp,mode="k",vmt,mit,l=0.004) {
  if (temp >= vmt) {
    eq <- exp(1137/temp^2-0.4156/temp-0.00207)
    value <- eq
  }
  else if (temp <= mit) {
    eq <- exp(11.839/temp-0.028224)
    Si <- 1-l*(temp-273.15)
    diff16_18 <- 1.0285
    k <- Si/(eq*diff16_18*(Si-1)+1)
    if (mode=="no k") {k <- 1}
    value <- k*eq
  }
  else {
    eq_273 <- exp(1137/vmt^2-0.4156/vmt-0.00207)
    eq_274 <- exp(1137/(vmt+1)^2-0.4156/(vmt+1)-0.00207)
    eq_250 <- exp(11.839/mit-0.028224)
    eq_249 <- exp(11.839/(mit-1)-0.028224)
    Si_250 <- 1-l*(mit-273.15)
    Si_249 <- 1-l*((mit-1)-273.15)
    diff16_18 <- 1.0285
    k_250 <- Si_250/(eq_250*diff16_18*(Si_250-1)+1)
    k_249 <- Si_249/(eq_249*diff16_18*(Si_249-1)+1)
    if (mode=="no k") {k_250 <- k_249 <- 1}
    x1 <- mit
    x2 <- vmt
    k1 <- eq_250*k_250-eq_249*k_249
    k2 <- eq_274-eq_273
    y1 <- eq_250*k_250
    y2 <- eq_273
    answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
    value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
  }
  return(value)
}

a17eff <- function(temp,mode="k",vmt,mit,l=0.004) {
  if (temp >= vmt) {
    eq <- exp(1137/temp^2-0.4156/temp-0.00207)^0.529
    value <- eq
  }
  else if (temp <= mit) {
    eq <- exp(11.839/temp-0.028224)^0.529
    Si <- 1-l*(temp-273.15)
    diff16_17 <- 1.0285^0.518
    k <- Si/(eq*diff16_17*(Si-1)+1)
    if (mode=="no k") {k <- 1}
    value <- k*eq
  }
  else {
    eq_273 <- exp(1137/vmt^2-0.4156/vmt-0.00207)^0.529
    eq_274 <- exp(1137/(vmt+1)^2-0.4156/(vmt+1)-0.00207)^0.529
    eq_250 <- exp(11.839/mit-0.028224)^0.529
    eq_249 <- exp(11.839/(mit-1)-0.028224)^0.529
    Si_250 <- 1-l*(mit-273.15)
    Si_249 <- 1-l*((mit-1)-273.15)
    diff16_17 <- 1.0285^0.518
    k_250 <- Si_250/(eq_250*diff16_17*(Si_250-1)+1)
    k_249 <- Si_249/(eq_249*diff16_17*(Si_249-1)+1)
    if (mode=="no k") {k_250 <- k_249 <- 1}
    x1 <- mit
    x2 <- vmt
    k1 <- eq_250*k_250-eq_249*k_249
    k2 <- eq_274-eq_273
    y1 <- eq_250*k_250
    y2 <- eq_273
    answer <- solve(matrix(c(x1^3,x2^3,3*x1^2,3*x2^2,x1^2,x2^2,2*x1,2*x2,x1,x2,1,1,1,1,0,0),ncol=4),c(y1,y2,k1,k2))
    value <- answer[1]*temp^3+answer[2]*temp^2+answer[3]*temp+answer[4]
  }
  return(value)
}

a2eff <- function(temp,mode="k",vmt,mit,l=0.004) {
  if (temp >= vmt) {
    eq <- exp(24844/temp^2-76.248/temp+0.05261)
    value <- eq
  }
  else if (temp <= mit) {
    eq <- exp(16289/temp^2-0.0945)
    Si <- 1-l*(temp-273.15)
    diff1_2 <- 1.0251
    k <- Si/(eq*diff1_2*(Si-1)+1)
    if (mode=="no k") {k <- 1}
    value <- k*eq
  }
  else {
    eq_273 <- exp(24844/vmt^2-76.248/vmt+0.05261)
    eq_274 <- exp(24844/(vmt+1)^2-76.248/(vmt+1)+0.05261)
    eq_250 <- exp(16289/mit^2-0.0945)
    eq_249 <- exp(16289/(mit-1)^2-0.0945)
    Si_250 <- 1-l*(mit-273.15)
    Si_249 <- 1-l*((mit-1)-273.15)
    diff1_2 <- 1.0251
    k_250 <- Si_250/(eq_250*diff1_2*(Si_250-1)+1)
    k_249 <- Si_249/(eq_249*diff1_2*(Si_249-1)+1)
    if (mode=="no k") {k_250 <- k_249 <- 1}
    x1 <- mit
    x2 <- vmt
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

Rayleigh_computer_1 <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,step,advection=T,vmt,mit,lambda){
  if (advection==T) {times <- 1} else (times <- 0.5)
  td_list <- seq(initial_td,final_td,length.out=step)
  f_list <- CC_relation(td_list)
  f_list <- f_list/f_list[1]
  
  d18Ov_list <- d18Ov
  d2Hv_list <- dxsv+8*d18Ov
  d17Ov_list <- get_d17O(d18Ov,D17Ov)
  
  d18Op_list <- (d18Ov_list+1000)*a18eff(initial_td+273.15,"k",vmt,mit,lambda)-1000
  d2Hp_list <- (d2Hv_list+1000)*a2eff(initial_td+273.15,"k",vmt,mit,lambda)-1000
  d17Op_list <- (d17Ov_list+1000)*a17eff(initial_td+273.15,"k",vmt,mit,lambda)-1000
  
  for (i in 2:step){
    d18Ov_list[i] <- (d18Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a18eff(td_list[i]+273.15,"k",vmt,mit,lambda)^times-1)-1000
    d17Ov_list[i] <- (d17Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a17eff(td_list[i]+273.15,"k",vmt,mit,lambda)^times-1)-1000
    d2Hv_list[i] <- (d2Hv_list[i-1]+1000)*(f_list[i]/f_list[i-1])^(a2eff(td_list[i]+273.15,"k",vmt,mit,lambda)^times-1)-1000
    
    d18Op_list[i] <- (d18Ov_list[i]+1000)*a18eff(td_list[i]+273.15,"k",vmt,mit,lambda)-1000
    d17Op_list[i] <- (d17Ov_list[i]+1000)*a17eff(td_list[i]+273.15,"k",vmt,mit,lambda)-1000
    d2Hp_list[i] <- (d2Hv_list[i]+1000)*a2eff(td_list[i]+273.15,"k",vmt,mit,lambda)-1000
  }
  
  return(c(d18Op_list[step],d17Op_list[step],d2Hp_list[step],
           d18Ov_list[step],d17Ov_list[step],d2Hv_list[step]))
}

mixing_computer <- function(d18O_input,dxs_input,D17O_input,td_input,td_final,vmt,mit,lambda,n,td_step){
  repeat{
    d17O_input <- get_d17O(d18O_input,D17O_input)
    d2H_input <- dxs_input+8*d18O_input

    if ((td_input-td_step) < td_final) {
      return(c((d18O_input+1000)*a18eff(td_input+273.15,"k",vmt,mit,lambda)-1000,
               (d17O_input+1000)*a17eff(td_input+273.15,"k",vmt,mit,lambda)-1000,
               (d2H_input+1000)*a2eff(td_input+273.15,"k",vmt,mit,lambda)-1000))}
    
    result <- Rayleigh_computer_1(d18O_input,dxs_input,D17O_input,td_input,td_input-td_step,n,T,vmt,mit,lambda)
    d18O_mix <- weighted.mean(c(d18O_input,result[4]),
                              c(CC_relation(td_input),CC_relation(td_input-td_step)))
    d17O_mix <- weighted.mean(c(d17O_input,result[5]),
                              c(CC_relation(td_input),CC_relation(td_input-td_step)))
    d2H_mix <- weighted.mean(c(d2H_input,result[6]),
                             c(CC_relation(td_input),CC_relation(td_input-td_step)))
    dxs_mix <- d2H_mix-8*d18O_mix
    D17O_mix <- get_D17O(d18O_mix,d17O_mix)
    td_mix <- weighted.mean(c(td_input,td_input-td_step),
                            c(CC_relation(td_input),CC_relation(td_input-td_step)))
    
    d18O_input <- d18O_mix
    dxs_input <- dxs_mix
    D17O_input <- D17O_mix
    td_input <- td_mix
  }
}

Rayleigh_computer_1b <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,initial_t,final_t,step,advection=T,vmt,mit,lambda,dset=F,dd=NA){
  if (advection==T) {times <- 1} else (times <- 0.5)
  td_list <- seq(initial_td,final_td,length.out=step)
  f_list <- CC_relation(td_list)
  f_list <- f_list/f_list[1]
  t_list <- seq(initial_t,final_t,length.out=step)
  if (dset==T) {d <- dd} else {d <- sample(1:12,1)}
  
  rh_list <- a18rr <- a17rr <- a2rr <- e <- vector()
  for (i in 1:step){
    rh_list[i] <- infer_rh(t_list[i],td_list[i])
    t_id <- which.min(abs(t_list[i]-seq(0,44,0.25)))
    rh_id <- which.min(abs(rh_list[i]-seq(0.4,1,0.0025)))
    e[i] <- rr_result_matrix[t_id,rh_id,d,1]
    a18rr[i] <- rr_result_matrix[t_id,rh_id,d,2]
    a17rr[i] <- rr_result_matrix[t_id,rh_id,d,3]
    a2rr[i] <- rr_result_matrix[t_id,rh_id,d,4]
  }
  if (max(e)==100){return(NA)}
  
  d18Ov_list <- d18Ov
  d2Hv_list <- dxsv+8*d18Ov
  d17Ov_list <- get_d17O(d18Ov,D17Ov)
  
  d18Op_list <- (d18Ov_list+1000)*a18eff(initial_td+273.15,"k",vmt,mit,lambda)*a18rr[1]-1000
  d2Hp_list <- (d2Hv_list+1000)*a2eff(initial_td+273.15,"k",vmt,mit,lambda)*a2rr[1]-1000
  d17Op_list <- (d17Ov_list+1000)*a17eff(initial_td+273.15,"k",vmt,mit,lambda)*a17rr[1]-1000
  
  for (i in 2:step){
    d18Ov_list[i] <- (d18Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^((a18eff(td_list[i]+273.15,"k",vmt,mit,lambda)*a18rr[i])^times-1)-1000
    d17Ov_list[i] <- (d17Ov_list[i-1]+1000)*(f_list[i]/f_list[i-1])^((a17eff(td_list[i]+273.15,"k",vmt,mit,lambda)*a17rr[i])^times-1)-1000
    d2Hv_list[i] <- (d2Hv_list[i-1]+1000)*(f_list[i]/f_list[i-1])^((a2eff(td_list[i]+273.15,"k",vmt,mit,lambda)*a2rr[i])^times-1)-1000
    
    d18Op_list[i] <- (d18Ov_list[i]+1000)*a18eff(td_list[i]+273.15,"k",vmt,mit,lambda)*a18rr[i]-1000
    d17Op_list[i] <- (d17Ov_list[i]+1000)*a17eff(td_list[i]+273.15,"k",vmt,mit,lambda)*a17rr[i]-1000
    d2Hp_list[i] <- (d2Hv_list[i]+1000)*a2eff(td_list[i]+273.15,"k",vmt,mit,lambda)*a2rr[i]-1000
  }
  
  return(c(d18Op_list[step],d17Op_list[step],d2Hp_list[step],
           d18Ov_list[step],d17Ov_list[step],d2Hv_list[step],
           d18Op_list[1],d17Op_list[1],d2Hp_list[1]))
}

# 1
# d18Op_result <- d17Op_result <- d2Hp_result <- td_sink_result <- vector()
# repeat{
#   mn <- sample(1:24,1,replace=T,rowSums(d2m_land_freq_slice))
#   rhs_source <- sample(rhs_unit,1,replace=T,prob=rhs_freq_slice[mn,])
#   td_source <- sample(d2m_ocean_unit,1,replace=T,prob=td_freq_per_rhs_slice[mn,which(rhs_source==rhs_unit),])
#   td_sink <- sample(d2m_land_unit,1,replace=T,prob=d2m_land_freq_slice[mn,])
#   if ((td_sink <= td_source)&(td_sink >= -20)){
#     if (td_source==28){
#       d18Ov_source <- sample(d18O_lmdz_unit,1,replace=T,prob=d18O_freq_per_td[74,])
#     } else {
#       d18Ov_source <- sample(d18O_lmdz_unit,1,replace=T,prob=d18O_freq_per_td[which(td_source==td_lmdz_unit),])
#     }
#     dxsv_source <- rnorm(1,-0.54*rhs_source+48.2,4.35)
#     D17Ov_source <- rnorm(1,-0.74*rhs_source+65.3,5.6)
#     result <- Rayleigh_computer_1(d18Ov_source,dxsv_source,D17Ov_source,td_source,td_sink,100,T,273.15,-23+273.15,0.004)
#     d18Op_result <- c(d18Op_result,result[1])
#     d17Op_result <- c(d17Op_result,result[2])
#     d2Hp_result <- c(d2Hp_result,result[3])
#     td_sink_result <- c(td_sink_result,td_sink)
#     if (length(d18Op_result)==100000){
#       break
#     }
#     print(length(d18Op_result))
#   }
# }
# dxsp_result <- d2Hp_result-8*d18Op_result
# D17Op_result <- get_D17O(d18Op_result,d17Op_result)
# D17Op_result_wer <- vector()
# for (i in 1:length(D17Op_result)){
#   D17Op_result_wer[i] <- rnorm(1,D17Op_result[i],5)
# }
# save.image("C:/Users/zhyxi/Desktop/17O project/SEP simulation 1.RData")

load("C:/Users/zhyxi/Desktop/17O project/SEP simulation 1.RData")
# a
plot(NA,NA,xlim=c(-25,5),ylim=c(-21,42),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-25,5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,40,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.7)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.7)

polygon(c(-25,0,5,5,0,-25),c(22,22,10,-20,-2,-2),border=NA,col=alpha("cornflowerblue",0.4))

d18O_slice <- seq(-24.5,4.5,1)
for (i in 1:length(d18O_slice)){
  id <- which((d18Op_result>(d18O_slice[i]-0.5))&(d18Op_result<=(d18O_slice[i]+0.5))&
                (dxsp_result< 40)&(dxsp_result> -20))
  mean <- mean(dxsp_result[id])
  sd <- sd(dxsp_result[id])
  lines(rep(d18O_slice[i],2),c(mean-2*sd,mean+2*sd),lwd=0.5)
  points(d18O_slice[i],mean,pch=21,bg="white",cex=1)
}

text(-26,45.5,"(a) Stochastic simulation #1",cex=1,pos=4)
text(-25,38,"Rayleigh distillation",cex=1,pos=4)

lines(c(-25,0,5),c(10,10,-5),col="cornflowerblue",lwd=1)

# b
plot(NA,NA,xlim=c(-25,5),ylim=c(-42,82),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-25,5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,80,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.7)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.7)

d18O_slice <- seq(-24.5,4.5,1)
for (i in 1:length(d18O_slice)){
  id <- which((d18Op_result>(d18O_slice[i]-0.5))&(d18Op_result<=(d18O_slice[i]+0.5))&
                (D17Op_result_wer< 80)&(D17Op_result_wer> -40))
  mean <- mean(D17Op_result_wer[id])
  sd <- sd(D17Op_result_wer[id])
  lines(rep(d18O_slice[i],2),c(mean-2*sd,mean+2*sd),lwd=0.5)
  points(d18O_slice[i],mean,pch=21,bg="white",cex=1)
}

lines(c(-25,0,5),c(38,18,-22),col="cornflowerblue",lwd=1)

# c
plot(NA,NA,xlim=c(-21,42),ylim=c(-42,82),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-20,40,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,80,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("d-excess (","\u2030",")")),1,padj=2,cex=0.7)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.7)

id <- which((dxsp_result> -20)&(dxsp_result< 40)&(D17Op_result_wer< 80)&(D17Op_result_wer> -40))
freq <- array(0,dim=c(length(seq(-21,42,3)),length(seq(-42,82,4))))
for (i in 1:length(id)){
  x <- which.min(abs(dxsp_result[id][i]-seq(-21,42,3)))
  y <- which.min(abs(D17Op_result_wer[id][i]-seq(-42,82,4)))
  freq[x,y] <- freq[x,y] + 1
}
contour(seq(-21,42,3),seq(-42,82,4),freq,levels=c(0.005/100*sum(freq),0.05/100*sum(freq),0.5/100*sum(freq),2/100*sum(freq)),labels=c("0.005%","0.05%","0.5%","2%"),vfont=NULL,add=T)
coef <- odregress(dxsp_result[id],D17Op_result_wer[id])$coef
lines((c(-42,82)-coef[2])/coef[1],c(-42,82))
dxs_output1 <- dxsp_result[id]
D17O_output1 <- D17Op_result_wer[id]

# 2
# d18Op_result <- d17Op_result <- d2Hp_result <- td_sink_result <- vector()
# repeat{
#   mn <- sample(1:24,1,replace=T,rowSums(d2m_land_freq_slice))
#   rhs_source <- sample(rhs_unit,1,replace=T,prob=rhs_freq_slice[mn,])
#   td_source <- sample(d2m_ocean_unit,1,replace=T,prob=td_freq_per_rhs_slice[mn,which(rhs_source==rhs_unit),])
#   td_sink <- sample(d2m_land_unit,1,replace=T,prob=d2m_land_freq_slice[mn,])
#   if ((td_sink <= td_source)&(td_sink >= -20)){
#     if (td_source==28){
#       d18Ov_source <- sample(d18O_lmdz_unit,1,replace=T,prob=d18O_freq_per_td[74,])
#     } else {
#       d18Ov_source <- sample(d18O_lmdz_unit,1,replace=T,prob=d18O_freq_per_td[which(td_source==td_lmdz_unit),])
#     }
#     dxsv_source <- rnorm(1,-0.54*rhs_source+48.2,4.35)
#     D17Ov_source <- rnorm(1,-0.74*rhs_source+65.3,5.6)
#     t_source <- sample(t2m_ocean_unit,1,replace=T,prob=t2m_freq_per_d2m_ocean_slice[mn,which(td_source==d2m_ocean_unit),])
#     t_sink <- sample(t2m_land_unit,1,replace=T,prob=t2m_freq_per_d2m_land_slice[mn,which(td_sink==d2m_land_unit),])
#     result <- Rayleigh_computer_1b(d18Ov_source,dxsv_source,D17Ov_source,td_source,td_sink,t_source,t_sink,100,T,273.15,-23+273.15,0.004)
#     if (length(result)==1){
#       print("skipped")
#       next
#     }
#     d18Op_result <- c(d18Op_result,result[1])
#     d17Op_result <- c(d17Op_result,result[2])
#     d2Hp_result <- c(d2Hp_result,result[3])
#     td_sink_result <- c(td_sink_result,td_sink)
#     if (length(d18Op_result)==100000){
#       break
#     }
#     print(length(d18Op_result))
#   }
# }
# dxsp_result <- d2Hp_result-8*d18Op_result
# D17Op_result <- get_D17O(d18Op_result,d17Op_result)
# D17Op_result_wer <- vector()
# for (i in 1:length(D17Op_result)){
#   D17Op_result_wer[i] <- rnorm(1,D17Op_result[i],5)
# }
# save.image("C:/Users/zhyxi/Desktop/17O project/SEP simulation 3.RData")

load("C:/Users/zhyxi/Desktop/17O project/SEP simulation 3.RData")
# a
plot(NA,NA,xlim=c(-25,5),ylim=c(-21,42),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-25,5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,40,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.7)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.7)

polygon(c(-25,0,5,5,0,-25),c(22,22,10,-20,-2,-2),border=NA,col=alpha("cornflowerblue",0.4))

d18O_slice <- seq(-24.5,4.5,1)
for (i in 1:length(d18O_slice)){
  id <- which((d18Op_result>(d18O_slice[i]-0.5))&(d18Op_result<=(d18O_slice[i]+0.5))&
                (dxsp_result< 40)&(dxsp_result> -20))
  mean <- mean(dxsp_result[id])
  sd <- sd(dxsp_result[id])
  lines(rep(d18O_slice[i],2),c(mean-2*sd,mean+2*sd),lwd=0.5)
  points(d18O_slice[i],mean,pch=21,bg="white",cex=1)
}

text(-26,45.5,"(b) Stochastic simulation #2",cex=1,pos=4)
text(-25,35,"Rayleigh distillation +\nraindrop re-evaporation",cex=1,pos=4)

lines(c(-25,0,5),c(10,10,-5),col="cornflowerblue",lwd=1)

# b
plot(NA,NA,xlim=c(-25,5),ylim=c(-42,82),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-25,5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,80,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.7)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.7)

d18O_slice <- seq(-24.5,4.5,1)
for (i in 1:length(d18O_slice)){
  id <- which((d18Op_result>(d18O_slice[i]-0.5))&(d18Op_result<=(d18O_slice[i]+0.5))&
                (D17Op_result_wer< 80)&(D17Op_result_wer> -40))
  mean <- mean(D17Op_result_wer[id])
  sd <- sd(D17Op_result_wer[id])
  lines(rep(d18O_slice[i],2),c(mean-2*sd,mean+2*sd),lwd=0.5)
  points(d18O_slice[i],mean,pch=21,bg="white",cex=1)
}

lines(c(-25,0,5),c(38,18,-22),col="cornflowerblue",lwd=1)

# c
plot(NA,NA,xlim=c(-21,42),ylim=c(-42,82),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-20,40,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,80,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("d-excess (","\u2030",")")),1,padj=2,cex=0.7)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.7)

id <- which((dxsp_result> -20)&(dxsp_result< 40)&(D17Op_result_wer< 80)&(D17Op_result_wer> -40))
freq <- array(0,dim=c(length(seq(-21,42,3)),length(seq(-42,82,4))))
for (i in 1:length(id)){
  x <- which.min(abs(dxsp_result[id][i]-seq(-21,42,3)))
  y <- which.min(abs(D17Op_result_wer[id][i]-seq(-42,82,4)))
  freq[x,y] <- freq[x,y] + 1
}
contour(seq(-21,42,3),seq(-42,82,4),freq,levels=c(0.005/100*sum(freq),0.05/100*sum(freq),0.5/100*sum(freq),2/100*sum(freq)),labels=c("0.005%","0.05%","0.5%","2%"),vfont=NULL,add=T)
coef <- odregress(dxsp_result[id],D17Op_result_wer[id])$coef
lines((c(-42,82)-coef[2])/coef[1],c(-42,82))
dxs_output2 <- dxsp_result[id]
D17O_output2 <- D17Op_result_wer[id]

# 4
# d18Op_result <- d17Op_result <- d2Hp_result <- td_sink_result <- vector()
# repeat{
#   mn <- sample(1:24,1,replace=T,rowSums(d2m_land_freq_slice))
#   rhs_source <- sample(rhs_unit,1,replace=T,prob=rhs_freq_slice[mn,])
#   td_source <- sample(d2m_ocean_unit,1,replace=T,prob=td_freq_per_rhs_slice[mn,which(rhs_source==rhs_unit),])
#   td_sink <- sample(d2m_land_unit,1,replace=T,prob=d2m_land_freq_slice[mn,])
#   if ((td_sink <= td_source)&(td_sink >= -20)){
#     if (td_source==28){
#       d18Ov_source <- sample(d18O_lmdz_unit,1,replace=T,prob=d18O_freq_per_td[74,])
#     } else {
#       d18Ov_source <- sample(d18O_lmdz_unit,1,replace=T,prob=d18O_freq_per_td[which(td_source==td_lmdz_unit),])
#     }
#     dxsv_source <- rnorm(1,-0.54*rhs_source+48.2,4.35)
#     D17Ov_source <- rnorm(1,-0.74*rhs_source+65.3,5.6)
#     t_source <- sample(t2m_ocean_unit,1,replace=T,prob=t2m_freq_per_d2m_ocean_slice[mn,which(td_source==d2m_ocean_unit),])
#     t_sink <- sample(t2m_land_unit,1,replace=T,prob=t2m_freq_per_d2m_land_slice[mn,which(td_sink==d2m_land_unit),])
#     result <- mixing_computer_1b(d18Ov_source,dxsv_source,D17Ov_source,td_source,td_sink,t_source,t_sink,273.15,-23+273.15,0.004,100,6)
#     if (length(result)==1){
#       print("skipped")
#       next
#     }
#     d18Op_result <- c(d18Op_result,result[1])
#     d17Op_result <- c(d17Op_result,result[2])
#     d2Hp_result <- c(d2Hp_result,result[3])
#     td_sink_result <- c(td_sink_result,td_sink)
#     if (length(d18Op_result)==100000){
#       break
#     }
#     print(length(d18Op_result))
#   }
# }
# dxsp_result <- d2Hp_result-8*d18Op_result
# D17Op_result <- get_D17O(d18Op_result,d17Op_result)
# D17Op_result_wer <- vector()
# for (i in 1:length(D17Op_result)){
#   D17Op_result_wer[i] <- rnorm(1,D17Op_result[i],5)
# }
# save.image("C:/Users/zhyxi/Desktop/17O project/SEP simulation 4.RData")

load("C:/Users/zhyxi/Desktop/17O project/SEP simulation 4.RData")
# a
plot(NA,NA,xlim=c(-25,5),ylim=c(-21,42),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-25,5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,40,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.7)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.7)

polygon(c(-25,0,5,5,0,-25),c(22,22,10,-20,-2,-2),border=NA,col=alpha("cornflowerblue",0.4))

d18O_slice <- seq(-24.5,4.5,1)
for (i in 1:length(d18O_slice)){
  id <- which((d18Op_result>(d18O_slice[i]-0.5))&(d18Op_result<=(d18O_slice[i]+0.5))&
                (dxsp_result< 40)&(dxsp_result> -20))
  mean <- mean(dxsp_result[id])
  sd <- sd(dxsp_result[id])
  lines(rep(d18O_slice[i],2),c(mean-2*sd,mean+2*sd),lwd=0.5)
  points(d18O_slice[i],mean,pch=21,bg="white",cex=1)
}

text(-26,45.5,"(c) Stochastic simulation #3",cex=1,pos=4)
text(-25,32,"Rayleigh distillation +\nraindrop re-evaporation +\nvapor mixing",cex=1,pos=4)

lines(c(-25,0,5),c(10,10,-5),col="cornflowerblue",lwd=1)

# b
plot(NA,NA,xlim=c(-25,5),ylim=c(-42,82),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-25,5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,80,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.7)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.7)

d18O_slice <- seq(-24.5,4.5,1)
for (i in 1:length(d18O_slice)){
  id <- which((d18Op_result>(d18O_slice[i]-0.5))&(d18Op_result<=(d18O_slice[i]+0.5))&
                (D17Op_result_wer< 80)&(D17Op_result_wer> -40))
  mean <- mean(D17Op_result_wer[id])
  sd <- sd(D17Op_result_wer[id])
  lines(rep(d18O_slice[i],2),c(mean-2*sd,mean+2*sd),lwd=0.5)
  points(d18O_slice[i],mean,pch=21,bg="white",cex=1)
}

lines(c(-25,0,5),c(38,18,-22),col="cornflowerblue",lwd=1)

# c
plot(NA,NA,xlim=c(-21,42),ylim=c(-42,82),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,seq(-20,40,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-40,80,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("d-excess (","\u2030",")")),1,padj=2,cex=0.7)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.7)

id <- which((dxsp_result> -20)&(dxsp_result< 40)&(D17Op_result_wer< 80)&(D17Op_result_wer> -40))
freq <- array(0,dim=c(length(seq(-21,42,3)),length(seq(-42,82,4))))
for (i in 1:length(id)){
  x <- which.min(abs(dxsp_result[id][i]-seq(-21,42,3)))
  y <- which.min(abs(D17Op_result_wer[id][i]-seq(-42,82,4)))
  freq[x,y] <- freq[x,y] + 1
}
contour(seq(-21,42,3),seq(-42,82,4),freq,levels=c(0.005/100*sum(freq),0.05/100*sum(freq),0.5/100*sum(freq),2/100*sum(freq)),labels=c("0.005%","0.05%","0.5%","2%"),vfont=NULL,add=T)
coef <- odregress(dxsp_result[id],D17Op_result_wer[id])$coef
lines((c(-42,82)-coef[2])/coef[1],c(-42,82))
dxs_output3 <- dxsp_result[id]
D17O_output3 <- D17Op_result_wer[id]
# 
# rm(list=setdiff(ls(), c("dxs_output1","dxs_output2","dxs_output3",
#                         "D17O_output1","D17O_output2","D17O_output3")))
# save.image("C:/Users/zhyxi/Desktop/17O project/random1_output.RData")