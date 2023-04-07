rm(list=ls())

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

Rayleigh_computer4 <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,step,advection=T,lambda){
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
              d18Ov_list,d17Ov_list,d2Hv_list,f_list))
}

### d-excess/D17O vs f ###
result1 <- Rayleigh_computer4(-10,10,10,15,0,1000,T,0.004)
result2 <- Rayleigh_computer4(-15,10,10,20,0,1000,T,0.004)
result3 <- Rayleigh_computer4(-20,10,10,25,0,1000,T,0.004)
result4 <- Rayleigh_computer4(-15,10,10,-15,-50,1000,T,0.004)
result5 <- Rayleigh_computer4(-25,10,10,-10,-50,1000,T,0.004)
result6 <- Rayleigh_computer4(-35,10,10,-5,-50,1000,T,0.004)

par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
# 1
plot(NA,NA,xlim=c(1,0.2),ylim=c(0,30),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(1,0.2,-0.2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,30,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("f in Rayleigh distillation",1,padj=2.2,cex=0.8)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.8)

lines(result1[[7]],result1[[3]]-8*result1[[1]],type='l',col="blue")
lines(result1[[7]],result1[[6]]-8*result1[[4]],type='l',col="blue",lty=2)
lines(result2[[7]],result2[[3]]-8*result2[[1]],type='l',col="black")
lines(result2[[7]],result2[[6]]-8*result2[[4]],type='l',col="black",lty=2)
lines(result3[[7]][1:980],result3[[3]][1:980]-8*result3[[1]][1:980],type='l',col="red")
lines(result3[[7]][1:980],result3[[6]][1:980]-8*result3[[4]][1:980],type='l',col="red",lty=2)
text(0.98,31.3,"(a)",cex=0.9)
legend("topright",c(expression(paste("initial T"[d]," = 15\u00B0C, vapor ",delta^18,"O = -10\u2030")),
                   expression(paste("initial T"[d]," = 20\u00B0C, vapor ",delta^18,"O = -15\u2030")),
                   expression(paste("initial T"[d]," = 25\u00B0C, vapor ",delta^18,"O = -20\u2030")),"vapor"),
       col=c("blue","black","red","black"),lwd=1,lty=c(1,1,1,2),cex=0.9,text.width=0.58)

# 2
plot(NA,NA,xlim=c(1,0.2),ylim=c(0,30),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(1,0.2,-0.2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,30,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("f in Rayleigh distillation",1,padj=2.2,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

lines(result1[[7]],get_D17O(result1[[1]],result1[[2]]),type='l',col="blue")
lines(result1[[7]],get_D17O(result1[[4]],result1[[5]]),type='l',col="blue",lty=2)
lines(result2[[7]],get_D17O(result2[[1]],result2[[2]]),type='l',col="black")
lines(result2[[7]],get_D17O(result2[[4]],result2[[5]]),type='l',col="black",lty=2)
lines(result3[[7]][1:980],get_D17O(result3[[1]],result3[[2]])[1:980],type='l',col="red")
lines(result3[[7]][1:980],get_D17O(result3[[4]],result3[[5]])[1:980],type='l',col="red",lty=2)
text(0.98,31.3,"(b)",cex=0.9)

# 3
plot(NA,NA,xlim=c(1,0.2),ylim=c(-30,70),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(1,0.2,-0.2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,60,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("f in Rayleigh distillation",1,padj=2.2,cex=0.8)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.8)

lines(result4[[7]][1:514],result4[[3]][1:514]-8*result4[[1]][1:514],type='l',col="blue")
lines(result4[[7]][1:514],result4[[6]][1:514]-8*result4[[4]][1:514],type='l',col="blue",lty=2)
lines(result5[[7]][1:469],result5[[3]][1:469]-8*result5[[1]][1:469],type='l',col="black")
lines(result5[[7]][1:469],result5[[6]][1:469]-8*result5[[4]][1:469],type='l',col="black",lty=2)
lines(result6[[7]][1:434],result6[[3]][1:434]-8*result6[[1]][1:434],type='l',col="red")
lines(result6[[7]][1:434],result6[[6]][1:434]-8*result6[[4]][1:434],type='l',col="red",lty=2)
text(0.98,74.5,"(c)",cex=0.9)
legend("topright",c(expression(paste("initial T"[d]," = -15\u00B0C, vapor ",delta^18,"O = -15\u2030")),
                   expression(paste("initial T"[d]," = -10\u00B0C, vapor ",delta^18,"O = -25\u2030")),
                   expression(paste("initial T"[d]," = -5\u00B0C, vapor ",delta^18,"O = -35\u2030")),"vapor"),
       col=c("blue","black","red","black"),lwd=1,lty=c(1,1,1,2),cex=0.9,text.width=0.6)

# 4
plot(NA,NA,xlim=c(1,0.2),ylim=c(-30,70),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(1,0.2,-0.2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(-20,60,20),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext("f in Rayleigh distillation",1,padj=2.2,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

lines(result4[[7]][1:514],get_D17O(result4[[1]],result4[[2]])[1:514],type='l',col="blue")
lines(result4[[7]][1:393],get_D17O(result4[[4]],result4[[5]])[1:393],type='l',col="blue",lty=2)
lines(result5[[7]][1:469],get_D17O(result5[[1]],result5[[2]])[1:469],type='l',col="black")
lines(result5[[7]][1:443],get_D17O(result5[[4]],result5[[5]])[1:443],type='l',col="black",lty=2)
lines(result6[[7]][1:434],get_D17O(result6[[1]],result6[[2]])[1:434],type='l',col="red")
lines(result6[[7]][1:434],get_D17O(result6[[4]],result6[[5]])[1:434],type='l',col="red",lty=2)
text(0.98,74.5,"(d)",cex=0.9)