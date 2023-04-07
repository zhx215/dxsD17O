rm(list=ls())
load("C:\\Users\\zhyxi\\Desktop\\new\\RR and rh.RData")

#
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

# ae
a18eeq <- function(temp){
  if (temp>273.15){
    eq <- exp(1137/temp^2-0.4156/temp-0.00207)
    return(1/eq)
  } else {
    eq <- exp(1137/273.15^2-0.4156/273.15-0.00207)
    return(1/eq)
  }
}

a17eeq <- function(temp){
  if (temp>273.15){
    eq <- exp(1137/temp^2-0.4156/temp-0.00207)^0.529
    return(1/eq)
  } else {
    eq <- exp(1137/273.15^2-0.4156/273.15-0.00207)^0.529
    return(1/eq)
  }
}

a2eeq <- function(temp){
  if (temp>273.15){
    eq <- exp(24844/temp^2-76.248/temp+0.05261)
    return(1/eq)
  } else {
    eq <- exp(24844/273.15^2-76.248/273.15+0.05261)
    return(1/eq)
  }
}

rayleigh_ET <- function(nd,t,h,delta_initial,alpha,k,smow){
  frac <-  1/(1/nd+1)
  r_initial <- (delta_initial/1000+1)*smow
  ae <- ((1-t)*((alpha*k*r_initial)/(1-h+(1-t)*k*h))+t*r_initial*(1/(1+(1-t)*k*(h/(1-h)))))/r_initial
  delta_ET <- ((delta_initial+1000)*(1-frac)^ae-1000*(1-frac)-delta_initial)/(1-frac-1)*(1-t)+delta_initial*t
  return(c((delta_initial-frac*delta_ET)/(1-frac),delta_ET))
}

smow18 <- 0.0020052
smow17 <- 0.00037990
smow2 <- 0.00015576

RTM_computer_rs <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,initial_t,final_t,tet_number,ae_number,nd_number,step,advection=T,vmt,mit,lambda,d18Ors,d17Ors,d2Hrs,ps){
  if (advection==T) {times <- 1} else (times <- 0.5)
  td_list <- seq(initial_td,final_td,length.out=step)
  f_list <- CC_relation(td_list)
  f_list <- f_list/f_list[1]
  t_list <- seq(initial_t,final_t,length.out=step)
  if (nd_number>0) 
    {nd_number_prime <- nd_number/(nd_number+1)*ps} else
    {nd_number_prime <- 1/((nd_number+1)/(nd_number*ps)+(nd_number+1)/nd_number-1)}
  nd_ET <- rep(nd_number_prime,length.out=step)
  nd <- rep(nd_number,length.out=step)
  nd[t_list< 0] <- 0.001
  nd_ET[t_list< 0] <- 0.001
  
  n <- ae_number
  k18 <- 0.9723^n
  k17 <- (0.9723^0.518)^n
  k2 <- 0.9755^n
  tet <- tet_number
  
  rh_list <- a18 <- a17 <- a2 <- a18e <- a17e <- a2e <- vector()
  for (i in 1:step){
    rh_list[i] <- infer_rh(t_list[i],td_list[i])
    a18[i] <- a18eff(td_list[i]+273.15,"k",vmt,mit,lambda)
    a17[i] <- a17eff(td_list[i]+273.15,"k",vmt,mit,lambda)
    a2[i] <- a2eff(td_list[i]+273.15,"k",vmt,mit,lambda)
    a18e[i] <- a18eeq(t_list[i]+273.15)
    a17e[i] <- a17eeq(t_list[i]+273.15)
    a2e[i] <- a2eeq(t_list[i]+273.15)
  }

  d18Ov_list <- d18Ov
  d2Hv_list <- dxsv+8*d18Ov
  d17Ov_list <- get_d17O(d18Ov,D17Ov)
  
  d18Op_list <- (d18Ov_list+1000)*a18[1]-1000
  d2Hp_list <- (d2Hv_list+1000)*a2[1]-1000
  d17Op_list <- (d17Ov_list+1000)*a17[1]-1000
  
  expo18_list <- (a18+nd*a18)^times-1
  expo17_list <- (a17+nd*a17)^times-1
  expo2_list <- (a2+nd*a2)^times-1
  
  ET18_model <- rayleigh_ET(nd_ET[1],tet,rh_list[1],d18Op_list[1]/(1+1/ps)+d18Ors[1]/(1+ps),a18e[1],k18,smow18)
  ET17_model <- rayleigh_ET(nd_ET[1],tet,rh_list[1],d17Op_list[1]/(1+1/ps)+d17Ors[1]/(1+ps),a17e[1],k17,smow17)
  ET2_model <- rayleigh_ET(nd_ET[1],tet,rh_list[1],d2Hp_list[1]/(1+1/ps)+d2Hrs[1]/(1+ps),a2e[1],k2,smow2)
  d18Os_list <- ET18_model[1]; d18Oe_list <- ET18_model[2]
  d17Os_list <- ET17_model[1]; d17Oe_list <- ET17_model[2]
  d2Hs_list <- ET2_model[1]; d2He_list <- ET2_model[2]
  
  d18inf_list <- (nd[1]*d18Oe_list[1]-(1+nd[1])*(a18[1]-1)*1000)/(a18[1]+nd[1]*a18[1]-1)
  d17inf_list <- (nd[1]*d17Oe_list[1]-(1+nd[1])*(a17[1]-1)*1000)/(a17[1]+nd[1]*a17[1]-1)
  d2inf_list <- (nd[1]*d2He_list[1]-(1+nd[1])*(a2[1]-1)*1000)/(a2[1]+nd[1]*a2[1]-1)
  
  for (i in 2:step){
    d18Ov_list[i] <- (d18Ov_list[i-1]-d18inf_list[i-1])*((f_list[i]/f_list[i-1])^expo18_list[i-1])+d18inf_list[i-1]
    d17Ov_list[i] <- (d17Ov_list[i-1]-d17inf_list[i-1])*((f_list[i]/f_list[i-1])^expo17_list[i-1])+d17inf_list[i-1]
    d2Hv_list[i] <- (d2Hv_list[i-1]-d2inf_list[i-1])*((f_list[i]/f_list[i-1])^expo2_list[i-1])+d2inf_list[i-1]
    
    d18Op_list[i] <- (d18Ov_list[i]+1000)*a18[i]-1000
    d17Op_list[i] <- (d17Ov_list[i]+1000)*a17[i]-1000
    d2Hp_list[i] <- (d2Hv_list[i]+1000)*a2[i]-1000
    
    ET18_model <- rayleigh_ET(nd_ET[i],tet,rh_list[i],d18Op_list[i]/(1+1/ps)+d18Ors[i]/(1+ps),a18e[i],k18,smow18)
    ET17_model <- rayleigh_ET(nd_ET[i],tet,rh_list[i],d17Op_list[i]/(1+1/ps)+d17Ors[i]/(1+ps),a17e[i],k17,smow17)
    ET2_model <- rayleigh_ET(nd_ET[i],tet,rh_list[i],d2Hp_list[i]/(1+1/ps)+d2Hrs[i]/(1+ps),a2e[i],k2,smow2)
    d18Os_list[i] <- ET18_model[1]; d18Oe_list[i] <- ET18_model[2]
    d17Os_list[i] <- ET17_model[1]; d17Oe_list[i] <- ET17_model[2]
    d2Hs_list[i] <- ET2_model[1]; d2He_list[i] <- ET2_model[2]
    
    d18inf_list[i] <- (nd[i]*d18Oe_list[i]-(1+nd[i])*(a18[i]-1)*1000)/(a18[i]+nd[i]*a18[i]-1)
    d17inf_list[i] <- (nd[i]*d17Oe_list[i]-(1+nd[i])*(a17[i]-1)*1000)/(a17[i]+nd[i]*a17[i]-1)
    d2inf_list[i] <- (nd[i]*d2He_list[i]-(1+nd[i])*(a2[i]-1)*1000)/(a2[i]+nd[i]*a2[i]-1)
  }
  
  return(list(d18Op_list,d17Op_list,d2Hp_list,
              d18Os_list,d17Os_list,d2Hs_list,
              d18Oe_list,d17Oe_list,d2He_list))
}

repeat_RTM <- function(d18Ov,dxsv,D17Ov,initial_td,final_td,initial_t,final_t,tet_number,ae_number,nd_number,step,advection=T,vmt,mit,lambda,ps){
  d18Ors_0 <- d17Ors_0 <- d2Hrs_0 <- rep(0,step)
  result <- RTM_computer_rs(d18Ov,dxsv,D17Ov,initial_td,final_td,initial_t,final_t,tet_number,ae_number,nd_number,step,advection,vmt,mit,lambda,d18Ors_0,d17Ors_0,d2Hrs_0,ps)
  dxs_final <- (result[[3]]-8*result[[1]])[step]
  D17O_final <- (get_D17O(result[[1]],result[[2]]))[step]
  repeat{
    result <- RTM_computer_rs(d18Ov,dxsv,D17Ov,initial_td,final_td,initial_t,final_t,tet_number,ae_number,nd_number,step,advection,vmt,mit,lambda,result[[4]],result[[5]],result[[6]],ps)
    if ((abs((result[[3]]-8*result[[1]])[step]-dxs_final)<0.01)&(abs(get_D17O(result[[1]],result[[2]])[step]-D17O_final)<0.01)){
      print("done")
      return(result)
    } else {
      print((c(dxs_final,(result[[3]]-8*result[[1]])[step])))
      dxs_final <- (result[[3]]-8*result[[1]])[step]
      D17O_final <- get_D17O(result[[1]],result[[2]])[step]
    }
  }
}

#
d18Ors_0 <- d17Ors_0 <- d2Hrs_0 <- rep(0,100)
nd_number <- c(0.3,0.6,0.9)/(1-c(0.3,0.6,0.9))
tet_number <- 1-c(0.8,0.5,0.2)
ae_number <- c(0.5,0.75,1)
ps_number <- 1/0.2
Rayleigh_result <- Rayleigh_computer4b(-11,10,10,22,10,100,advection=T,0.004)

#
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
color=c("#0072B2","#009E73","#E69F00")
code <- c(24,22,25)
linecolor <- c("darkgray","black","brown3")

# 1
plot(NA,NA,xlim=c(-1,-12),ylim=c(9.5,20.5),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-12,-2,2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(10,20,2),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.8)

for (i in 1:3){
  for (j in 1:3){
    for (k in 1:3){
      RTM_result <- repeat_RTM(-11,10,10,22,10,29.5,16.7,tet_number[j],ae_number[k],nd_number[i],100,T,273.15,250.15,0.004,ps_number)
      points(RTM_result[[1]][100],RTM_result[[3]][100]-8*RTM_result[[1]][100],bg=color[i],col=linecolor[j],pch=code[k],cex=1.5,lwd=1.5)
    }
  }
}
points(Rayleigh_result[[1]][100],Rayleigh_result[[3]][100]-8*Rayleigh_result[[1]][100],pch=8,cex=1.5,lwd=1.5)
lines(Rayleigh_result[[1]],Rayleigh_result[[3]]-8*Rayleigh_result[[1]])
text(-6,10.5,"Rayleigh curve",cex=0.9,srt=-5)
text(-1.27,21,"(a)",cex=0.9)
text(-1.7,9.9,"start",cex=0.9)

legend("topleft",c("ET/P = 0.3","ET/P = 0.6","ET/P = 0.9","T/ET = 20%","T/ET = 50%","T/ET = 80%","m = 0.5","m = 0.75","m = 1","Rayleigh"),
       pch=c(rep(22,6),code,8),col=c(rep("black",3),linecolor,rep("black",4)),
       pt.lwd=1.5,cex=0.9,pt.cex=1.5,pt.bg=c(color,rep("white",3),rep("gray",3),"black"))

# 2
plot(NA,NA,xlim=c(-1,-12),ylim=c(18,57),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-12,-2,2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(20,55,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

for (i in 1:3){
  for (j in 1:3){
    for (k in 1:3){
      RTM_result <- repeat_RTM(-11,10,10,22,10,29.5,16.7,tet_number[j],ae_number[k],nd_number[i],100,T,273.15,250.15,0.004,ps_number)
      points(RTM_result[[1]][100],get_D17O(RTM_result[[1]],RTM_result[[2]])[100],bg=color[i],col=linecolor[j],pch=code[k],cex=1.5,lwd=1.5)
    }
  }
}
points(Rayleigh_result[[1]][100],get_D17O(Rayleigh_result[[1]],Rayleigh_result[[2]])[100],pch=8,cex=1.5,lwd=1.5)
lines(Rayleigh_result[[1]],get_D17O(Rayleigh_result[[1]],Rayleigh_result[[2]]))
text(-1.27,58.7,"(b)",cex=0.9)

# 3
plot(NA,NA,xlim=c(9.5,20.5),ylim=c(21,57),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(10,20,2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(25,55,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("d-excess (","\u2030",")")),1,padj=1.9,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

for (i in 1:3){
  for (j in 1:3){
    for (k in 1:3){
      RTM_result <- repeat_RTM(-11,10,10,22,10,29.5,16.7,tet_number[j],ae_number[k],nd_number[i],100,T,273.15,250.15,0.004,ps_number)
      points(RTM_result[[3]][100]-8*RTM_result[[1]][100],get_D17O(RTM_result[[1]],RTM_result[[2]])[100],bg=color[i],col=linecolor[j],pch=code[k],cex=1.5,lwd=1.5)
    }
  }
}
points(Rayleigh_result[[3]][100]-8*Rayleigh_result[[1]][100],
       get_D17O(Rayleigh_result[[1]],Rayleigh_result[[2]])[100],pch=8,cex=1.5,lwd=1.5)
x0 <- Rayleigh_result[[3]][100]-8*Rayleigh_result[[1]][100]
y0 <- get_D17O(Rayleigh_result[[1]],Rayleigh_result[[2]])[100]

lines(c(x0,20.5),1*c(x0,20.5)+y0-1*x0,lty=2)
lines(c(x0,20.5),2*c(x0,20.5)+y0-2*x0,lty=2)
lines(c(x0,20.5),3*c(x0,20.5)+y0-3*x0,lty=2)
lines(c(x0,(57-y0+4*x0)/4),4*c(x0,(57-y0+4*x0)/4)+y0-4*x0,lty=2)
text(20,31,"1:1",cex=0.9,srt=10)
text(20,41,"2:1",cex=0.9,srt=27)
text(20,51.1,"3:1",cex=0.9,srt=38)
text(18.25,54,"4:1",cex=0.9,srt=48)
text(9.77,58.7,"(c)",cex=0.9)

# 4
plot(NA,NA,xlim=c(0,30),ylim=c(-0.2,4.2),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(0,4,0.5),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      RTM_result <- repeat_RTM(-11,10,10,22,10,29.5,16.7,tet_number[j],ae_number[k],nd_number[i],100,T,273.15,250.15,0.004,ps_number)
      points(number,(get_D17O(RTM_result[[1]],RTM_result[[2]])[100]-y0)/
               (RTM_result[[3]][100]-8*RTM_result[[1]][100]-x0),
             bg=color[i],col=linecolor[j],pch=code[k],cex=1.5,lwd=1.5)
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(-0.2,4.2))
lines(c(20,20),c(-0.2,4.2))
text(5,-0.04,"T/ET = 20%",cex=0.9)
text(15,-0.04,"T/ET = 50%",cex=0.9)
text(25,-0.04,"T/ET = 80%",cex=0.9)

text(-1,4.4,expression(paste("(d) ",Delta,"(",Delta*minute^17,"O)/",
                              Delta,"(d-excess)")),cex=0.9,pos=4)