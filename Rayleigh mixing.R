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

# vapor input, T > 0C
d18O_input <- -12
dxs_input <- 10
D17O_input <- 10
d17O_input <- get_d17O(d18O_input,D17O_input)
d2H_input <- dxs_input+8*d18O_input
t_input <- 25

### case 1 ###
d18O_mix <- d18O_input
dxs_mix <- dxs_input
D17O_mix <- D17O_input
d17O_mix <- get_d17O(d18O_mix,D17O_mix)
d2H_mix <- dxs_mix+8*d18O_mix
temp_mix <- t_input
d18O_mix_p <- (d18O_mix+1000)*a18eff(273.15+temp_mix)-1000
d17O_mix_p <- (d17O_mix+1000)*a17eff(273.15+temp_mix)-1000
d2H_mix_p <- (d2H_mix+1000)*a2eff(273.15+temp_mix)-1000
dxs_mix_p <- d2H_mix_p-8*d18O_mix_p
D17O_mix_p <- get_D17O(d18O_mix_p,d17O_mix_p)
n <- 100
t_step <- 2
i <- 1
repeat{
  result <- Rayleigh_computer4(d18O_mix[i],dxs_mix[i],D17O_mix[i],temp_mix[i],temp_mix[i]-t_step,n,T,0.004)

  d18O_mix[i+1] <- weighted.mean(c(d18O_mix[i],result[[4]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d17O_mix[i+1] <- weighted.mean(c(d17O_mix[i],result[[5]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d2H_mix[i+1] <- weighted.mean(c(d2H_mix[i],result[[6]][n]),
                                c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  dxs_mix[i+1] <- d2H_mix[i+1]-8*d18O_mix[i+1]
  D17O_mix[i+1] <- get_D17O(d18O_mix[i+1],d17O_mix[i+1])
  temp_mix[i+1] <- weighted.mean(c(temp_mix[i],temp_mix[i]-t_step),
                             c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  
  d18O_mix_p[i+1] <- (d18O_mix[i+1]+1000)*a18eff(273.15+temp_mix[i+1])-1000
  d17O_mix_p[i+1] <- (d17O_mix[i+1]+1000)*a17eff(273.15+temp_mix[i+1])-1000
  d2H_mix_p[i+1] <- (d2H_mix[i+1]+1000)*a2eff(273.15+temp_mix[i+1])-1000
  dxs_mix_p[i+1] <- d2H_mix_p[i+1]-8*d18O_mix_p[i+1]
  D17O_mix_p[i+1] <- get_D17O(d18O_mix_p[i+1],d17O_mix_p[i+1])
  
  if ((temp_mix[i]-t_step)<0) {break}
  i <- i+1
}
d18O_1 <- d18O_mix
D17O_1 <- D17O_mix
dxs_1 <- dxs_mix
d18O_1_p <- d18O_mix_p
D17O_1_p <- D17O_mix_p
dxs_1_p <- dxs_mix_p
#

### case 2 ###
d18O_mix <- d18O_input
dxs_mix <- dxs_input
D17O_mix <- D17O_input
d17O_mix <- get_d17O(d18O_mix,D17O_mix)
d2H_mix <- dxs_mix+8*d18O_mix
temp_mix <- t_input
d18O_mix_p <- (d18O_mix+1000)*a18eff(273.15+temp_mix)-1000
d17O_mix_p <- (d17O_mix+1000)*a17eff(273.15+temp_mix)-1000
d2H_mix_p <- (d2H_mix+1000)*a2eff(273.15+temp_mix)-1000
dxs_mix_p <- d2H_mix_p-8*d18O_mix_p
D17O_mix_p <- get_D17O(d18O_mix_p,d17O_mix_p)
n <- 100
t_step <- 6
i <- 1
repeat{
  result <- Rayleigh_computer4(d18O_mix[i],dxs_mix[i],D17O_mix[i],temp_mix[i],temp_mix[i]-t_step,n,T,0.004)
  
  d18O_mix[i+1] <- weighted.mean(c(d18O_mix[i],result[[4]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d17O_mix[i+1] <- weighted.mean(c(d17O_mix[i],result[[5]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d2H_mix[i+1] <- weighted.mean(c(d2H_mix[i],result[[6]][n]),
                                c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  dxs_mix[i+1] <- d2H_mix[i+1]-8*d18O_mix[i+1]
  D17O_mix[i+1] <- get_D17O(d18O_mix[i+1],d17O_mix[i+1])
  temp_mix[i+1] <- weighted.mean(c(temp_mix[i],temp_mix[i]-t_step),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  
  d18O_mix_p[i+1] <- (d18O_mix[i+1]+1000)*a18eff(273.15+temp_mix[i+1])-1000
  d17O_mix_p[i+1] <- (d17O_mix[i+1]+1000)*a17eff(273.15+temp_mix[i+1])-1000
  d2H_mix_p[i+1] <- (d2H_mix[i+1]+1000)*a2eff(273.15+temp_mix[i+1])-1000
  dxs_mix_p[i+1] <- d2H_mix_p[i+1]-8*d18O_mix_p[i+1]
  D17O_mix_p[i+1] <- get_D17O(d18O_mix_p[i+1],d17O_mix_p[i+1])
  
  if ((temp_mix[i]-t_step)<0) {break}
  i <- i+1
}
d18O_2 <- d18O_mix
D17O_2 <- D17O_mix
dxs_2 <- dxs_mix
d18O_2_p <- d18O_mix_p
D17O_2_p <- D17O_mix_p
dxs_2_p <- dxs_mix_p
#

### case 3 ###
d18O_mix <- d18O_input
dxs_mix <- dxs_input
D17O_mix <- D17O_input
d17O_mix <- get_d17O(d18O_mix,D17O_mix)
d2H_mix <- dxs_mix+8*d18O_mix
temp_mix <- t_input
d18O_mix_p <- (d18O_mix+1000)*a18eff(273.15+temp_mix)-1000
d17O_mix_p <- (d17O_mix+1000)*a17eff(273.15+temp_mix)-1000
d2H_mix_p <- (d2H_mix+1000)*a2eff(273.15+temp_mix)-1000
dxs_mix_p <- d2H_mix_p-8*d18O_mix_p
D17O_mix_p <- get_D17O(d18O_mix_p,d17O_mix_p)
n <- 100
t_step <- 10
i <- 1
repeat{
  result <- Rayleigh_computer4(d18O_mix[i],dxs_mix[i],D17O_mix[i],temp_mix[i],temp_mix[i]-t_step,n,T,0.004)
  
  d18O_mix[i+1] <- weighted.mean(c(d18O_mix[i],result[[4]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d17O_mix[i+1] <- weighted.mean(c(d17O_mix[i],result[[5]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d2H_mix[i+1] <- weighted.mean(c(d2H_mix[i],result[[6]][n]),
                                c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  dxs_mix[i+1] <- d2H_mix[i+1]-8*d18O_mix[i+1]
  D17O_mix[i+1] <- get_D17O(d18O_mix[i+1],d17O_mix[i+1])
  temp_mix[i+1] <- weighted.mean(c(temp_mix[i],temp_mix[i]-t_step),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  
  d18O_mix_p[i+1] <- (d18O_mix[i+1]+1000)*a18eff(273.15+temp_mix[i+1])-1000
  d17O_mix_p[i+1] <- (d17O_mix[i+1]+1000)*a17eff(273.15+temp_mix[i+1])-1000
  d2H_mix_p[i+1] <- (d2H_mix[i+1]+1000)*a2eff(273.15+temp_mix[i+1])-1000
  dxs_mix_p[i+1] <- d2H_mix_p[i+1]-8*d18O_mix_p[i+1]
  D17O_mix_p[i+1] <- get_D17O(d18O_mix_p[i+1],d17O_mix_p[i+1])
  
  if ((temp_mix[i]-t_step)<0) {break}
  i <- i+1
}
d18O_3 <- d18O_mix
D17O_3 <- D17O_mix
dxs_3 <- dxs_mix
d18O_3_p <- d18O_mix_p
D17O_3_p <- D17O_mix_p
dxs_3_p <- dxs_mix_p
#

### Rayleigh ###
rayleigh <- Rayleigh_computer4(d18O_input,dxs_input,D17O_input,t_input,0,n,T,0.004)

#
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
#
plot(NA,NA,xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(-2,-18),
     ylim=c(7,11))
axis(1,seq(-18,-2,2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(7,11,1),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.8)

lines(rayleigh[[1]],rayleigh[[3]]-rayleigh[[1]]*8)
points(d18O_1_p[2:length(d18O_1_p)],dxs_1_p[2:length(dxs_1_p)],pch=21,bg="gray")
points(d18O_2_p[2:length(d18O_2_p)],dxs_2_p[2:length(dxs_2_p)],pch=22,bg="white")
points(d18O_3_p[2:length(d18O_3_p)],dxs_3_p[2:length(dxs_3_p)],pch=25,bg="brown1")
legend("bottomleft",c("Rayleigh curve","cooling step = 2\u00B0C","cooling step = 6\u00B0C","cooling step = 10\u00B0C"),
       pt.bg=c(NA,"gray","white","brown1"),pch=c(NA,21,22,25),lwd=c(1,NA,NA,NA),pt.cex=1,cex=0.9)
text(-2.4,11.17,"(a)",cex=0.9)
text(-8.25,10.8,expression(paste("initial T"[d]," = 25\u00B0C, vapor ",delta^18,"O = -12\u2030")),cex=0.9)
text(-3,8.35,"start",cex=0.9)
text(-17.2,9.55,"end",cex=0.9)

#
plot(NA,NA,xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(-2,-18),
     ylim=c(10,30))
axis(1,seq(-18,-2,2),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(10,30,5),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

lines(rayleigh[[1]],get_D17O(rayleigh[[1]],rayleigh[[2]]))
points(d18O_1_p[2:length(d18O_1_p)],D17O_1_p[2:length(D17O_1_p)],pch=21,bg="gray")
points(d18O_2_p[2:length(d18O_2_p)],D17O_2_p[2:length(D17O_2_p)],pch=22,bg="white")
points(d18O_3_p[2:length(d18O_3_p)],D17O_3_p[2:length(D17O_3_p)],pch=25,bg="brown1")
text(-2.4,30.85,"(b)",cex=0.9)



# vapor input, T < 0C
d18O_input <- -25
dxs_input <- 10
D17O_input <- 10
d17O_input <- get_d17O(d18O_input,D17O_input)
d2H_input <- dxs_input+8*d18O_input
t_input <- 0

### case 1 ###
d18O_mix <- d18O_input
dxs_mix <- dxs_input
D17O_mix <- D17O_input
d17O_mix <- get_d17O(d18O_mix,D17O_mix)
d2H_mix <- dxs_mix+8*d18O_mix
temp_mix <- t_input
d18O_mix_p <- (d18O_mix+1000)*a18eff(273.15+temp_mix)-1000
d17O_mix_p <- (d17O_mix+1000)*a17eff(273.15+temp_mix)-1000
d2H_mix_p <- (d2H_mix+1000)*a2eff(273.15+temp_mix)-1000
dxs_mix_p <- d2H_mix_p-8*d18O_mix_p
D17O_mix_p <- get_D17O(d18O_mix_p,d17O_mix_p)
n <- 100
t_step <- 2
i <- 1
repeat{
  result <- Rayleigh_computer4(d18O_mix[i],dxs_mix[i],D17O_mix[i],temp_mix[i],temp_mix[i]-t_step,n,T,0.004)
  
  d18O_mix[i+1] <- weighted.mean(c(d18O_mix[i],result[[4]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d17O_mix[i+1] <- weighted.mean(c(d17O_mix[i],result[[5]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d2H_mix[i+1] <- weighted.mean(c(d2H_mix[i],result[[6]][n]),
                                c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  dxs_mix[i+1] <- d2H_mix[i+1]-8*d18O_mix[i+1]
  D17O_mix[i+1] <- get_D17O(d18O_mix[i+1],d17O_mix[i+1])
  temp_mix[i+1] <- weighted.mean(c(temp_mix[i],temp_mix[i]-t_step),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  
  d18O_mix_p[i+1] <- (d18O_mix[i+1]+1000)*a18eff(273.15+temp_mix[i+1])-1000
  d17O_mix_p[i+1] <- (d17O_mix[i+1]+1000)*a17eff(273.15+temp_mix[i+1])-1000
  d2H_mix_p[i+1] <- (d2H_mix[i+1]+1000)*a2eff(273.15+temp_mix[i+1])-1000
  dxs_mix_p[i+1] <- d2H_mix_p[i+1]-8*d18O_mix_p[i+1]
  D17O_mix_p[i+1] <- get_D17O(d18O_mix_p[i+1],d17O_mix_p[i+1])
  
  if ((temp_mix[i]-t_step)< -20) {break}
  i <- i+1
}
d18O_1 <- d18O_mix
D17O_1 <- D17O_mix
dxs_1 <- dxs_mix
d18O_1_p <- d18O_mix_p
D17O_1_p <- D17O_mix_p
dxs_1_p <- dxs_mix_p
#

### case 2 ###
d18O_mix <- d18O_input
dxs_mix <- dxs_input
D17O_mix <- D17O_input
d17O_mix <- get_d17O(d18O_mix,D17O_mix)
d2H_mix <- dxs_mix+8*d18O_mix
temp_mix <- t_input
d18O_mix_p <- (d18O_mix+1000)*a18eff(273.15+temp_mix)-1000
d17O_mix_p <- (d17O_mix+1000)*a17eff(273.15+temp_mix)-1000
d2H_mix_p <- (d2H_mix+1000)*a2eff(273.15+temp_mix)-1000
dxs_mix_p <- d2H_mix_p-8*d18O_mix_p
D17O_mix_p <- get_D17O(d18O_mix_p,d17O_mix_p)
n <- 100
t_step <- 6
i <- 1
repeat{
  result <- Rayleigh_computer4(d18O_mix[i],dxs_mix[i],D17O_mix[i],temp_mix[i],temp_mix[i]-t_step,n,T,0.004)
  
  d18O_mix[i+1] <- weighted.mean(c(d18O_mix[i],result[[4]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d17O_mix[i+1] <- weighted.mean(c(d17O_mix[i],result[[5]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d2H_mix[i+1] <- weighted.mean(c(d2H_mix[i],result[[6]][n]),
                                c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  dxs_mix[i+1] <- d2H_mix[i+1]-8*d18O_mix[i+1]
  D17O_mix[i+1] <- get_D17O(d18O_mix[i+1],d17O_mix[i+1])
  temp_mix[i+1] <- weighted.mean(c(temp_mix[i],temp_mix[i]-t_step),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  
  d18O_mix_p[i+1] <- (d18O_mix[i+1]+1000)*a18eff(273.15+temp_mix[i+1])-1000
  d17O_mix_p[i+1] <- (d17O_mix[i+1]+1000)*a17eff(273.15+temp_mix[i+1])-1000
  d2H_mix_p[i+1] <- (d2H_mix[i+1]+1000)*a2eff(273.15+temp_mix[i+1])-1000
  dxs_mix_p[i+1] <- d2H_mix_p[i+1]-8*d18O_mix_p[i+1]
  D17O_mix_p[i+1] <- get_D17O(d18O_mix_p[i+1],d17O_mix_p[i+1])
  
  if ((temp_mix[i]-t_step)< -20) {break}
  i <- i+1
}
d18O_2 <- d18O_mix
D17O_2 <- D17O_mix
dxs_2 <- dxs_mix
d18O_2_p <- d18O_mix_p
D17O_2_p <- D17O_mix_p
dxs_2_p <- dxs_mix_p
#

### case 3 ###
d18O_mix <- d18O_input
dxs_mix <- dxs_input
D17O_mix <- D17O_input
d17O_mix <- get_d17O(d18O_mix,D17O_mix)
d2H_mix <- dxs_mix+8*d18O_mix
temp_mix <- t_input
d18O_mix_p <- (d18O_mix+1000)*a18eff(273.15+temp_mix)-1000
d17O_mix_p <- (d17O_mix+1000)*a17eff(273.15+temp_mix)-1000
d2H_mix_p <- (d2H_mix+1000)*a2eff(273.15+temp_mix)-1000
dxs_mix_p <- d2H_mix_p-8*d18O_mix_p
D17O_mix_p <- get_D17O(d18O_mix_p,d17O_mix_p)
n <- 100
t_step <- 10
i <- 1
repeat{
  result <- Rayleigh_computer4(d18O_mix[i],dxs_mix[i],D17O_mix[i],temp_mix[i],temp_mix[i]-t_step,n,T,0.004)
  
  d18O_mix[i+1] <- weighted.mean(c(d18O_mix[i],result[[4]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d17O_mix[i+1] <- weighted.mean(c(d17O_mix[i],result[[5]][n]),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  d2H_mix[i+1] <- weighted.mean(c(d2H_mix[i],result[[6]][n]),
                                c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  dxs_mix[i+1] <- d2H_mix[i+1]-8*d18O_mix[i+1]
  D17O_mix[i+1] <- get_D17O(d18O_mix[i+1],d17O_mix[i+1])
  temp_mix[i+1] <- weighted.mean(c(temp_mix[i],temp_mix[i]-t_step),
                                 c(CC_relation(temp_mix[i]),CC_relation(temp_mix[i]-t_step)))
  
  d18O_mix_p[i+1] <- (d18O_mix[i+1]+1000)*a18eff(273.15+temp_mix[i+1])-1000
  d17O_mix_p[i+1] <- (d17O_mix[i+1]+1000)*a17eff(273.15+temp_mix[i+1])-1000
  d2H_mix_p[i+1] <- (d2H_mix[i+1]+1000)*a2eff(273.15+temp_mix[i+1])-1000
  dxs_mix_p[i+1] <- d2H_mix_p[i+1]-8*d18O_mix_p[i+1]
  D17O_mix_p[i+1] <- get_D17O(d18O_mix_p[i+1],d17O_mix_p[i+1])
  
  if ((temp_mix[i]-t_step)< -20) {break}
  i <- i+1
}
d18O_3 <- d18O_mix
D17O_3 <- D17O_mix
dxs_3 <- dxs_mix
d18O_3_p <- d18O_mix_p
D17O_3_p <- D17O_mix_p
dxs_3_p <- dxs_mix_p
#

### Rayleigh ###
rayleigh <- Rayleigh_computer4(d18O_input,dxs_input,D17O_input,t_input,-20,n,T,0.004)

#
plot(NA,NA,xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(-12,-33),
     ylim=c(4,12))
axis(1,seq(-30,-15,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(4,12,1),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste("d-excess (","\u2030",")")),2,padj=-1.7,cex=0.8)

lines(rayleigh[[1]],rayleigh[[3]]-rayleigh[[1]]*8)
points(d18O_1_p[2:length(d18O_1_p)],dxs_1_p[2:length(dxs_1_p)],pch=21,bg="gray")
points(d18O_2_p[2:length(d18O_2_p)],dxs_2_p[2:length(dxs_2_p)],pch=22,bg="white")
points(d18O_3_p[2:length(d18O_3_p)],dxs_3_p[2:length(dxs_3_p)],pch=25,bg="brown1")
text(-12.5,12.38,"(c)",cex=0.9)
text(-20,11.6,expression(paste("initial T"[d]," = 0\u00B0C, vapor ",delta^18,"O = -25\u2030")),cex=0.9)

#
plot(NA,NA,xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(-12,-33),
     ylim=c(10,60))
axis(1,seq(-30,-15,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(10,60,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)

lines(rayleigh[[1]],get_D17O(rayleigh[[1]],rayleigh[[2]]))
points(d18O_1_p[2:length(d18O_1_p)],D17O_1_p[2:length(D17O_1_p)],pch=21,bg="gray")
points(d18O_2_p[2:length(d18O_2_p)],D17O_2_p[2:length(D17O_2_p)],pch=22,bg="white")
points(d18O_3_p[2:length(d18O_3_p)],D17O_3_p[2:length(D17O_3_p)],pch=25,bg="brown1")
text(-12.5,62.2,"(d)",cex=0.9)