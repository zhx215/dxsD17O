rm(list=ls())

#
get_d17O <- function(d18O,D17O){
  1000*(exp(D17O/10^6+0.528*log(d18O/1000+1))-1)
}
get_D17O <- function(d18O,d17O){
  (log(d17O/1000+1)-0.528*log(d18O/1000+1))*10^6
}

# ET
a18_P_lv <- function(temp) {  # function for alpha of precipitation from vapor to liquid (or ice) phase (d18O) (>1)
  if (temp >= 273.15) {exp(1137/temp^2-0.4156/temp-0.00207)}
}
a17_P_lv <- function(temp) {  # function for alpha of precipitation from vapor to liquid (or ice) phase (d18O) (>1)
  if (temp >= 273.15) {exp(1137/temp^2-0.4156/temp-0.00207)^0.529}
}
a2_P_lv <- function(temp){  # function for alpha of precipitation from vapor to liquid (or ice) phase (dD) (>1)
  if (temp >= 273.15){exp(24844/temp^2-76.248/temp+0.05261)}
}
rayleigh_ET <- function(nd,t,h,delta_initial,alpha,k,smow){
  frac <-  1/(1/nd+1)
  r_initial <- (delta_initial/1000+1)*smow
  ae <- ((1-t)*((alpha*k*r_initial)/(1-h+(1-t)*k*h))+t*r_initial*(1/(1+(1-t)*k*(h/(1-h)))))/r_initial
  delta_ET <- ((delta_initial+1000)*(1-frac)^ae-1000*(1-frac)-delta_initial)/(1-frac-1)*(1-t)+delta_initial*t
  return(c((delta_initial-frac*delta_ET)/(1-frac),delta_ET))
}

#
diff18 <- 0.9723
diff17 <- (0.9723^0.518)
diff2 <- 0.9755
smow18 <- 0.0020052
smow17 <- 379.9*10^-6
smow2 <- 155.76*10^-6

#
et_p_ratio <- c(0.3,0.6,0.9)
nd <- et_p_ratio/(1-et_p_ratio)
y <- 1/(nd+1)
eet <- c(0.8,0.5,0.2)
ck_18 <- (1/(0.9723^0.5)-1)*1000
ck_17 <- (1/((0.9723^0.518)^0.5)-1)*1000
ck_2 <- (1/(0.9755^0.5)-1)*1000
epe_18 <- (a18_P_lv(22+273.15)-1)*1000
epe_17 <- (a17_P_lv(22+273.15)-1)*1000
epe_2 <- (a2_P_lv(22+273.15)-1)*1000
delta_p_18 <- -5
delta_p_17 <- get_d17O(delta_p_18,20)
delta_p_2 <- delta_p_18*8+10
h <- 0.7
color=c("#0072B2","#009E73","#E69F00")
code <- c(24,22,25)

#
d18O <- dxs <- D17O <- array(NA,dim=c(3,3,3))
# model A
for (i in 1:3){
  for (j in 1:3){
    x <- eet[j]*(1-y[i])
    z <- y[i]/(x+y[i])
    d18O_e <- delta_p_18-z*(1-h)/(1-z*h)*(ck_18+epe_18)
    d18O_t <- delta_p_18
    d18O_et <- (d18O_e*x+d18O_t*(1-x-y[i]))/(1-y[i])
    d17O_e <- delta_p_17-z*(1-h)/(1-z*h)*(ck_17+epe_17)
    d17O_t <- delta_p_17
    d17O_et <- (d17O_e*x+d17O_t*(1-x-y[i]))/(1-y[i])
    d2H_e <- delta_p_2-z*(1-h)/(1-z*h)*(ck_2+epe_2)
    d2H_t <- delta_p_2
    d2H_et <- (d2H_e*x+d2H_t*(1-x-y[i]))/(1-y[i])
    d18O[1,i,j] <- d18O_et
    dxs[1,i,j] <- d2H_et-8*d18O_et
    D17O[1,i,j] <- get_D17O(d18O_et,d17O_et)
  }  
}
# model B
for (i in 1:3){
  for (j in 1:3){
    x <- eet[j]*(1-y[i])
    z <- 1-x
    d18O_e <- delta_p_18-z*(1-h)/(1-z*h)*(ck_18+epe_18)
    d18O_t <- (1-z)*(1-h)*(ck_18+epe_18)/(1-z*h)+delta_p_18
    d18O_et <- (d18O_e*x+d18O_t*(1-x-y[i]))/(1-y[i])
    d17O_e <- delta_p_17-z*(1-h)/(1-z*h)*(ck_17+epe_17)
    d17O_t <- (1-z)*(1-h)*(ck_17+epe_17)/(1-z*h)+delta_p_17
    d17O_et <- (d17O_e*x+d17O_t*(1-x-y[i]))/(1-y[i])
    d2H_e <- delta_p_2-z*(1-h)/(1-z*h)*(ck_2+epe_2)
    d2H_t <- (1-z)*(1-h)*(ck_2+epe_2)/(1-z*h)+delta_p_2
    d2H_et <- (d2H_e*x+d2H_t*(1-x-y[i]))/(1-y[i])
    d18O[2,i,j] <- d18O_et
    dxs[2,i,j] <- d2H_et-8*d18O_et
    D17O[2,i,j] <- get_D17O(d18O_et,d17O_et)
  }  
}
# model C
for (i in 1:3){
  for (j in 1:3){
    x <- eet[j]*(1-y[i])
    p <- (2*x+y[i])/2
    z <- (p-x)/p
    d18O_e <- delta_p_18-z*(1-h)/(1-z*h)*(ck_18+epe_18)
    d18O_t <- delta_p_18
    d18O_et <- (d18O_e*x+d18O_t*(1-x-y[i]))/(1-y[i])
    d17O_e <- delta_p_17-z*(1-h)/(1-z*h)*(ck_17+epe_17)
    d17O_t <- delta_p_17
    d17O_et <- (d17O_e*x+d17O_t*(1-x-y[i]))/(1-y[i])
    d2H_e <- delta_p_2-z*(1-h)/(1-z*h)*(ck_2+epe_2)
    d2H_t <- delta_p_2
    d2H_et <- (d2H_e*x+d2H_t*(1-x-y[i]))/(1-y[i])
    d18O[3,i,j] <- d18O_et
    dxs[3,i,j] <- d2H_et-8*d18O_et
    D17O[3,i,j] <- get_D17O(d18O_et,d17O_et)
  }  
}

lines(c(10,10),c(-17.5,-4.5))
lines(c(20,20),c(-17.5,-4.5))
text(5,-17.05,"T/ET = 20%",cex=0.9)
text(15,-17.05,"T/ET = 50%",cex=0.9)
text(25,-8,"T/ET = 80%",cex=0.9)

legend("bottomright",c("n = 0.5","n = 0.75","n = 1","ET/P = 0.3","ET/P = 0.6","ET/P = 0.9","no WS"),
       pch=c(code,rep(22,4)),cex=0.9,
       pt.bg=c(rep("gray",3),color,"white"))

text(-1,-4,expression(paste("(a) ",delta^18,"O"," (\u2030)")),cex=0.9,pos=4)


#
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)

# 1
plot(NA,NA,xlim=c(0,30),ylim=c(-14.5,-4.5),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(-14,-6,2),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,d18O[k,i,j],bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}
lines(c(10,10),c(-14.5,-4.5))
lines(c(20,20),c(-14.5,-4.5))
text(5,-14.1,"T/ET = 20%",cex=0.9)
text(15,-14.1,"T/ET = 50%",cex=0.9)
text(25,-9.3,"T/ET = 80%",cex=0.9)

legend("bottomright",c("model A","model B","model C","ET/P = 0.3","ET/P = 0.6","ET/P = 0.9"),
       pch=c(code,rep(22,3)),cex=0.9,
       pt.bg=c(rep("gray",3),color))

text(-1,-4.05,expression(paste("(a) ",delta^18,"O"," (\u2030)")),cex=0.9,pos=4)

# 2
plot(NA,NA,xlim=c(0,30),ylim=c(7,50),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(10,50,10),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,dxs[k,i,j],bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(7,50))
lines(c(20,20),c(7,50))
text(5,8.6,"T/ET = 20%",cex=0.9)
text(15,8.6,"T/ET = 50%",cex=0.9)
text(25,8.6,"T/ET = 80%",cex=0.9)

text(-1,51.7,expression(paste("(b) d-excess (","\u2030",")")),cex=0.9,pos=4)


# 3
plot(NA,NA,xlim=c(0,30),ylim=c(15,120),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(20,120,20),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,D17O[k,i,j],bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(15,120))
lines(c(20,20),c(15,120))
text(5,19,"T/ET = 20%",cex=0.9)
text(15,19,"T/ET = 50%",cex=0.9)
text(25,19,"T/ET = 80%",cex=0.9)

text(-1,125,expression(paste("(c) ",Delta*minute^17,"O"," (per meg)")),cex=0.9,pos=4)


# 4
plot(NA,NA,xlim=c(0,30),ylim=c(1.5,3),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(1.5,3,0.5),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,(D17O[k,i,j]-20)/(dxs[k,i,j]-10),bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(1.5,3))
lines(c(20,20),c(1.5,3))
text(5,1.56,"T/ET = 20%",cex=0.9)
text(15,1.56,"T/ET = 50%",cex=0.9)
text(25,1.56,"T/ET = 80%",cex=0.9)

text(-1,3.07,expression(paste("(d) ",Delta,"(",Delta*minute^17,"O)/",
                              Delta,"(d-excess)")),cex=0.9,pos=4)