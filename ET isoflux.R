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
eet <- c(0.8,0.5,0.2)
ae <- c(0.5,0.75,1)

#
delta_p_18 <- -5
delta_p_17 <- get_d17O(delta_p_18,20)
delta_p_2 <- delta_p_18*8+10
h <- 0.7
color=c("#0072B2","#009E73","#E69F00")
code <- c(24,22,25)

#
dxs_ws <- D17O_ws <- d18O_ws <- d17O_ws <- d2H_ws <- array(NA,dim=c(3,3,3,2))
sp_ratio <- c(0.2,2)
for (i in 1:3){
  for (j in 1:3){
    for (k in 1:3){
      for (l in 1:2){
        num <- 1
        d18Os <- delta_p_18
        d17Os <- delta_p_17
        d2Hs <- delta_p_2
        repeat{
          d18O_mix <- d18Os[num]*(1-1/(sp_ratio[l]+1))+delta_p_18*(1/(sp_ratio[l]+1))
          result18 <- rayleigh_ET(et_p_ratio[i]/sp_ratio[l],1-eet[j],h,d18O_mix,1/a18_P_lv(22+273.15),diff18^ae[k],smow18)
          
          d17O_mix <- d17Os[num]*(1-1/(sp_ratio[l]+1))+delta_p_17*(1/(sp_ratio[l]+1))
          result17 <- rayleigh_ET(et_p_ratio[i]/sp_ratio[l],1-eet[j],h,d17O_mix,1/a17_P_lv(22+273.15),diff17^ae[k],smow17)
          
          d2H_mix <- d2Hs[num]*(1-1/(sp_ratio[l]+1))+delta_p_2*(1/(sp_ratio[l]+1))
          result2 <- rayleigh_ET(et_p_ratio[i]/sp_ratio[l],1-eet[j],h,d2H_mix,1/a2_P_lv(22+273.15),diff2^ae[k],smow2)
          
          num <- num+1
          d18Os[num] <- result18[1]
          d17Os[num] <- result17[1]
          d2Hs[num] <- result2[1]
          dxs <- d2Hs-8*d18Os
          D17O <- get_D17O(d18Os,d17Os)
          if ((abs(dxs[num]-dxs[num-1])<0.01)&(abs(D17O[num]-D17O[num-1])<0.01)){
            break
          }
        }
        dxs_ws[i,j,k,l] <- result2[2]-8*result18[2]
        D17O_ws[i,j,k,l] <- get_D17O(result18[2],result17[2])
        d18O_ws[i,j,k,l] <- result18[2]
        d17O_ws[i,j,k,l] <- result17[2]
        d2H_ws[i,j,k,l] <- result2[2]
      }
    }
  }
}

#
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)

# 1
plot(NA,NA,xlim=c(0,30),ylim=c(-16,-4),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(-16,-4,2),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,d18O_ws[i,j,k,2],col=color[i],pch=code[k])
      points(number,d18O_ws[i,j,k,1],bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(-16,-4))
lines(c(20,20),c(-16,-4))
text(5,-15.5,"T/ET = 20%",cex=0.9)
text(15,-15.5,"T/ET = 50%",cex=0.9)
text(25,-8,"T/ET = 80%",cex=0.9)

legend("bottomright",c("ET/P = 0.3","ET/P = 0.6","ET/P = 0.9","S/P = 0.2","S/P = 2","m = 0.5","m = 0.75","m = 1"),
       pch=c(rep(22,5),code),cex=0.9,
       pt.bg=c(color,"black","white",rep("gray",3)))

text(-1,-3.48,expression(paste("(a) ",delta^18,"O"," (\u2030)")),cex=0.9,pos=4)

# 2
plot(NA,NA,xlim=c(0,30),ylim=c(7,45),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(10,45,5),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,dxs_ws[i,j,k,2],col=color[i],pch=code[k])
      points(number,dxs_ws[i,j,k,1],bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(7,45))
lines(c(20,20),c(7,45))
text(5,8.6,"T/ET = 20%",cex=0.9)
text(15,8.6,"T/ET = 50%",cex=0.9)
text(25,8.6,"T/ET = 80%",cex=0.9)

text(-1,46.4,expression(paste("(b) d-excess (","\u2030",")")),cex=0.9,pos=4)

# 3
plot(NA,NA,xlim=c(0,30),ylim=c(15,70),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(20,70,10),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,D17O_ws[i,j,k,2],col=color[i],pch=code[k])
      points(number,D17O_ws[i,j,k,1],bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(15,70))
lines(c(20,20),c(15,70))
text(5,17.2,"T/ET = 20%",cex=0.9)
text(15,17.2,"T/ET = 50%",cex=0.9)
text(25,17.2,"T/ET = 80%",cex=0.9)

text(-1,72.5,expression(paste("(c) ",Delta*minute^17,"O"," (per meg)")),cex=0.9,pos=4)

# 4
plot(NA,NA,xlim=c(0,30),ylim=c(-0.3,2.7),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(2,seq(0,2.5,0.5),cex.axis=0.9,padj=1.1,tck=-0.02)

number <- 1
for (j in 1:3){
  for (i in 1:3){
    for (k in 1:3){
      points(number,(D17O_ws[i,j,k,2]-20)/(dxs_ws[i,j,k,2]-10),col=color[i],pch=code[k])
      points(number,(D17O_ws[i,j,k,1]-20)/(dxs_ws[i,j,k,1]-10),bg=color[i],pch=code[k])
      number <- number +1
    }
  }
  number <- number+1
}

lines(c(10,10),c(-0.3,2.7))
lines(c(20,20),c(-0.3,2.7))
text(5,-0.175,"T/ET = 20%",cex=0.9)
text(15,-0.175,"T/ET = 50%",cex=0.9)
text(25,-0.175,"T/ET = 80%",cex=0.9)

text(-1,2.83,expression(paste("(d) ",Delta,"(",Delta*minute^17,"O)/",
                                 Delta,"(d-excess)")),cex=0.9,pos=4)