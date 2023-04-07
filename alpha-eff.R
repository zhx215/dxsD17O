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

temp <- seq(-50,30,0.1)
a18 <- a17 <- a2 <- vector()
a18_2 <- a17_2 <- a2_2 <- vector()
a18_6 <- a17_6 <- a2_6 <- vector()
a18_d <- a17_d <- a2_d <- vector()
for (i in 1:length(temp)){
  a18[i] <- a18eff(temp[i]+273.15,"k",273.15,250.15)
  a17[i] <- a17eff(temp[i]+273.15,"k",273.15,250.15)
  a2[i] <- a2eff(temp[i]+273.15,"k",273.15,250.15)
  a18_2[i] <- a18eff(temp[i]+273.15,"k",273.15,250.15,0.002)
  a17_2[i] <- a17eff(temp[i]+273.15,"k",273.15,250.15,0.002)
  a2_2[i] <- a2eff(temp[i]+273.15,"k",273.15,250.15,0.002)
  a18_6[i] <- a18eff(temp[i]+273.15,"k",273.15,250.15,0.006)
  a17_6[i] <- a17eff(temp[i]+273.15,"k",273.15,250.15,0.006)
  a2_6[i] <- a2eff(temp[i]+273.15,"k",273.15,250.15,0.006)
  a18_d[i] <- a18eff(temp[i]+273.15,"k",263.15,233.15)
  a17_d[i] <- a17eff(temp[i]+273.15,"k",263.15,233.15)
  a2_d[i] <- a2eff(temp[i]+273.15,"k",263.15,233.15)
}

#
par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
plot(temp,a18,type="l",xlim=c(-50,30),ylim=c(1.0088,1.021),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-50,30,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(1.01,1.02,0.002),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste(""^18,alpha[eff])),2,padj=-1.2,cex=0.8)
lines(temp,a18_2,lty=2)
lines(temp,a18_6,lty=3)
lines(temp,a18_d,lty=4)
text(-48,1.0216,"(a)",cex=0.9)

plot(temp,a17,type="l",xlim=c(-50,30),ylim=c(1.0047,1.011),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-50,30,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(1.005,1.011,0.001),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste(""^17,alpha[eff])),2,padj=-1.2,cex=0.8)
lines(temp,a17_2,lty=2)
lines(temp,a17_6,lty=3)
lines(temp,a17_d,lty=4)
text(-48,1.01129,"(b)",cex=0.9)

plot(temp,a2,type="l",xlim=c(-50,30),ylim=c(1.074,1.23),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-50,30,10),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(1.08,1.22,0.02),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste("T"[d]," (\u00B0C)")),1,padj=1.9,cex=0.8)
mtext(expression(paste(""^2,alpha[eff])),2,padj=-1.2,cex=0.8)
lines(temp,a2_2,lty=2)
lines(temp,a2_6,lty=3)
lines(temp,a2_d,lty=4)
text(-48,1.238,"(c)",cex=0.9)

plot.new()
legend("topleft",c("standard",
                   expression(paste(lambda," = 0.002")),
                   expression(paste(lambda," = 0.006")),
                   "delayed ice-vapor fractionation"),
       lty=c(1,2,3,4))