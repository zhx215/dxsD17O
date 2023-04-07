# 17O converter
get_D17O <- function(d18O,d17O){
  (log(d17O/1000+1)-0.528*log(d18O/1000+1))*10^6
}

get_d17O <- function(d18O,D17O){
  1000*(exp(D17O/10^6+0.528*log(d18O/1000+1))-1)
}

par(mfrow=c(2,2),mar=c(3,3,2,2),xpd=T)
#
plot(NA,NA,xlim=c(-35,-5),ylim=c(5,55),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(-35,-5,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(10,50,10),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(delta^18,"O (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(paste(Delta*minute^17,"O"," (per meg)")),2,padj=-1.3,cex=0.8)
text(-36,57.2,"(a)",cex=0.9,pos=4)

# case 1
d18O_1 <- -10
D17O_1 <- 20
d17O_1 <- get_d17O(d18O_1,D17O_1)
d18O_2 <- -30
D17O_2 <- 20
d17O_2 <- get_d17O(d18O_2,D17O_2)
v1 <- seq(0,1,0.01)
d18O_mix <- v1*d18O_1+(1-v1)*d18O_2
d17O_mix <- v1*d17O_1+(1-v1)*d17O_2
D17O_mix <- get_D17O(d18O_mix,d17O_mix)

lines(d18O_mix,D17O_mix)
lines(d18O_mix,v1*D17O_1+(1-v1)*D17O_2,col="gray",lty=2)
points(d18O_2,D17O_2,pch=21,bg="white",cex=1.2)
text(d18O_2-1.5,D17O_2,"2a",cex=0.9)

# case 2
d18O_1 <- -10
D17O_1 <- 20
d17O_1 <- get_d17O(d18O_1,D17O_1)
d18O_2 <- -30
D17O_2 <- 50
d17O_2 <- get_d17O(d18O_2,D17O_2)
v1 <- seq(0,1,0.01)
d18O_mix <- v1*d18O_1+(1-v1)*d18O_2
d17O_mix <- v1*d17O_1+(1-v1)*d17O_2
D17O_mix <- get_D17O(d18O_mix,d17O_mix)

lines(d18O_mix,D17O_mix)
lines(d18O_mix,v1*D17O_1+(1-v1)*D17O_2,col="gray",lty=2)
points(d18O_2,D17O_2,pch=21,bg="white",cex=1.2)
text(d18O_2-1.5,D17O_2,"2b",cex=0.9)

# case 3
d18O_1 <- -10
D17O_1 <- 20
d17O_1 <- get_d17O(d18O_1,D17O_1)
d18O_2 <- -20
D17O_2 <- 50
d17O_2 <- get_d17O(d18O_2,D17O_2)
v1 <- seq(0,1,0.01)
d18O_mix <- v1*d18O_1+(1-v1)*d18O_2
d17O_mix <- v1*d17O_1+(1-v1)*d17O_2
D17O_mix <- get_D17O(d18O_mix,d17O_mix)

lines(d18O_mix,D17O_mix)
lines(d18O_mix,v1*D17O_1+(1-v1)*D17O_2,col="gray",lty=2)
points(c(d18O_1,d18O_2),c(D17O_1,D17O_2),pch=21,bg="white",cex=1.2)
text(-9,20,"1");text(d18O_2-1.5,D17O_2,"2c",cex=0.9)

text(-12.5,12,"mixing line",srt=45,cex=0.9)

# offset plot
d18O_1 <- -10
D17O_1 <- 20
d18O_offset <- seq(0,30,0.5)
d18O_2 <- d18O_1-d18O_offset
D17O_2 <- 20
v_1 <- seq(0,1,0.01)

D17O_offset <- array(NA,dim=c(length(d18O_2),length(v_1)))
for (i in 1:length(d18O_2)){
  for (j in 1:length(v_1)){
    d18O_mix <- d18O_1*v_1[j]+d18O_2[i]*(1-v_1[j])
    d17O_mix <- get_d17O(d18O_1,D17O_1)*v_1[j]+get_d17O(d18O_2[i],D17O_2)*(1-v_1[j])
    D17O_offset[i,j] <- get_D17O(d18O_mix,d17O_mix)-(D17O_1*v_1[j]+D17O_2*(1-v_1[j]))
  }
}

plot(NA,NA,xlim=c(0,30),ylim=c(0,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
axis(1,seq(0,30,5),cex.axis=0.9,padj=-1.3,tck=-0.02)
axis(2,seq(0,1,0.2),cex.axis=0.9,padj=1.1,tck=-0.02)
mtext(expression(paste(Delta,"(",delta^18,"O) (","\u2030",")")),1,padj=1.5,cex=0.8)
mtext(expression(M[1]/(M[1]+M[2])),2,padj=-1.5,cex=0.8)
contour(d18O_offset,v_1,D17O_offset,levels=c(-1,-3,-5,-10,-15,-20,-25),labcex=0.8,vfont=NULL,add=T)
text(-1.1,1.042,expression(paste("(b) vapor mixing ",Delta,"(",Delta*minute^17,"O)"," (per meg)")),cex=0.9,pos=4)