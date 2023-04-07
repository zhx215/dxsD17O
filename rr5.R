rm(list=ls())
setwd("C:\\Users\\zhyxi\\Desktop\\17O project")
source("raindrop re-evaporation model.R")
#
rh <- seq(0.4,1,0.0025)
temp <- seq(0,44,0.25)

d <- seq(0.4,2.6,0.2)
rr_result_matrix <- array(NA,dim=c(length(temp),length(rh),length(d),5))
count <- 0
for (i in 128:length(temp)){
  for (j in 1:length(rh)){
    for (k in 1:length(d)){
      rr_result_matrix[i,j,k,] <- reevap(101325,temp[i]+273.15,rh[j],d[k],0,0,0,3,0.004)
      count <- count + 1
      print(count/(length(temp)*length(rh)*length(d))*100)
    }
  }
}
print(Sys.time())

#
load("C:/Users/zhyxi/Desktop/17O project/rr5_v1.RData")

#
rh <- seq(0.4,1,0.0025)
temp <- seq(-64,40,0.25)
td_matrix <- array(NA,dim=c(length(temp),length(rh)))
for (i in 1:length(temp)){
  for (j in 1:length(rh)){
    td_matrix[i,j] <- lcl(101325,temp[i]+273.15,rhl=rh[j])[2]-273.15
  }
}

infer_rh <- function(t,td){
  id_t <- which.min(abs(seq(-64,40,0.25)-t))
  id_rh <- which.min(abs((td_matrix[id_t,]-td)))
  return(seq(0.4,1,0.0025)[id_rh])
}

#
rm(list=setdiff(ls(), c("rr_result_matrix","infer_rh","td_matrix")))
save.image("C:/Users/zhyxi/Desktop/17O project/RR and rh_5.RData")
