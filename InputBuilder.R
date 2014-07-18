setwd("C:/Users/KingJ/Documents/Projects/MEMLS_Active/memls_TVC/input")
ssa <-read.csv(file="2LayerSSA.csv",head=TRUE,sep=",")
trench <-read.csv(file="Trench7Data.csv",head=TRUE,sep=",")

slab = ssa[ssa$Layer == 1,]
hoar = ssa[ssa$Layer == 2,]

numLayer = 2
members = 1000
nmembers = 0

depth.total = rweibull(members, 6.5539951, 31.1303054)
depth.dhFract = rweibull(members, 2.2526008, 29.2109527)

#Assuming a dry terrestrial environment build constants
snowSal <- rep(0, times = numLayer)
snowWet <- rep(0, times = numLayer)

snowTemp <- rep(250.00, times =numLayer) #not so important for TVC
snowRho <- c(362.00,250.00)

snowRhoV <- c((snowRho[1]/917.00), (snowRho[2]/917.00))


#Assemble D_0 distrobutions D_0 = 6/rho_Ice*SSA
hoar.d0 = 6/(0.917*rlnorm(members, 2.18378849, 0.17641227))
slab.d0 = 6/(0.917*rlnorm(members, 3.11093875, 0.44167987))


#Assmble P_ex input distrobutions
slab.pex = 0.5*slab.d0*(1-snowRhoV[1])
hoar.pex = 0.5*hoar.d0*(1-snowRhoV[2])



for (i in 1:members) {
  snowPEx <-c(slab.pex[i],  hoar.pex[i])
  snowThick <-c(depth.total[i]-(depth.total[i]*(depth.dhFract[i]/100)),depth.total[i]*(depth.dhFract[i]/100))

memls_input_top = paste(1,snowTemp[1], snowWet[1], snowRho[1], snowThick[1], snowSal[1], snowPEx[1], sep=" ")
memls_input_bottom = paste(2,snowTemp[2], snowWet[2], snowRho[2], snowThick[2], snowSal[2], snowPEx[2], sep=" ")

print(memls_input_top)
print(memls_input_bottom)



if(dh_fraction.lnormal[i] > 0 && dh_fraction.lnormal[i] < 100){
  nmembers = nmembers  + 1
  filename = paste("TVC",nmembers,".txt", sep ="")
  fileConn<-file(filename)
  writeLines(c(memls_input_top, memls_input_bottom), sep = "\n",  fileConn)
  close(fileConn)
 }

}


par(mfrow=c(2,2))

plot(density(trench$Depth))
lines(density(depth.total), col = "red")

plot(density(trench$DH_Fract))
lines(density(depth.dhFract/100), col = "red")

plot(density(hoar$D_0))
lines(density(hoar.d0), col = "red")

plot(density(slab$D_0))
lines(density(slab.d0), col = "red")


plot(depth.total, type = "l")

library(MASS)
x <- 1:length(depth.total)
y <- depth.total
lo <- loess(y~x)
plot(x,y, type = "l")
lines(predict(lo), col='red', lwd=2)