temp <- seq(from=1.029e+01, to=3.872e+01, length=50) #get from raster
library(rjags)


tmax <- max(temp) #get from raster
tmin<-min(temp) #get from raster

m<-2
a<-0.5
c_est <- 5.605e-04
#Ellie's Briere equation
y_ellie <- c_est*(temp-tmin)*(tmax-temp)^(1/m)

#Briere equation from paper
y_paper <- c_est*temp*((temp-tmin)*sqrt((tmax-temp))*(tmax>temp))
plot(y_paper/max(y_paper)~temp)


plot(y_ellie/max(y_ellie)~temp)

mu.temp[i] <- c*T[i]*(T[i]-T0)*sqrt((Tm-T[i])*(Tm>T[i]))





#just substitute my raster for the temp in this
#but we want to use data to estimate this shape
#and find whatever the values of m and a should be
#using data from the ecology paper
#temp would be t col, y would be the trait column (e2a)
#use the briere_1999 function that we've found to give this a wee go
#divide by max value of whatever y comes out of and should produce 0-1
