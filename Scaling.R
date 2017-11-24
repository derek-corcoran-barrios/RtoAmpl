alphamax <- function(m){(3.5 * 7.5  * (10^-4) * (m^-0.25) )* 3650}

DistMax <- function(m){3.31*(m^0.65)*10 * 1000}
DistMedian <- function(m){1.45*(m^0.54)}
Prd <- function(m, d){exp(-(d/DistMedian(m)))}

#peso 0.199 kgs

DistMedian(0.199)
0.606

DistMax(0.199)
1.156

Prd(m = 0.199, d = 0.606)


Prd(m = 0.199, d = 0)

Hendricks
#k = -0.05
#Gamma =91

Density <- function(m) {m*4.23*(m^-0.75)}


Density(50)
