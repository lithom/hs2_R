#
# Here we test the functions in hrstuff_01.r
#



#install.packages("rgl")
#install.packages("pracma")

library(rgl)
library(pracma)

#C:\Users\lithomas\
source("hrstuff_01.r")


#
# TESTS
#
G1 = rbind(diag(1,3),-diag(1,3))
h1 = rbind(1,2,3,4,5,6);

X1 = cbind(c(0.5,0.5,0.5),c(0.7,0.7,0.7))
V1 = cbind(c(1,1,1),c(1,1,1.5))

xx1 <- intersectLineWithPlanes(G1,h1,X1,V1)

xx1[2]

hrx1 = hrstep_unif(G1,h1,X1)

hrx2 = hr_unif(G1,h1,X1,c(5))
hrx3 = hr_unif(G1,h1,X1,1:500)



hrx3_j = Reduce(cbind,hrx3)

plot3d(hrx3_j[1,],hrx3_j[2,],hrx3_j[3,])




#
# TEstS NOn-Unif:
#

G2 = rbind(diag(1,2),-diag(1,2))
h2 = rbind(1,2,3,4);


fa = function(x) { return( dnorm(x[1])*dnorm(x[2]) ) }

x2_a = hrstep_01( G2,h2,fa, cbind(c(0.1,0.1),c(0.2,0.2)) ,40)

hrgx_a = hr(G2,h2,fa,cbind(c(0.1,0.1),c(0.2,0.2)),200,40)


