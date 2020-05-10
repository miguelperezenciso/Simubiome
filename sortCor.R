#------------------------------
sortCorr <- function(x, y, rho) {
#------------------------------
# returns list of permuted positions in y s.t. cor(x,y[sortCorr(x,y,rho)])=rho
   z = sort(y) + rnorm(length(y), 0, sd=sqrt(var(y) * (1./(rho**2)-1.)) )
   iy = order(y)
   iy = iy[order(z)]
   iy = iy[rank(x)]
   return(iy)
}

set.seed(13)
x = rnorm(1000)
y = rnorm(1000)
z = rgamma(1000,shape=2,rate=5)
print(cor(x,y))
print(cor(x,y[sortCorr(x,y,.3)]))
print(cor(x,z[sortCorr(x,z,.6)]))
