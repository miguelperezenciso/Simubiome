# Simubiome
### Simulation of complex phenotypes mediated by genome and/or microbiome.
#### M Pérez-Enciso with help from LM Zingaretti
#### miguel.perez@uab.es
#### CITATION: This is accompanying code to
Opportunities and limits of using microbiome data for complex trait prediction. M. Pérez-Enciso, L.M. Zingaretti, Y. Ramayo-Caldas, G. de los Campos (submitted)

***

### Sortcorr
This function takes two vectors x and y and allows reordering y such that cor(x,reordered(y)) = rho, where rho is the desired correlation between x and y. This function is critical to reorder abundances such that a correlation between genotypic values and abundances is established. In so doing, a genetic basis of given abundances is mimicked without modifying the observed data. OTU reordering can be done within cluster of correlated abundances in order to preserve - to an extent - the covariance between different OTU abundances.

    #------------------------------
    sortCorr <- function(x, y, rho) {
    #------------------------------
    # returns list of permuted positions of y s.t. cor(x,y[sortCorr(x,y,rho)])=rho
      z = sort(y) 
      z = z + rnorm(length(y), 0, sd=sqrt(var(y) * (1./(rho**2)-1.)))
      iy = order(y)
      iy = iy[order(z)]
      iy = iy[rank(x)]
      return(iy)
    }
  
    # repeatable results
    set.seed(3)
    x = rnorm(1000)
    y = rgamma(1000, 1, 1)
    # -0.08020344
    cor(x,y)
    # 0.5104921
    cor(x,y[sortCorr(x,y,rho=0.50)]) 
    #  0.09558201
    cor(x,y[sortCorr(x,y,rho=0.10)])  

