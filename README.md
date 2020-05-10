# Simubiome
### Simulation of complex phenotypes mediated by genome and/or microbiome.
#### M Pérez-Enciso (miguel.perez@uab.es) with help from LM Zingaretti

#### CITATION: This is accompanying code to
Opportunities and limits of using microbiome data for complex trait prediction. M. Pérez-Enciso, L.M. Zingaretti, Y. Ramayo-Caldas, G. de los Campos (submitted)

***

### Input
As input, the script requires a nxp **X** matrix containing genotypes and nxq **B** matrix containing OTU abundances. Parameters are
* h2: heritability
* b2: microbiability (b2+h2<1)
* Nqtl_y: number of causative SNPs
* Notu_y: number of causative OTUs
* Notu_y_g: number of causative OTUs with a genetic basis (Notu_y_g <= Notu_y). Each of those is assumed to be determined by Nqtl_otu=10 SNPs (independently sampled of Nqtl_y SNPs). Heritabilities are sampled from gamma distributions (see code).

Effects are sampled from gamma distributions (see publication and code for details).

### Models and usage
Simubiome allows five generic causality models (Fig. 1). 

* Joint model: Both genome **G** and microbiome **B** act independently on **y**. It is perhaps the most widely assumed model, implicitly or explicitly, in the literature.
* Recursive model:
* Microbiome model: Only **B** has an effect on **y**, and **G** is only noise.
* Indirect model:
* Genome model: Only **G** has an effect on **y**, and **B** is only noise.

* . The amount of variance explained by G and B are determined by the heritability (h2) and microbiability (b2), respectively. This is the simplest model once we know the trait is under the influence of both heredity and microbiome. It is perhaps the most widely assumed model, implicitly or explicitly, in the literature, e.g., [11,18,20]. Simplifications of the Joint model lead to the Genome and Microbiome models when either the microbiome or the genome do not affect the phenotype, respectively. In any case, we take that SNP and microbiome data are available even if they may simply add noise to the analysis. The Recursive model is similar to the Joint model, except that allows for some causative OTU abundances to be under partial genetic control. Therefore, the influence of the genome is both direct and indirect (microbiome mediated). Note that we do not assume that a given SNP has simultaneously direct and indirect effects, we rather hypothesize that the genes controlling OTU abundance are different from those with a direct effect on the phenotype. We also do not assume that all OTU abundances are under genetic control. It may be that non causative OTUs are under genetic control but these are irrelevant for our purposes and are not modeled here. In the Indirect model, finally, the genome affects the phenotype but only indirectly. In this case, the phenotype is controlled by a number of causative OTUs and part of these OTUs’ abundances are under genetic control. This would be equivalent to a scenario where a phenotype is controlled by gene expression levels and expression is in turn controlled genetically.

    s = SimuBiome(X, B, Bclust=Bclust, h2=h2, b2=b2, Nqtl_y=Nqtl_y, Notu_y=Notu_y, Notu_y_g=Notu_y_g)

### Output

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

