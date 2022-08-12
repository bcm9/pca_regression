# pca_regression
principal component analysis w/regression in MATLAB


pca_regression Function runs principal component analysis and conducts PC
regression with nPC components, including plots


   Arguments:
   
   
   Xdata = observations x variables input data matrix
   
  
   Ydata = response variable array
   
   
   nPC = number of PCs to regress



   Outputs:
   
   
   coeff_loadings = Principal component coefficients are recipe for counting any given PC
       each column of coeff contains coefficients for one principal component.
       columns are in order of descending component variance, latent.


   score = how each individual observation is composed of the PCs. matrix of PCs x observations.


   latent = principal component variances, the eigenvalues of the covariance matrix, returned as a column vector
       an eigenvalue is the total amount of variance in the variables in the dataset explained by the common factor.
       PCs are typically retained for further analysis when their eigenvalues (of the transformation matrix) are larger than 1.


   explained = the contribution of each PC to variability in data


   rsq_PCR = r2 for PCR
