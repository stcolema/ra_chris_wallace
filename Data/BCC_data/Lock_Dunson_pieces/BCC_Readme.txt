The function "BCC(...)" performs Bayesian consensus clustering.

To use BCC, first load the script 'BCC.r' (e.g., enter "source('BCC.r')" if the BCC.r file is in your working directory).

The function inputs are 
BCC(X,K=2,a=1,b=1,IndivAlpha = FALSE, mu0 = list(),a0=list(),b0=list(),Concentration=1,NumDraws = 1000)
where
-X is a list of quantitative data sources, X = list(X1,...,XM),  Where each Xi is a D_i X N data matrix (D_i variables on N columns).  
The columns should match for each data source.
-K is the (maximum) number of clusters
-a and b are hyperparameters for the Beta(a,b) prior distribution on alpha
-IndivAlpha indicates whether the alpha should be separate for each data source ("TRUE" or "FALSE")
-mu0 is a list of the prior mean vectors for each data source (i.e. mu0_i is the prior mean for the ith data source)
-a0 and b0 are lists of the Gamma parameters (to determine variance) for each data source
-Concentration is the Dirichlet concentration parameter for the overall cluster sized
-NumDraws is the number of MCMC draws (NumDraws/2 is used as the "burn-in")


The outputs of BCC are 
-Alpha: the average adherence (alpha).  If IndivAlpha = TRUE, Alpha[m] gives average adherence for data source m.
-AlphaBounds:95% credible interval for Alpha
-CBest: The "hard" overall clustering, as a binary matrix.  Cbest[i,j] =1 if sample i is in cluster j, 0 otherwise.
-Lbest: A list of the separate clusterings.  Lbest[[m]][i,j] = 1 if sample i is in cluster j for source m, 0 otherwise.
-AlphaVec: Vector of alpha values over MCMC draws, to assess mixing.  


For example, 
Clusts = BCC(list(X1,X2,X3),K=4) performs BCC with 4 clusters and 3 data sources, using the default settings. 
Then, Clusts$Cbest, gives overall clustering, Clusts$Lbest[[1]] gives clustering specific to X1, etc... 

