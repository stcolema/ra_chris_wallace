
BCC <- function(X,K=2,a=1,b=1,IndivAlpha = FALSE, mu0 = list(),a0=list(),b0=list(),Concentration=1,NumDraws = 1000){
#Initialize parameters
Gamma = rep(1,K) #Posterior Dirichlet concentration
Lbest = list() #List of best clustering for each data source
L = list() #Current clustering for each data source (updates each MCMC draw)
Llist = list() #saves all source clustering realizations
Clist = list() #saves all overall clustering realizations
mu = list() #list of mean vector for each data source (updates each draw)
logF = list() #log likeihood matrix for each data source
logL = list() #log clustering probabilities for each data source
Lp = list() #clustering probabilities for each data source
M = length(X) #Number of data sources
N = dim(X[[1]])[2] #sample sive
d = sapply(X,nrow) #Vector of dimension for each data source
S = list() #sample variance for each data source
A = list() #Posterior gamma shape parameter
B = list() #Posterior gamma rate parameter
Tau = list() #1/Tau = posterior variance
Sigma = list() #Posterior variance
Cprobs = matrix(nrow= N,ncol= K) #overall clustering probabilities
for(m in 1:M){ 
	Llist[[m]] = list()
	S[[m]] =  matrix(nrow = d[m],ncol = K)
	A[[m]] =  matrix(nrow = d[m],ncol = K)
	B[[m]] = matrix(nrow = d[m],ncol = K)
	Tau[[m]] = matrix(nrow = d[m],ncol = K)
	Sigma[[m]] = matrix(nrow = d[m],ncol = K)
	StDev = apply(X[[m]],1,sd) #Determine a0, b0 based on overall sample variance
	a0[[m]] = rep(1,d[m])
	b0[[m]] = StDev^2
	if(d[m]>1) mu0[[m]] = rowMeans(X[[m]]) #Determine mu0 by the overall mean
	if(d[m]==1) mu0[[m]] = mean(X[[m]])
	InitL = kmeans(t(X[[m]]),K) ##Initialize via K-means clustering of each data source
	if(m==1) prevClusters = InitL$cluster
	if(m>1){InitL$cluster = AlignClusters(prevClusters,InitL$cluster) ###Quick, imperfect step that helps to align cluster indices
		prevClusters = InitL$cluster}
	mu[[m]] = matrix(nrow = d[m],ncol= K)
	for(k in 1:K){
		if(d[m]>1){ Sigma[[m]][,k] = apply(X[[m]][,InitL$cluster==k],MARGIN=1,FUN='sd')
			mu[[m]][,k] = apply(X[[m]][,InitL$cluster==k],MARGIN=1,FUN='mean')}
		if(d[m]==1){ Sigma[[m]][,k] = sd(X[[m]][,InitL$cluster==k])
			mu[[m]][,k] = mean(X[[m]][,InitL$cluster==k])} }
	logF[[m]] = matrix(nrow = N, ncol = K)
	logL[[m]] = matrix(nrow = N, ncol = K)
	L[[m]] = matrix(nrow=N,ncol=K)
}
alphaVec = c() #Vector of alpha (adherence) parameters
C = matrix(nrow=N,ncol=K) #overall clustering
nu = array(rep(1/K,N*M*K),dim = c(N,M,K)) #source-specific cluster probabiities (given C)
Pi = rep(1/K,K) #cluster probabilites
n = matrix(nrow = M,ncol = K)#size of each source-specific cluster
alpha = rep(0,M) #initialize alpha
if(IndivAlpha)for(m in 1:M){ while(alpha[m] < 1/K) alpha[m] = rbeta(1,a,b) }
if(!IndivAlpha) while(alpha[1]<= 1/K) alpha[] = rbeta(1,a,b)

for(w in 1:NumDraws) { #w is the current MCMC iteration 
for(m in 1:M){
	for(k in 1:K)
  		logF[[m]][,k] = log(nu[,m,k])-sum(log(Sigma[[m]][,k]))+(-d[m]*log(2*pi)-colSums((((X[[m]]-mu[[m]][,k])/Sigma[[m]][,k])^2)))/2 #Determine logF from normal density
  logL[[m]] = logF[[m]]-apply(logF[[m]],1,logSum) #Normalize
}
	Lp = lapply(logL, exp)  
for(m in 1:M){	
	for(i in 1:N) L[[m]][i,] = rmultinom(1,1,Lp[[m]][i,]) #Generate L from Lp 
	if(w>1) L[[m]] = AlignClusters(C,L[[m]], type = 'mat') #Helps to align indices
	n[m,] = colSums(L[[m]])
	for(k in 1:K){ ###Update cluster parameters based on normal-gamma distribution
		if(d[m]==1&n[m,k]>1){
			S[[m]][,k] = sd(X[[m]][,L[[m]][,k]==1])^2
			PostMean = sum(X[[m]][,L[[m]][,k]==1])/(n[m,k]+1)
			B[[m]][,k] = b0[[m]]+0.5*(n[m,k]*S[[m]][,k]+n[m,k]*(mean(X[[m]][,L[[m]][,k]==1])-mu0[[m]])^2/(1+n[m,k]))}
		if(d[m]>1&n[m,k]>1){
			PostMean = (mu0[[m]]+rowSums(X[[m]][,L[[m]][,k]==1]))/(n[m,k]+1)
			S[[m]][,k] = apply(X[[m]][,L[[m]][,k]==1],MARGIN=1,FUN='sd')^2
	    B[[m]][,k] = b0[[m]]+0.5*(n[m,k]*S[[m]][,k]+n[m,k]*(rowMeans(X[[m]][,L[[m]][,k]==1])-mu0[[m]])^2/(1+n[m,k]))}
		if(n[m,k]==1){
			PostMean = (mu0[[m]]+X[[m]][,L[[m]][,k]==1])/2
		B[[m]][,k] = b0[[m]]+0.5*(X[[m]][,L[[m]][,k]==1]-mu0[[m]])^2/2}
		if(n[m,k]==0){
			PostMean = mu0[[m]]
			B[[m]][,k] = b0[[m]]}
		Lambda = 1+n[m,k]
		A[[m]][,k] = a0[[m]]+n[m,k]/2
		Tau[[m]][,k] = rgamma(d[m],shape=A[[m]][,k],rate=B[[m]][,k])
		mu[[m]][,k] = rnorm(d[m],PostMean,sqrt(1/(Tau[[m]][,k]*Lambda)))
		Sigma[[m]][,k] = sqrt(1/Tau[[m]][,k])}
}

if(w>1){ ###Updatae alpha
alpha = rep(1/K,M)
if(!IndivAlpha){
	NumEq = 0;
	for(m in 1:M) NumEq = NumEq+sum(L[[m]]==1 & C==1)
	for(Count in 1:10){ alphaTemp = rbeta(1,a+NumEq,b+M*N-NumEq)
#generate from beta dist until result >1/K (or set to 1/K after 10 trys)
		if(alphaTemp > 1/K){
 			alpha[] = alphaTemp
		break} } }
if(IndivAlpha){
	alpha = rep(1/K,M) 
	for(m in 1:M){
	 NumEq = sum(L[[m]]==1 & C==1)
	for(Count in 1:10){ alphaTemp = rbeta(1,a+NumEq,b+N-NumEq)
#generate from beta dist until result >1/K (or set to 1/K after 10 trys)
		if(alphaTemp > 1/K){ 
			alpha[m] = alphaTemp
			break; }
}}}
}
##Update C
for(k in 1:K){
	Cprobs[,k] = Pi[k]
	for(m in 1:M){
	Cprobs[,k] = Cprobs[,k]*alpha[m]^(L[[m]][,k])*((1-alpha[m])/(K-1))^(1-L[[m]][,k])}
}
Denom = rowSums(Cprobs)
for(i in 1:N){
	Cprobs[i,] = Cprobs[i,]/Denom[i]
C[i,] = rmultinom(1,1,Cprobs[i,])}
#update Pi
Gamma = Concentration+colSums(C)
Pi = rdirichlet(1,Gamma)

if(IndivAlpha) alphaVec = cbind(alphaVec,alpha)
if(!IndivAlpha) alphaVec[w] = alpha[1]

#update tau
for(m in 1:M){
	for(k in 1:K){
		nu[as.logical(C[,k]),m,k] = alpha[m]
	nu[!as.logical(C[,k]),m,k] = (1-alpha[m])/(K-1)
}}

for (m in 1:M)	Llist[[m]][[w]] = L[[m]]
Clist[[w]] = C
}
 ##Find Alpha by averaging over iteration (and 95% cred interval)
if(IndivAlpha){
AlphaBounds = matrix(nrow = length(alpha),ncol = 2)
for(m in 1:M){
	AlphaBounds[m,1] = quantile(alphaVec[m,floor(NumDraws/5):NumDraws],.025)
	AlphaBounds[m,2] = quantile(alphaVec[m,floor(NumDraws/5):NumDraws],.975)}
Alpha = rowMeans(alphaVec[,floor(NumDraws/5):NumDraws])
}
if(!IndivAlpha){
AlphaBounds = c()
AlphaBounds[1] = quantile(alphaVec[floor(NumDraws/5):NumDraws],.025)
AlphaBounds[2] = quantile(alphaVec[floor(NumDraws/5):NumDraws],.975)
Alpha = mean(alphaVec[floor(NumDraws/5):NumDraws])
}

#Choose hard clustering by least squares as in (Dahl,2006)
for(m in 1:M){
	Lkern = Llist[[m]][[floor(NumDraws/5)]]%*%t(Llist[[m]][[floor(NumDraws/5)]])
for(w in floor(NumDraws/5+1):NumDraws) Lkern = Lkern + Llist[[m]][[w]]%*%t(Llist[[m]][[w]])
Lkern = Lkern/(NumDraws-floor(NumDraws/5)+1)
CountLbest = N^2+1
for(w in floor(NumDraws/5):NumDraws){
	CountL = norm(Lkern-Llist[[m]][[w]]%*%t(Llist[[m]][[w]]), 'F')^2
	if(CountL<CountLbest){
		Lbest[[m]] = Llist[[m]][[w]]
	    CountLbest = CountL
	}		
}}

#Choose overall clustering
Ckern = Clist[[floor(NumDraws/5)]]%*%t(Clist[[floor(NumDraws/5)]])
for(w in floor(NumDraws/5+1):NumDraws) Ckern = Ckern + Clist[[w]]%*%t(Clist[[w]])
Ckern = Ckern/(NumDraws-floor(NumDraws/5)+1)
CountCbest = N^2+1
for(w in floor(NumDraws/5):NumDraws){
CountC = norm(Ckern-Clist[[w]]%*%t(Clist[[w]]), 'F')^2
	if(CountC<CountCbest){
		Cbest = Clist[[w]]
	    CountCbest = CountC
	}
}
return(list(Alpha = Alpha,AlphaBounds=AlphaBounds, Cbest = Cbest, Lbest = Lbest,AlphaVec = alphaVec))
}

logSum <- function(l){max(l) + log(sum(exp(l-max(l))))}

AlignClusters <- function(Z1,Z2, type = 'vec'){
if(type == 'vec'){
	for(k in 1:length(unique(Z1))){
		Max = sum(Z1==k &Z2==k)/(.01+sum(Z2==k)+sum(Z1==k));
		for(tempk in  1:length(unique(Z2))){
		if(sum(Z1==k &Z2==tempk)/(.01+sum(Z2==tempk)+sum(Z1==k)) > Max){
			Max = sum(Z1==k &Z2==tempk)/(.01+sum(Z2==tempk)+sum(Z1==k))
		dummy = Z2==k
		Z2[Z2==tempk] = k
		Z2[dummy] = tempk
}}}}
if(type == 'mat'){
	for(k in 1:dim(Z1)[2]){
		for(tempk in  1:dim(Z2)[2]){
		Max = sum(Z1==Z2)
		Z2dummy = Z2
		Z2dummy[,k] = Z2[,tempk]
		Z2dummy[,tempk] = Z2[,k]
		if(sum(Z1==Z2dummy)>Max) Z2 = Z2dummy
}}}
return(Z2)
}