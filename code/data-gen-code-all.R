# Generate y-X data 
# M replicates are created using Monte Carlo replication 
# M replicates for homoskedastic errors case and 
# M replicates for heteroskedastic errors case
# Once generated, the design matrix is kept fixed for all replications
# the design matrix columns are centered and have norm 1
# two types of errors are considered: N(0,1) and centered and scaled chi.square(df = 1) errors
# M replicates are created for both error types, for both homoskedastic and heteroskedastic cases
# the node-wise Lasso is run for the covariate matrix (for future Debiased Lasso computations)
# and the resulting quantity is stored for future Monte Carlo replications
# see van de Geer et al (2014, AoS) for details on debiased lasso


# two choices of (n,p): (n,p) = (150, 500) and (300, 500)
# p_0 = # of non-zero coefficients in true.beta.0 = 10 for both cases
# M = 500

# two different signal strengths are considered for (n = 150, p = 500) case
# this is particularly required for the heteroskedastic case, for use with PB
# in case (a): true.beta[j] = 5
# in case (b): true.beta[j] = 10

# for (n,p) = (300,500), we use true.beta[j] = 5, for both homoskedastic and heteroskedastic

rm(list = ls())
library(parallel)
library(mvtnorm)
library(hdi)

# generates error variables
error.fn = function(i, n, err.type)
{
	if(err.type == 1){
		err.vec = rnorm(n, sd = 1)
	}
	else{
		err.vec = (rchisq(n, df = 1) - 1)/sqrt(2) 
	}
	return(err.vec) 
}

# each data set (for the m-th Monte Carlo replication) are kept in a matrix format
# each is of dimension n x (1+p), where the first column contains y observations
# the irrepresentible condition for Lasso is verified (for the underlying design matrix)

yx.gen = function(n, beta.0.true, p0, rho, M)
{
	p = length(beta.0.true)	
	
	# Sigma is the cov-mat of a multivariate normal from which rows of X are generated #
	# elements are of the form sigma[i,j] = rho^(|i-j|).
	Sigma = matrix(0, nrow = p, ncol = p)
	for(i in 1:p)
	{
		Sigma[i, ] = (rho)^(abs(1:p - i))	
	}

	# design matrix columns centered and scaled, this matrix is not changed through the Monte-Carlo replications
	# this is because we are working in a non-stochastic X setup
	
	X.prelim.1 = rmvnorm(n = n, mean = rep(0, p), sigma = Sigma) 
	X.prelim.2 = apply(X.prelim.1, 2, function(u){u-mean(u)}) # columns are centered 
	X.fix = apply(X.prelim.2, 2, function(u){u/sqrt(sum(u^2))}) # columns are scaled to have norm 1 
	
	signal.true = X.fix%*%beta.0.true	
	
	# checking irrepresentible condition for Lasso #
	C.mat = t(X.fix)%*%X.fix/n 
	C11.mat <- C.mat[1:p0, 1:p0]
	C21.mat <- C.mat[-(1:p0), 1:p0]
	
	# because, in our setup all non-zero beta.j's are positive, i.e., sign(beta.j) = + 1, for j = 1,..,p_0
	a = as.vector(C21.mat%*%solve(C11.mat)%*%rep(1,p0)) 
	
	
	# heteroskedastic error variance generation
	# heterogenous variance terms depend only first two covariate values X_1 and X_2 
	# we use sigma.square = 3*sqrt(abs(X.fix[,1])) + abs(X.fix[,2])) + 0.2 
	  
	sigma.sqr.het <- 3*sqrt(abs(X.fix[,1]) + abs(X.fix[,2])) + 0.2 # this is a n-vector 
	
	# with above choice of X.fix, 
	# for n = 150: this fixes the sigma.sqr.het range between roughly (0.5, 2)
	# for n = 300: the range becomes roughly (0.3, 2)
		
	# y.matrix for both types of data are generated, M Monte Carlo replications are used	
	y.array <- y.array.het <- array(0, dim = c(2, n, M)) 
	for(j in 1:2)
	{
		v.j <- sapply(1:M, error.fn, n = n, err.type = j) # n x M matrix of errors
		
		# the homoskedastic errors are scaled by the changing standard deviations in the heteroskedastic case
		v.j.het <- apply(v.j, 2, function(e.vec, scale.fac){e.vec*scale.fac}, scale.fac = sqrt(sigma.sqr.het)) 
		
		y.array[j, , ] = as.vector(signal.true) + v.j 
		y.array.het[j, , ] = as.vector(signal.true) + v.j.het 
		# different errors are added to the same true signal 
	}
	yx.list = list() # list storing all relevant data used for Monte Carlo simulation
	yx.list[[1]] <- y.array # homogenous errors y data
	
	yx.list[[2]] <- X.fix # X matrix
	yx.list[[3]] <- beta.0.true # 
	
	# irrepresentible condition value is stored #
	yx.list[[4]] <- c(min(a), max(a)) 
	
	
	# finding Z - see hdi manual # used for debiased Lasso stage #
	# this does not require use of y values
	# each column of X is predicted by Lasso using other columns of X 
	# the resulting output is saved, for future use in Monte Carlo simulation stage
	A = lasso.proj(X.fix, y.array[1, ,1], parallel = TRUE, ncores = 12, return.Z = TRUE, suppress.grouptesting = TRUE)
	yx.list[[5]] = A$Z
	
	yx.list[[6]] = y.array.het # heterogenous errors y data
	
	my.file.name = paste("yx-all-n-",n,"-p-",p,"-true-beta-",beta.0.true[1],".Rdata", sep = "")
	save(yx.list, file = my.file.name)
}
beta.0.true = c(rep(5,10), rep(0,490)) 
yx.gen(n = 150, beta.0.true, p0 = 10, rho = 0.6, M = 500)
