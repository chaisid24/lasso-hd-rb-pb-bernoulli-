# description of the code
# this code deals with the homogenous errors case
# computes empirical coverages (EC) and average lengths (ALen) of confidence intervals (CIs) for a target parameter
# the target parameter is a true regression coefficient beta_{j,0}, for j = 1,...,p_0 ( = 10)
# in our simulations we focused on two different targets: beta.0[1] and beta.0[10]
# the EC and ALen are found at various choices of lambda = k*lambda_{CV}
# where k > 0, is a constant factor that is changed and
# lambda_{CV} is found by 5 fold CV
# the EC and ALen are found/plotted over changing values of k

# the code is run by using (see at the bottom of this code)
# lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 750, err.type = 1, k.fac = k0.seq)
# where, n = sample size, p = number of variables, CV.k = 5 (for 5 fold CV)
# coef.index = 1 (or 10), depending on whether we target beta_{1,0} (or beta_{10,0}) (this can be changed)
# alpha = 0.1, for 90% nominal coverage levels
# B.boot = 750 = no. of bootstrap iterations within each Monte Carlo replication
# err.type = 1 (or 2), for using the N(0,1) errors data set (or using the chi.square(1) errors data set)
# k.fac = a vector of positive numbers (in decreasing order) for use in place of k in k*lambda_{CV}
# we used k.fac = k0.seq = seq(4, 0.375, by = -0.125)

rm(list = ls())
library(parallel)
library(glmnet)
library(hdi)


# this code computes lasso estimates using glmnet package
# also computes lambda_{CV} for each Monte Carlo data-set
lasso.fn = function(y, X.scale, CV.list) # use this code when cv.lam.search is required 
{
	# y is a n-vector, a single response vector for the m-th Monte Carlo replication
	
	n = nrow(X.scale) # X.scale is the design matrix
	p = ncol(X.scale)
	
	CV.k = CV.list[[1]] # the k-fold parameter (in our simulation it is always equal to 5) 
    	k.fac = CV.list[[2]] # this is the k.fac vector of length R

	# finds the optimal lambda_{CV} value 
    	CV.lam = cv.glmnet(X.scale, y, nfolds = CV.k, intercept = F, standardize = F)$lambda.min
    	lam1 = k.fac*CV.lam 									# creates the k.fac*lambda_{CV} vector 
    	A1 = glmnet(X.scale, y, intercept = F, standardize = F, lambda = lam1) 			# computes Lasso solution at these lambda values
    	beta.lasso.1 = A1$beta 									# puts these Lasso estimates in a pxR matrix format
    
  	out.all = list(beta.lasso.1, lam1)							# returning the Lasso solutions and lam.1 values
	return(out.all)
}

# computes Lasso solutions using glmnet, but at a fixed choice of lambda, used in the bootstrap stage where lambda is provided
lasso.within.boot <- function(b, y.boot.mat, X, lam.fix)
{
	y.star = y.boot.mat[ ,b] 							# b-th bootstrapped response vector (either RB or PB)
    	A = glmnet(X, y.star, intercept = F, standardize = F, lambda = lam.fix)		# lam.fix is a scalar
	beta.star = as.vector(A$beta)							# a p-vector, bootstrapped Lasso estimate
	return(beta.star)				
}

# intermediate step in RB based CI computation
# works on the b-th bootstrapped data set
# d.vec is a p-vector with all components = 0, but d.vec[coef.index] <- 1
# beta.plols is the post Lasso ordinary least squares (OLS) estimate
rb.mid.step.vareq = function(b, y.star.mat, X, beta.star.mat, lam, d.vec, beta.plols)
{
    	y.star = y.star.mat[ ,b]
    	beta.star = beta.star.mat[ ,b]
        
	n = nrow(X)
    	p = ncol(X)
    	T.star = sqrt(n)*sum(d.vec*(beta.star - beta.plols)) 	# this finds sqrt(n)*(beta.star[coef.index] - beta.plols[coef.index])
    
    	beta.dot.star = rep(0, p)				# initial definition

	# constructing relevant quantities for computing the RB based pivotal quantity
	A.star = (1:p)[beta.star!=0]
	C11.star = t(X[ ,A.star])%*%(X[ ,A.star])/n
    	c.star = sign(beta.star[A.star])
    	term.1 = as.vector(matrix(d.vec[A.star], nrow = 1)%*%solve(C11.star)%*%matrix(c.star, ncol = 1))
    	b0.star = -sqrt(n)*lam*term.1
    
	beta.dot.star[A.star] = beta.star[A.star] + lam*as.vector(solve(C11.star)%*%matrix(c.star, ncol = 1))
    	sigma.star = sqrt(mean((y.star - X%*%beta.dot.star)^2))
    	R.check.star = (T.star - b0.star)/sqrt(sigma.star)	# bias corrected and studentized RB based pivot
    	return(R.check.star)
}

# main RB related code for computing relavent CIs and if they contain the true parameter or not and the CI lengths
# nz.index is a select set of values from 1:R, where R = length of k.fac vector (see details below)
# beta.0 is the true regression coefficient
# init.info.list contains objects needed to run this code
# this code provides RB related output for a single Monte-Carlo iteration
rbci.vareq.fn = function(nz.index, y, X, coef.index, B, alpha, beta.0, init.info.list)
{
    	n = nrow(X)
	p = ncol(X)
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	# init.info.list[[1]] was originally a p x R matrix, where R is described above  
	# containing p-dimensional Lasso estimates for all R choice of lambda = k.fac*lambda_{CV}
	beta.hat = init.info.list[[1]][ ,nz.index]	# a specific column of Lasso estimates for a specific lambda   	
	lam = init.info.list[[2]][nz.index]       	# that specific lambda value 
	theta.hat = init.info.list[[3]][nz.index] 	# target parameter estimate extracted from beta.hat; in our simulation it is simply beta.hat[coef.index] 
	theta.0 = init.info.list[[4]]			# true target parameter; in our simulation it is beta.0[coef.index]
    	rboot.imat = init.info.list[[5]]             	# nxB matrix of sampled indices

	# RB related preiliminary computations 
    	A.hat = (1:p)[beta.hat!=0]
	d.hat = d.vec[A.hat]
	C11.hat = t(X[ ,A.hat])%*%(X[ ,A.hat])/n
	beta.plols = rep(0, p)
	beta.plols[A.hat] = as.vector(lm(y~-1+X[ ,A.hat])$coefficients)
	res.check = y - X%*%beta.plols
	sigma.check.sqr = mean(res.check^2)
	s.vec = sign(beta.hat[A.hat])
	term.1 = as.vector(matrix(d.hat, nrow = 1)%*%solve(C11.hat)%*%s.vec)
	b0 = -sqrt(n)*lam*term.1
	c.res.check = res.check - mean(res.check)
	
	# RB part #
	y.boot.rb = as.vector(X%*%beta.plols) + sapply(1:B, function(i, b1.mat, c.res){c.res[b1.mat[ ,i]]}, b1.mat = rboot.imat, c.res = c.res.check) 
	beta.star.mat = sapply(1:B, lasso.within.boot, y.boot.rb, X, lam) # pxB mat # 
	R.check.star = sapply(1:B, rb.mid.step.vareq, y.boot.rb, X, beta.star.mat, lam, d.vec, beta.plols) # B vector #
	
    	quant.R.check.star = quantile(R.check.star, probs = c(alpha/2, alpha, 1-alpha, 1 - alpha/2))
    	quant.R.check.abs.star = quantile(abs(R.check.star), probs = 1-alpha)

    	lower.rbci.ep = theta.hat - (b0 + sqrt(sigma.check.sqr)*quant.R.check.star[3])/sqrt(n)
    	upper.rbci.ep = theta.hat - (b0 + sqrt(sigma.check.sqr)*quant.R.check.star[2])/sqrt(n)
    	two.rbci.ep = rep(theta.hat - b0/sqrt(n), 2) - quant.R.check.star[c(4, 1)]*sqrt(sigma.check.sqr)/sqrt(n)
    	two.symm.rbci.ep = rep(theta.hat - b0/sqrt(n), 2) + c(-1, 1)*quant.R.check.abs.star*sqrt(sigma.check.sqr)/sqrt(n)

    	logic.L = ifelse(lower.rbci.ep <= theta.0, 1, 0)
    	logic.U = ifelse(theta.0 <= upper.rbci.ep, 1, 0)
    	logic.2 = ifelse(two.rbci.ep[1] <= theta.0 & theta.0 <= two.rbci.ep[2], 1, 0)
    	logic.2s = ifelse(two.symm.rbci.ep[1] <= theta.0 & theta.0 <= two.symm.rbci.ep[2], 1, 0)
    	len.2 = abs(two.rbci.ep[1] - two.rbci.ep[2])
    	len.2s = abs(two.symm.rbci.ep[1] - two.symm.rbci.ep[2])

	# RB related output for one Monte-Carlo dataset 
    	rb.ci.out.nz.index = c(logic.L, logic.U, logic.2, logic.2s, len.2, len.2s)  # vector of length 6
    	return(rb.ci.out.nz.index)
}

# intermediate steps for PB based CI computation
pb.mid.step = function(b, beta.ss.mat, X, y, lam, d.vec, beta.plols, sigma.check, Sigma.tilde)
{
    	n = nrow(X)
    	p = ncol(X)
    
    	beta.ss = beta.ss.mat[ ,b]
    	A.ss = (1:p)[beta.ss!=0]
    	d.ss = d.vec[A.ss]
    
    	C11.ss = t(X[ , A.ss])%*%(X[ , A.ss])/n
    	s.vec.ss = sign(beta.ss[A.ss])
    	term.2 = as.vector(lam*solve(C11.ss)%*%s.vec.ss)
    	b0.ss = -sqrt(n)*lam*as.vector(matrix(d.ss, nrow = 1)%*%solve(C11.ss)%*%matrix(s.vec.ss, ncol = 1))
    
    	beta.dot.ss = rep(0, p)
    	beta.dot.ss[A.ss] = beta.ss[A.ss] + term.2
    	sigma.ss.sqr = mean((y - X%*%beta.dot.ss)^2)
    
	T.ss = sqrt(n)*sum(d.vec*(beta.ss - beta.plols))
    	Rt.ss = (1/sqrt(sigma.ss.sqr))*sigma.check*(Sigma.tilde^(-1/2))*(T.ss - b0.ss) 
	return(Rt.ss)
}

# this code plays a similar roles as the rbci.vareq.fn(.) code, expect for the case of PB
# see rbci.vareq.fn(.) and the explanations/comments within that code
pbci.vareq.fn = function(nz.index, y, X, coef.index, B, alpha, beta.0, init.info.list)
{
    	n = nrow(X)
	p = ncol(X)
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.hat = init.info.list[[1]][ ,nz.index]
	lam = init.info.list[[2]][nz.index]
	theta.hat = init.info.list[[3]][nz.index]
	theta.0 = init.info.list[[4]]
    	pboot.imat = init.info.list[[6]]             # nxB matrix of indices 
   
        A.hat = (1:p)[beta.hat!=0]
	d.hat = d.vec[A.hat]
	C11.hat = t(X[ ,A.hat])%*%(X[ ,A.hat])/n
	s.vec = sign(beta.hat[A.hat])
	term.1 = d.hat%*%solve(C11.hat) # this is a 1 x p0.hat matrix
	
	beta.plols = rep(0, p)
	beta.plols[A.hat] = as.vector(lm(y~-1+X[ ,A.hat])$coefficients)
	res.check = y - X%*%beta.plols
	sigma.check.sqr = mean(res.check^2)
	Sigma.hat = as.vector(term.1%*%matrix(d.hat, ncol = 1))
	b0 = -sqrt(n)*lam*as.vector(term.1%*%matrix(s.vec, ncol = 1))
	xi0.vec = as.vector(term.1%*%t(X[ ,A.hat])) # n vector 
	Sigma.tilde = mean((xi0.vec^2)*(res.check^2))

	# PB part #
	z.mat = as.vector(X%*%beta.plols) + sapply(1:B, function(b, u, r){r*u[ ,b]}, u = pboot.imat, r = res.check)  # nxB matrix
        beta.ss.mat = sapply(1:B, lasso.within.boot, z.mat, X, lam) # pxB matrix #
        Rt.ss.vec = sapply(1:B, pb.mid.step, beta.ss.mat, X, y, lam, d.vec, beta.plols, sqrt(sigma.check.sqr), Sigma.tilde) # B-vec
        quant.Rt.ss = quantile(Rt.ss.vec, probs = c(alpha, 1-alpha, alpha/2, 1-alpha/2))
    
        # Symmetric PB calculations #
    	z.alpha = qnorm(1 - alpha/2)
    	w2 = -mean((xi0.vec^2)*(res.check^4))/(sigma.check.sqr*Sigma.tilde) + mean(res.check^4)/(sigma.check.sqr^2)
    	w4 = 2*mean((xi0.vec^4)*(res.check^4))/(Sigma.tilde^2) + (4/(sigma.check.sqr*Sigma.tilde))*mean((xi0.vec^2)*(res.check^4)) - (3/(sigma.check.sqr^2))*mean(res.check^4) + 1
    	Cn.p = -(z.alpha/n)*(w2/2 + w4*(z.alpha^2 - 3)/(24))
    	h.dag = quantile(abs(Rt.ss.vec), probs = 1-alpha)
    	h.check = h.dag + Cn.p
    
    	# CI construction 
	ci.ep.low = theta.hat - b0/sqrt(n) - sqrt(sigma.check.sqr*Sigma.hat/n)*quant.Rt.ss[2]
    	ci.ep.upp = theta.hat - b0/sqrt(n) - sqrt(sigma.check.sqr*Sigma.hat/n)*quant.Rt.ss[1]
    	logic.L = ifelse(ci.ep.low <= theta.0, 1, 0)
    	logic.U = ifelse(theta.0 <= ci.ep.upp, 1, 0)
    
    	ci.ep.2 = rep(theta.hat - b0/sqrt(n), 2) - sqrt(sigma.check.sqr*Sigma.hat/n)*quant.Rt.ss[c(4,3)]
    	logic.2 = ifelse(ci.ep.2[1] <= theta.0 & theta.0 <= ci.ep.2[2], 1, 0) 
    	len.2 = abs(ci.ep.2[2] - ci.ep.2[1])
    
    	ci.ep.2s = rep(theta.hat - b0/sqrt(n), 2) + c(-1,1)*sqrt(sigma.check.sqr*Sigma.hat/n)*h.check
    	logic.2s = ifelse(ci.ep.2s[1] <= theta.0 & theta.0 <= ci.ep.2s[2], 1, 0) 
    	len.2s = abs(ci.ep.2s[2] -ci.ep.2s[1])

    	pb.ci.out.nz.index = c(logic.L, logic.U, logic.2, logic.2s, len.2, len.2s)  # vector of length 6
    	return(pb.ci.out.nz.index)
}

# this code finds Debiased Lasso based CIs for a given (y,X) utiliing the A$Z information stored earlier in the data generation phase
dblasso.ci = function(y, X, coef.index, alpha, A.Z)
{
    	A = lasso.proj(X, y, Z = A.Z, suppress.grouptesting = TRUE)
    	u1 = confint(A, level = 1-alpha)
    	return(as.vector(u1[coef.index, ]))
}

# code that handles a single Monte-Carlo iteration (m-th iteration)
data.set.m = function(m, y.mat, X, CV.k, coef.index, B, alpha, beta.0, k.fac, A.Z)
{
	y = y.mat[ ,m]			# m-th Monte Carlo response vector
	n = nrow(X)			
	p = ncol(X)
	R = length(k.fac)		# length of k.fac vector
	p0 = length(beta.0[beta.0!=0])	# no. of non-zero components in beta.0.true vector (here p_0 = 10)
	
	d.vec = numeric(p)
	d.vec[coef.index] <- 1

	C = lasso.fn(y, X, list(CV.k, k.fac))  		# finds a pxR matrix of Lasso estimates as well as the lambda = k.fac*lambda_{CV} vector of length R
	beta.mat = C[[1]] 		       		# a pxR matrix of Lasso estimates
    	lam = C[[2]] 			       		# R vector of lambda = k.fac*lambda_{CV} values 
	theta.vec = as.vector(d.vec%*%beta.mat) 	# R-vector of target estimates; in our case these are simply beta.mat[coef.index, ] values

	# those index (=r) values, for which theta.vec[r] is non-zero 
	# as we are targeting inference for underlying non-zero coefficients, 
	# the simulation is able to proceed only if a Lasso estimate for a target parameter turns out to be non-zero
	# essentially inz.theta identifies those lambda values which give non-zero Lasso estimates of the target theta.0 = beta.0[coef.index]
	inz.theta = (1:R)[theta.vec!=0] 		
	theta.0 = sum(d.vec*beta.0) # true parameter 

	## Debiased Lasso based CIs ##
	# The Debiased Lasso based CI's are always usable 
	dbl.ci = dblasso.ci(y, X, coef.index, alpha, A.Z)
    	dbl.logic = ifelse(dbl.ci[1]<= theta.0 & theta.0 <= dbl.ci[2], 1, 0)
	dbl.len = abs(dbl.ci[2] - dbl.ci[1])
	dlasso.out.m <- c(dbl.logic, dbl.len) 
	

	## Oracle (Normal), RB and PB based CIs ##
	# simulation proceeds if among all R choices of lambda = k.fac*lambda_{CV} there is at least one r, which gave theta.hat[r] non-zero
	# infact we drop all r, for which theta.vec[r] = 0
	# effectively to find empirical coverages and average lengths at a certain choice of k.fac[r], we divide not by M = 500 (no. of Monte-Carlo replications)
	# but, by M.hat[r] = number of Monte-Carlo replications at which lambda[r] = k.fac[r]*lambda_{CV} gave a non-zero Lasso estimate of the target parameter
	# these M.hat[r] values for r = 1,...,R, are available in our final output
	# for smaller lambda[r] values, M.hat[r] is close to M = 500, but for larger values of lambda[r], the M.hat[r] count decreases
	# dividing by M.hat[r] makes sense, as it is only in those replications, where we are identifying the target parameter as a non-zero value
	
	
	if(length(inz.theta)>0)
	{
		# print(c(m, 1)) # an indication that at least some theta.hat[r] is non-zero
		
		# Creating bootstrap indices #
		rboot.imat = sapply(1:B, function(i, z){sample(z, size = length(z), replace = T)}, z = 1:n)
		mu.star = (1/2)/(1/2 + 3/2)
		G.mat = sapply(1:B, function(i, n.size){rbeta(n.size, 1/2, 3/2)}, n.size = n) # n x B matrix 
		pboot.imat = apply(G.mat, 2, function(v2, c1){(1/c1)*(v2-c1)}, c1 = mu.star) # n x B transformed from G.mat#
    
		init.info.list = list(beta.mat, lam, theta.vec, theta.0, rboot.imat, pboot.imat)
	
		# Oracle CI 
		# Approximates Lasso by a Normal limit law assuming true non-zero coefficient indices (relevant variables) are known 
    		oracle.out.m <- matrix(0, nrow = R, ncol = 2) # allow R rows to allow the possibility that all theta.vec[r] are non-zero
		C11 = t(X[ ,1:p0])%*%(X[ ,1:p0])/n
		C11.inv = solve(C11)
		d1.vec = d.vec[1:p0] # p0 dimensional vector #
		rho.sqr = as.numeric(matrix(d.vec[1:p0], nrow = 1)%*%C11.inv%*%matrix(d.vec[1:p0], ncol = 1))
		sigma.check = numeric(R)
	
		for(r in inz.theta)
		{
		    sigma.check[r] = sqrt(mean((y - X%*%beta.mat[ ,r])^2))
		    oracle.ci.r = rep(theta.vec[r], 2) - qnorm(c(1-alpha/2, alpha/2))*sigma.check[r]*sqrt(rho.sqr)/sqrt(n)
		    oracle.out.m[r, 1] = ifelse(oracle.ci.r[1] <= theta.0 & theta.0 <= oracle.ci.r[2], 1, 0)
		    oracle.out.m[r, 2] = abs(oracle.ci.r[2] - oracle.ci.r[1])	      
		}
		
    		# RB and PB #
    		RB.out.mat.m <- PB.out.mat.m <- matrix(0, nrow = R, ncol = 6) # each is a Rx6 matrix
    		RB.out.mat.m[inz.theta, ] = t(sapply(inz.theta, rbci.vareq.fn, y, X, coef.index, B, alpha, beta.0, init.info.list)) # length(inz.theta)x6 matrix
    		PB.out.mat.m[inz.theta, ] = t(sapply(inz.theta, pbci.vareq.fn, y, X, coef.index, B, alpha, beta.0, init.info.list)) # length(inz.theta)x6 matrix

		# the first element in this list c(1,m) is used denote that at least one theta.vec[r], for r=1,...,R, is non-zero
		# essentially in the m-th Monte-Carlo iteration, we have atleast one choice of lambda[r], which gave a non-zero estimate of the target parameter
    		return(list(c(1,m), inz.theta, oracle.out.m, RB.out.mat.m, PB.out.mat.m, dlasso.out.m))
    	}
    	else 
    	{
		# c(2,m) is used to denote that in the m-th Monte-Carlo iteration, all theta.vec[r] = 0, for r = 1,...,R.
    		return(list(c(2,m), inz.theta, dlasso.out.m)) # only relevant output is the Debiased Lasso based output
    	}
	# processing of m-th Monte Carlo dataset finishes #
}


# processing of all Monte-Carlo datasets
# we supply the following
# n = sample size
# p = number of variables
# CV.k = number of folds (k) used for k-fold CV based choice of lambda_{CV}, in our case we always used CV.k = 5
# alpha = 0.1, for 90% nominal coverage 
# B.boot = number of bootstrap iterations within each Monte-Carlo replication (= 750, for our simulation)
# err.type = 1 or 2, depending on F = F_1 = N(0,1), or F = F_2 = chi.square(1) errors 
# k.fac = the R length vector of constants multiplied to lambda_{CV}

lasso.all <- function(n, p, CV.k, coef.index, alpha, B.boot, err.type, k.fac)
{
	my.file.name = paste("yx-all-n-",n,"-p-",p,"-true-beta-5.Rdata", sep = "")
	load(my.file.name) 	# loads the relavent data set, which has been generated and saved by the data generation code
		
	y.mat = yx.list[[1]][err.type, , ]	# loads n x M response matrix for homogenous-case, depending on err.type (1 or 2) 
	X.fix = yx.list[[2]]			# design matrix
	beta.0 = yx.list[[3]]			# true beta.0 parameter
	A.Z = yx.list[[5]]			# nodewise Lasso used for Debiased Lasso related computations
	M = ncol(y.mat)				# M = 500 for our simulation
	R = length(k.fac)			# length of k.fac vector

	# processing all Monte-Carlo datasets in parallel by calling the previous function
	# mc.cores can be changed
	# saves all longer/larger/initial output in a file prefixed with "long-n-..."
	B <- mclapply(1:M, data.set.m, mc.cores = 8, y.mat, X.fix, CV.k, coef.index, B.boot, alpha, beta.0, k.fac, A.Z)
	save(B, file = paste("long-n-",n,"-p-",p,"-coef-",coef.index,"-et-",err.type,".Rdata", sep = ""))
		
	
	# evaluating coverage probabilities 
	prelim.count <- matrix(0, nrow = M, ncol = 2) # an initial storage matrix for c(1,m) or c(2,m)

	# a R length vector, which will count for each r=1,..,R, the no. of Monte-Carlo data-sets where it provided theta.hat[r] not equal to zero
	count.vec <- numeric(R)	
	
	oracle.out.M <- matrix(0, nrow = R, ncol = 2)
	RB.out.mat.M <- PB.out.mat.M <- matrix(0, nrow = R, ncol = 6)
	dlasso.out.M <- rep(0, 2)

	# this part is used to count M.hat[r] values (see above)
	for(m in 1:M)
	{
		prelim.count[m, ] <- B[[m]][[1]]	# it is either c(1,m) or c(2,m)
		if(prelim.count[m, 1] == 1)		# if the first element is 1, then Oracle, RB and PB based CIs have been computed for m-th data-set
		{
			a.m <- B[[m]][[2]]				# identifying the inz.theta values for the m-th data set
			count.vec[a.m] <- count.vec[a.m] + 1		# increment those components of count.vec by +1
			oracle.out.M <- oracle.out.M + B[[m]][[3]]	# keep adding to the oracle.out.M matrix another Rx2 matrix
			RB.out.mat.M <- RB.out.mat.M + B[[m]][[4]]	# same as above, add a Rx6 matrix
			PB.out.mat.M <- PB.out.mat.M + B[[m]][[5]]	# ...
			dlasso.out.M <- dlasso.out.M + B[[m]][[6]]	# keep adding a 2-length vector
		}
		else
		{
			dlasso.out.M <- dlasso.out.M + B[[m]][[3]]	# only add for the dlasso part, a 2-length vector
		}
	}			
	
	for(r in 1:R)
	{
		if(count.vec[r]>0)	# coverages can only be found if count.vec[r] is positive
		{
			oracle.out.M[r, ] <- oracle.out.M[r, ]/count.vec[r]
			RB.out.mat.M[r, ] <- RB.out.mat.M[r, ]/count.vec[r]
			PB.out.mat.M[r, ] <- PB.out.mat.M[r, ]/count.vec[r]
		}
		else
		{
			oracle.out.M[r, ] <- RB.out.mat.M[r, ]<- PB.out.mat.M[r, ] <- 0
		}
	}
	dlasso.out.M <- matrix(dlasso.out.M/M, nrow = 1, ncol = 2)	# this is divided by M 
		
	result.mat <- cbind(k.fac, count.vec, oracle.out.M, RB.out.mat.M, PB.out.mat.M) # putting the main outputs together

	# prints the outputs
	colnames(result.mat) <- c("k.fac", "count.vec", "or.ec", "or.alen", "rbL", "rbU", "rb2", "rb2s", "rb2.alen", "rb2s.alen", "pbL", "pbU", "pb2", "pb2s", "pb2.alen", "pb2s.alen")
	print(signif(result.mat, 4))
	colnames(dlasso.out.M) <- c("dlasso.ec", "dlasso.alen")
	print(dlasso.out.M)

	# saves all output data which are found after processing the earlier main output in a file prefixed by "short-n-..."
	C <- list(prelim.count, k.fac, count.vec, oracle.out.M, RB.out.mat.M, PB.out.mat.M, dlasso.out.M, result.mat)
	save(C, file = paste("short-n-",n,"-p-",p,"-coef-",coef.index,"-et-",err.type,".Rdata", sep = ""))
		
}

# declaring the k.fac squence of constants
k0.seq= seq(4, 0.375, by = -0.125)

# In order to run this code, the .Rdata files generated by the data generation code should be placed in the same directory
# Currently these data files are stored in the data folder in this github repository

########################
# RUNNING n = 300 case #
# uncomment as per requirement #

#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 750, err.type = 1, k.fac = k0.seq)
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 750, err.type = 1, k.fac = k0.seq)
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 750, err.type = 2, k.fac = k0.seq)
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 750, err.type = 2, k.fac = k0.seq)

########################
# RUNNING n = 150 case #
# uncomment as per requirement #

#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 750, err.type = 1, k.fac = k0.seq)
#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 750, err.type = 1, k.fac = k0.seq)
#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 750, err.type = 2, k.fac = k0.seq)
#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 750, err.type = 2, k.fac = k0.seq)
