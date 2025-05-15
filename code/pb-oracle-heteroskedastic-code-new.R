rm(list = ls())
library(parallel)
library(hdi)
library(glmnet)

# computes lasso estimator
# finds lam.CV = optimal K-fold CV based choice of lambda
# supplies lasso solutions at c*lam.CV, for a choice of c values
lasso.fn = function(y, X, CV.info)  
{
		n = nrow(X)
		p = ncol(X)
		
		CV.K = CV.info[[1]] 	# the K parameter in K-fold CV 
		c.vec = CV.info[[2]]	# vector of constants 'c' in decreasing order for using in c*lam.CV # length = R
    
		# optimum CV based lambda
        lam.CV = cv.glmnet(X, y, nfolds = CV.K, intercept = F, standardize = F)$lambda.min  
        lam.vec = c.vec*lam.CV	# R length vector of lambda values where lasso solutions will be computed 	
    
        A = glmnet(X, y, intercept = F, standardize = F, lambda = lam.vec)
        beta.lasso = A$beta 	# pxR matrix of lasso coefficients
    
		out.all = list(beta.lasso, lam.vec, lam.CV)								
		return(out.all)
}

# computes lasso solution for one set of y and X at a fixed lambda 
# used in the bootstrap stage, lam.fix is supplied
lasso.within.boot <- function(b, y.boot.mat, X, lam.fix)
{
    	y.star = y.boot.mat[ ,b]
    	A = glmnet(X, y.star, intercept = F, standardize = F, lambda = lam.fix)
		beta.star = as.vector(A$beta)
		return(beta.star)
}

pb.mid.step = function(b, beta.ss.mat, X, y, lam, d.vec, beta.plols, pboot.imat, xi.0.vec)
{
		n = nrow(X)
		p = ncol(X)
    
    	beta.ss = beta.ss.mat[ ,b]
    	A.ss = (1:p)[beta.ss!=0]
    	d.ss = d.vec[A.ss]
    
    	C11.ss = t(X[ , A.ss])%*%(X[ , A.ss])/n
    	s.vec.ss = sign(beta.ss[A.ss])
    	term.2 = as.vector(lam*solve(C11.ss)%*%s.vec.ss) # originally lam is lam/(2*n)
		b0.ss = -sqrt(n)*lam*as.vector(matrix(d.ss, nrow = 1)%*%solve(C11.ss)%*%matrix(s.vec.ss, ncol = 1))
		# for b0.ss, original scaling factor was -1/(2*sqrt(n)) instead of -sqrt(n)

    	beta.dot.ss = rep(0, p)
    	beta.dot.ss[A.ss] = beta.ss[A.ss] + term.2
    	v1 = (y - X%*%beta.dot.ss)^2
    	v2 = v1*(pboot.imat[ ,b])^2
		Sigma.tilde.ss = mean((xi.0.vec^2)*v2)
    
    	T.ss = sqrt(n)*sum(d.vec*(beta.ss - beta.plols))
    	R.dot.s = ((Sigma.tilde.ss)^(-1/2))*(T.ss - b0.ss) 
		return(R.dot.s)
}
pbci.fn = function(nz.index, y, X, coef.index, B, alpha, beta.0, init.info.list)
{
    	#cat(nz.index,"\n")
		n = nrow(X)
		p = ncol(X)
		d.vec = numeric(p)
		d.vec[coef.index] = 1 

		beta.hat = init.info.list[[1]][ ,nz.index]
		lam = init.info.list[[2]][nz.index]
		theta.hat = init.info.list[[3]][nz.index]
		theta.0 = init.info.list[[4]]
    	pboot.imat = init.info.list[[5]]             # nxB matrix 
   
    	A.hat = (1:p)[beta.hat!=0]
		d.hat = d.vec[A.hat]
		C11.hat = t(X[ ,A.hat])%*%(X[ ,A.hat])/n
		c.vec = sign(beta.hat[A.hat])
		term.1 = d.hat%*%solve(C11.hat) # this is a 1 x p0.hat matrix
	
		beta.plols = rep(0, p)
		beta.plols[A.hat] = as.vector(lm(y~-1+X[ ,A.hat])$coefficients)
		theta.plols = sum(d.vec*beta.plols)
		res.check = y - X%*%beta.plols
		xi0.vec = as.vector(term.1%*%t(X[ ,A.hat])) # n vector 
		Sigma.tilde = mean((xi0.vec^2)*(res.check^2))
		b0 = -sqrt(n)*lam*as.vector(term.1%*%matrix(c.vec, ncol = 1)) # theory says divide by 2*sqrt(n) # changed due to package

		#sigma.check.sqr = mean(res.check^2)
		#Sigma.hat = as.vector(term.1%*%matrix(d.hat, ncol = 1))

		# PB part #
		z.mat = as.vector(X%*%beta.plols) + sapply(1:B, function(b, u, r){r*u[ ,b]}, u = pboot.imat, r = res.check)  # nxB matrix
   		beta.ss.mat = sapply(1:B, lasso.within.boot, z.mat, X, lam) # pxB matrix #

        R.dot.star.vec = sapply(1:B, pb.mid.step, beta.ss.mat, X, y, lam, d.vec, beta.plols, pboot.imat, xi0.vec) # B-vec
    	quant.R.dot.star = quantile(R.dot.star.vec, probs = c(alpha, 1-alpha, alpha/2, 1-alpha/2))
    
		# Symmetric PB calculations #
    	z.alpha = qnorm(1 - alpha/2)
    	Cn.p = ((-1)*(z.alpha*(z.alpha^2-3))/12)*mean((xi0.vec^4)*(res.check^4))/(Sigma.tilde^2)
    	h.dag = quantile(abs(R.dot.star.vec), probs = 1-alpha)
    	h.check = h.dag + Cn.p
    
    	# CI construction #
		ci.ep.low = theta.hat - b0/sqrt(n) - sqrt(Sigma.tilde/n)*quant.R.dot.star[2]
    	ci.ep.upp = theta.hat - b0/sqrt(n) - sqrt(Sigma.tilde/n)*quant.R.dot.star[1]
    	logic.L = ifelse(ci.ep.low <= theta.0, 1, 0)
    	logic.U = ifelse(theta.0 <= ci.ep.upp, 1, 0)
    
    	ci.ep.2 = rep(theta.hat - b0/sqrt(n), 2) - sqrt(Sigma.tilde/n)*quant.R.dot.star[c(4,3)]
    	logic.2 = ifelse(ci.ep.2[1] <= theta.0 & theta.0 <= ci.ep.2[2], 1, 0) 
    	len.2 = abs(ci.ep.2[2] - ci.ep.2[1])
    
    	ci.ep.2s = rep(theta.hat - b0/sqrt(n), 2)  + c(-1,1)*sqrt(Sigma.tilde/n)*h.check
    	logic.2s = ifelse(ci.ep.2s[1] <= theta.0 & theta.0 <= ci.ep.2s[2], 1, 0) 
    	len.2s = abs(ci.ep.2s[2] -ci.ep.2s[1])

    	pb.ci.out.nz.index = c(logic.L, logic.U, logic.2, logic.2s, len.2, len.2s)  # vector of length 6
    	return(pb.ci.out.nz.index)
}

dblasso.ci = function(y, X, coef.index, alpha, A.Z)
{
    	A = lasso.proj(X, y, Z = A.Z, suppress.grouptesting = TRUE)
    	u1 = confint(A, level = 1-alpha)
    	return(as.vector(u1[coef.index, ]))
}

data.set.m = function(m, y.mat, X, CV.k, coef.index, B, alpha, beta.0, k.fac, A.Z)
{
		y = y.mat[ ,m]
		n = nrow(X)
		p = ncol(X)
		R = length(k.fac)
		p0 = length(beta.0[beta.0!=0])
	
		d.vec = numeric(p)
		d.vec[coef.index] <- 1

		C = lasso.fn(y, X, list(CV.k, k.fac))  
		beta.mat = C[[1]] # pxR matrix 
    	lam = C[[2]] # R vector of lambda values #
    	CV.lam = C[[3]] 
    
		theta.vec = as.vector(d.vec%*%beta.mat) # R-vector
		inz.theta = (1:R)[theta.vec!=0] # indices with non zero theta.vec
		theta.0 = sum(d.vec*beta.0) # true parameter 
	
		if(length(inz.theta)>0)
		print(c(m, 1))
		else
		{
			print(c(m, 2))
		}
    	#if(m%%5 == 0){cat(m,"\n")}

	
		# the debiased Lasso based CI's are always usable # - Debiased Lasso CI - #
		dbl.ci = dblasso.ci(y, X, coef.index, alpha, A.Z)
        dbl.logic = ifelse(dbl.ci[1]<= theta.0 & theta.0 <= dbl.ci[2], 1, 0)
		dbl.len = abs(dbl.ci[2] - dbl.ci[1])
		dlasso.out.m <- c(dbl.logic, dbl.len) 
		
		if(length(inz.theta)>0)
		{
			# Creating bootstrap indices #
			mu.star = (1/2)/(1/2 + 3/2)
			G.mat = sapply(1:B, function(i, n.size){rbeta(n.size, 1/2, 3/2)}, n.size = n) # n x B matrix 
			pboot.imat = apply(G.mat, 2, function(v2, c1){(1/c1)*(v2-c1)}, c1 = mu.star) # n x B transformed from G.mat#
    
			init.info.list = list(beta.mat, lam, theta.vec, theta.0, pboot.imat)
			PB.out.mat.m <- matrix(0, nrow = R, ncol = 6)
		
				
			oracle.out.m <- matrix(0, nrow = R, ncol = 2)
			C11 = t(X[ ,1:p0])%*%(X[ ,1:p0])/n
			C11.inv = solve(C11)
			d1.vec = d.vec[1:p0] # p0 dimensional vector #
			rho.sqr = as.numeric(matrix(d1.vec, nrow = 1)%*%C11.inv%*%matrix(d1.vec, ncol = 1))
			sigma.check = numeric(R)
		
			for(r in inz.theta) # for these r, theta.vec[r] is non-zero 
			{
				#cat("\n")
			    sigma.check[r] = sqrt(mean((y - X%*%beta.mat[ ,r])^2))
			    oracle.ci.r = rep(theta.vec[r], 2) - qnorm(c(1-alpha/2, alpha/2))*sigma.check[r]*sqrt(rho.sqr)/sqrt(n)
		    
			    oracle.out.m[r, 1] = ifelse(oracle.ci.r[1] <= theta.0 & theta.0 <= oracle.ci.r[2], 1, 0)
			   	oracle.out.m[r, 2] = abs(oracle.ci.r[2] - oracle.ci.r[1])	  
			}
    	
    		# PB #
    		PB.out.mat.m[inz.theta, ] <- t(sapply(inz.theta, pbci.fn, y, X, coef.index, B, alpha, beta.0, init.info.list)) # length(inz.theta) x 6 matrix 
    		
    		return(list(c(1, m), inz.theta, oracle.out.m, PB.out.mat.m, dlasso.out.m))	
		}
		else # dlasso.out.m is the only output along with inz.theta, which should be a NULL 
		{
			return(list(c(2, m), inz.theta, dlasso.out.m))
		}	
}

lasso.all <- function(n, p, CV.k, coef.index, alpha, B.boot, err.type, k.fac)
{
		my.file.name = paste("yx-all-n-",n,"-p-",p,"-true-beta-5.Rdata", sep = "")
		load(my.file.name)	
		y.mat = yx.list[[6]][err.type, , ] # loads n x M response matrix for het-case 
		X.fix = yx.list[[2]]
		beta.0 = yx.list[[3]]
		A.Z = yx.list[[5]]
   		M = ncol(y.mat)
		R = length(k.fac)
		
		B <- mclapply(1:M, data.set.m, mc.cores = 12, y.mat, X.fix, CV.k, coef.index, B.boot, alpha, beta.0, k.fac, A.Z)
		
		save(B, file = paste("long-n-",n,"-p-",p,"-coef-",coef.index,"-et2.Rdata", sep = ""))
		# evaluating coverage probabilities #
		prelim.count <- matrix(0, nrow = M, ncol = 2)
		count.vec <- numeric(R)
		oracle.out.M <- matrix(0, nrow = R, ncol = 2)
		PB.out.mat.M <- matrix(0, nrow = R, ncol = 6)
		#colnames(PB.out.mat.M) <- c("ec.L", "ec.U", "ec.pb2", "ec.pb2s", "alen.pb2", "alen.pb2s")
		dlasso.out.M <- rep(0, 2)
		for(m in 1:M)
		{
			prelim.count[m, ] <- B[[m]][[1]]
			if(prelim.count[m, 1] == 1)
			{
				a.m <- B[[m]][[2]]
				count.vec[a.m] <- count.vec[a.m] + 1
				oracle.out.M <- oracle.out.M + B[[m]][[3]]
				PB.out.mat.M <- PB.out.mat.M + B[[m]][[4]]
				dlasso.out.M <- dlasso.out.M + B[[m]][[5]]
			}
			else
			{
				dlasso.out.M <- dlasso.out.M + B[[m]][[3]]
			}
		}			
		for(r in 1:R)
		{
			oracle.out.M[r, ] <- oracle.out.M[r, ]/count.vec[r]
			PB.out.mat.M[r, ] <- PB.out.mat.M[r, ]/count.vec[r]
		}
		dlasso.out.M <- dlasso.out.M/M 
		
		#print(signif(cbind(k.fac, count.vec, oracle.out.M, PB.out.mat.M),3))
		#print(dlasso.out.M)
		C1 <- list(prelim.count, k.fac, count.vec, oracle.out.M, PB.out.mat.M, dlasso.out.M)
		save(C1, file = paste("short-n-",n,"-p-",p,"-coef-",coef.index,"-et2.Rdata", sep = ""))
		
		
}
k0.seq= seq(5,0.375,by=-.125)
lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 700, err.type = 2, k.fac = k0.seq)

#cat("XXXX\nXXXX\n")
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 700, err.type = 2, k.fac = k0.seq)

