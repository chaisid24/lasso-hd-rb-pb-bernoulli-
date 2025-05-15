rm(list = ls())
library(parallel)
library(glmnet)
library(hdi)

lasso.fn = function(y, X.scale, CV.list) # use this code when cv.lam.search is required 
{
	n = nrow(X.scale)
	p = ncol(X.scale)
	CV.k = CV.list[[1]] # the k-fold parameter #
    k.fac = CV.list[[2]] # vector of length R # sequence should be in decreasing order
    
    CV.lam = cv.glmnet(X.scale, y, nfolds = CV.k, intercept = F, standardize = F)$lambda.min
    lam1 = k.fac*CV.lam # R length # 
    A1 = glmnet(X.scale, y, intercept = F, standardize = F, lambda = lam1)
    beta.lasso.1 = A1$beta # a pxR matrix
    
  	out.all = list(beta.lasso.1, lam1)
	#out.all = list(beta.lasso.1, lam1, beta.lasso.2, lam2) 
	return(out.all)
}

lasso.within.boot <- function(b, y.boot.mat, X, lam.fix)
{
    y.star = y.boot.mat[ ,b]
    A = glmnet(X, y.star, intercept = F, standardize = F, lambda = lam.fix)
	beta.star = as.vector(A$beta)
	return(beta.star)
}

rb.mid.step.vareq = function(b, y.star.mat, X, beta.star.mat, lam, d.vec, beta.plols)
{
    y.star = y.star.mat[ ,b]
    beta.star = beta.star.mat[ ,b]
        
    n = nrow(X)
    p = ncol(X)
    T.star = sqrt(n)*sum(d.vec*(beta.star - beta.plols))
    
    beta.dot.star = rep(0, p)
    
    A.star = (1:p)[beta.star!=0]
    C11.star = t(X[ ,A.star])%*%(X[ ,A.star])/n
    c.star = sign(beta.star[A.star])
    term.1 = as.vector(matrix(d.vec[A.star], nrow = 1)%*%solve(C11.star)%*%matrix(c.star, ncol = 1))
    b0.star = -sqrt(n)*lam*term.1
    
    beta.dot.star[A.star] = beta.star[A.star] + lam*as.vector(solve(C11.star)%*%matrix(c.star, ncol = 1))
    sigma.star = sqrt(mean((y.star - X%*%beta.dot.star)^2))
    R.check.star = (T.star - b0.star)/sqrt(sigma.star)
    return(R.check.star)
}

rbci.vareq.fn = function(nz.index, y, X, coef.index, B, alpha, beta.0, init.info.list)
{
    #cat(nz.index,"\n")
	n = nrow(X)
	p = ncol(X)
	
	d.vec = numeric(p)
	d.vec[coef.index] = 1 

	beta.hat = init.info.list[[1]][ ,nz.index]   # that column of beta.mat
	lam = init.info.list[[2]][nz.index]          # that element of lambda.vec 
	theta.hat = init.info.list[[3]][nz.index]    # that element of theta.vec
	theta.0 = init.info.list[[4]]
    rboot.imat = init.info.list[[5]]             # nxB matrix of sampled indices
    
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
    
    rb.ci.out.nz.index = c(logic.L, logic.U, logic.2, logic.2s, len.2, len.2s)  # vector of length 6
    return(rb.ci.out.nz.index)
}

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

pbci.vareq.fn = function(nz.index, y, X, coef.index, B, alpha, beta.0, init.info.list)
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
    pboot.imat = init.info.list[[6]]             # nxB matrix 
   
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
    w4 = 2*mean((xi0.vec^4)*(res.check^4))/(Sigma.tilde^2) + (4/(sigma.check.sqr*Sigma.tilde))*mean((xi0.vec^2)*(res.check^4)) 
        - (3/(sigma.check.sqr^2))*mean(res.check^4) + 1
    Cn.p = -(z.alpha/n)*(w2/2 + w4*(z.alpha^2 - 3)/(24))
    h.dag = quantile(abs(Rt.ss.vec), probs = 1-alpha)
    h.check = h.dag + Cn.p
    
    # CI construction #
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
	theta.vec = as.vector(d.vec%*%beta.mat) # R-vector
	inz.theta = (1:R)[theta.vec!=0] # indices with non zero theta.vec
	theta.0 = sum(d.vec*beta.0) # true parameter 
	
	
	if(length(inz.theta)>0)
	{
		print(c(m, 1))
		#print(cbind(inz.theta, theta.vec[inz.theta]))
		
	}
	else
	{
		print(c(m, 2))
	}
	
	# the debiased Lasso based CI's are always usable # - Debiased Lasso CI - #
	dbl.ci = dblasso.ci(y, X, coef.index, alpha, A.Z)
    dbl.logic = ifelse(dbl.ci[1]<= theta.0 & theta.0 <= dbl.ci[2], 1, 0)
	dbl.len = abs(dbl.ci[2] - dbl.ci[1])
	dlasso.out.m <- c(dbl.logic, dbl.len) 
	
	
	if(length(inz.theta)>0)
	{
		# Creating bootstrap indices #
		rboot.imat = sapply(1:B, function(i, z){sample(z, size = length(z), replace = T)}, z = 1:n)
		mu.star = (1/2)/(1/2 + 3/2)
		G.mat = sapply(1:B, function(i, n.size){rbeta(n.size, 1/2, 3/2)}, n.size = n) # n x B matrix 
		pboot.imat = apply(G.mat, 2, function(v2, c1){(1/c1)*(v2-c1)}, c1 = mu.star) # n x B transformed from G.mat#
    
		init.info.list = list(beta.mat, lam, theta.vec, theta.0, rboot.imat, pboot.imat)
	
		# oracle ci #
    	oracle.out.m <- matrix(0, nrow = R, ncol = 2)
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
    	RB.out.mat.m <- PB.out.mat.m <- matrix(0, nrow = R, ncol = 6)
    	RB.out.mat.m[inz.theta, ] = t(sapply(inz.theta, rbci.vareq.fn, y, X, coef.index, B, alpha, beta.0, init.info.list)) # length(inz.theta)x6 matrix
    	PB.out.mat.m[inz.theta, ] = t(sapply(inz.theta, pbci.vareq.fn, y, X, coef.index, B, alpha, beta.0, init.info.list)) 
    	
    	return(list(c(1,m), inz.theta, oracle.out.m, RB.out.mat.m, PB.out.mat.m, dlasso.out.m))
    }
    else
    {
    	return(list(c(2,m), inz.theta, dlasso.out.m))
    }
  
}

lasso.all <- function(n, p, CV.k, coef.index, alpha, B.boot, err.type, k.fac)
{
		my.file.name = paste("yx-all-n-",n,"-p-",p,"-true-beta-5.Rdata", sep = "")
		load(my.file.name)
		
		y.mat = yx.list[[1]][err.type, , ] # loads n x M response matrix for homogenous-case 
		X.fix = yx.list[[2]]
		beta.0 = yx.list[[3]]
		A.Z = yx.list[[5]]
   		M = ncol(y.mat)
		R = length(k.fac)
		
		B <- mclapply(1:M, data.set.m, mc.cores = 8, y.mat, X.fix, CV.k, coef.index, B.boot, alpha, beta.0, k.fac, A.Z)
		
		save(B, file = paste("long-n-",n,"-p-",p,"-coef-",coef.index,"-et-",err.type,".Rdata", sep = ""))
		# evaluating coverage probabilities #
		
		prelim.count <- matrix(0, nrow = M, ncol = 2)
		count.vec <- numeric(R)
		oracle.out.M <- matrix(0, nrow = R, ncol = 2)
		RB.out.mat.M <- PB.out.mat.M <- matrix(0, nrow = R, ncol = 6)
		dlasso.out.M <- rep(0, 2)
		for(m in 1:M)
		{
			prelim.count[m, ] <- B[[m]][[1]]
			if(prelim.count[m, 1] == 1)
			{
				a.m <- B[[m]][[2]]
				count.vec[a.m] <- count.vec[a.m] + 1
				oracle.out.M <- oracle.out.M + B[[m]][[3]]
				RB.out.mat.M <- RB.out.mat.M + B[[m]][[4]]
				PB.out.mat.M <- PB.out.mat.M + B[[m]][[5]]
				dlasso.out.M <- dlasso.out.M + B[[m]][[6]]
			}
			else
			{
				dlasso.out.M <- dlasso.out.M + B[[m]][[3]]
			}
		}			
		for(r in 1:R)
		{
			if(count.vec[r]>0)
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
		dlasso.out.M <- matrix(dlasso.out.M/M, nrow = 1, ncol = 2)
		
		result.mat <- cbind(k.fac, count.vec, oracle.out.M, RB.out.mat.M, PB.out.mat.M)
		colnames(result.mat) <- c("k.fac", "count.vec", "or.ec", "or.alen", "rbL", "rbU", "rb2", "rb2s", "rb2.alen", "rb2s.alen", "pbL", "pbU", "pb2", "pb2s", "pb2.alen", "pb2s.alen")
		print(signif(result.mat, 4))
		colnames(dlasso.out.M) <- c("dlasso.ec", "dlasso.alen")
		print(dlasso.out.M)
		C1 <- list(prelim.count, k.fac, count.vec, oracle.out.M, RB.out.mat.M, PB.out.mat.M, dlasso.out.M, result.mat)
		save(C1, file = paste("short-n-",n,"-p-",p,"-coef-",coef.index,"-et-",err.type,".Rdata", sep = ""))
		
}
k0.seq= seq(5,0.375,by=-.125)


#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 700, err.type = 1, k.fac = k0.seq)
#cat("####\n####\n")
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 700, err.type = 1, k.fac = k0.seq)
#cat("####\n####\n")
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 700, err.type = 2, k.fac = k0.seq)
#cat("####\n####\n")
#lasso.all(n = 300, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 700, err.type = 2, k.fac = k0.seq)
#cat("####\n####\n")
########################
# RUNNING n = 150 case #

#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 700, err.type = 1, k.fac = k0.seq)
#cat("####\n####\n####\n####\n")
#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 700, err.type = 1, k.fac = k0.seq)
#cat("####\n####\n####\n####\n")
#lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 1, alpha = 0.1, B.boot = 700, err.type = 2, k.fac = k0.seq)
#cat("####\n####\n####\n####\n")
lasso.all(n = 150, p = 500, CV.k = 5, coef.index = 10, alpha = 0.1, B.boot = 700, err.type = 2, k.fac = k0.seq)
#cat("####\n####\n####\n####\n")
