rm(list = ls())
# plotting code for heteroskedastic case with PB
plot.code <- function(n, p, coef.index)
{
    file.name.1 = paste("short-n-",n,"-p-",p,"-coef-",coef.index,"-et2.Rdata", sep = "")
    load(file.name.1)

    k.fac <- C1[[2]][-(1:8)]
    count.vec <- C1[[3]][-(1:8)]
    oracle.out <- C1[[4]][-(1:8), ]
    PB.out.mat <- C1[[5]][-(1:8), ]
    # removes the first 8 rows, as they correspond to k.fac \in (4,5], too high penalty
    
    # for our benefit we store the column names in oracle.out and PB.out.mat#
    # ec.orac and alen.orac
    # "ec.L", "ec.U", "ec.pb2", "ec.pb2s", "alen.pb2", "alen.pb2s" 
    # --------------------------------------------------------------------------- #
    
    C1.new.2sided <- cbind(oracle.out[ ,1], PB.out.mat[ ,3:4]) # ec.orac, ec.pb2, ec.pb2s
    u.2sided <- c(min(C1.new.2sided), max(C1.new.2sided))
    
    C1.new.2sided.len <- cbind(oracle.out[ ,2], PB.out.mat[ ,5:6]) # alen.orac, alen.pb2, alen.pb2s
    u.2sided.len <- c(min(C1.new.2sided.len), max(C1.new.2sided.len))
    
    C1.new.1sided <- cbind(PB.out.mat[ ,1:2])
    u.1sided <- c(min(C1.new.1sided), max(C1.new.1sided))
    
    # empirical coverage and avg. length of two sided CI #
    my.file.name.new.2 = paste("plots-all-heteroskedastic-n-",n,"-p-",p,"-coef-beta-",coef.index,"-et2.pdf", sep = "")
    pdf(my.file.name.new.2, width = 10, height = 15)

    par(mfrow = c(3,1))

    plot(k.fac, C1.new.2sided[ ,1], type = "l", xlab = "k", ylab = "empirical coverage", col = 1, lwd = 1.25, ylim = u.2sided)
    lines(k.fac, C1.new.2sided[ ,2], col = 4, lwd = 1.25)
    title(main = "empirical coverage of 90% two-sided CIs: heteroskedastic case", cex.main = 1)
    lines(k.fac, C1.new.2sided[ ,3], col = "yellow4", lwd = 1.25)
    abline(h = 0.9, col = "red", lty = 3)
    legend(x="bottomleft", lty = c(1,1,1), col = c("black","blue","yellow4"), legend = c("Orac.2", "PB.2", "PB.2Sym"), bty = "n", cex = 1, lwd = rep(1.25, 3))
    
 
    
    plot(k.fac, C1.new.2sided.len[ ,1], type = "l", xlab = "k", ylab = "average length", col = 1, lwd = 1.25, ylim = u.2sided.len)
    lines(k.fac, C1.new.2sided.len[ ,2], col = 4, lwd = 1.25)
    title(main = "average length of 90% two-sided CIs: heteroskedastic case", cex.main = 1)
    lines(k.fac, C1.new.2sided.len[ ,3], col = "yellow4", lwd = 1.25)
    legend(x="topleft", lty = c(1,1,1), col = c("black","blue","yellow4"), legend = c("Orac.2", "PB.2", "PB.2Sym"), bty = "n", cex = 1, lwd = rep(1.25, 3))
    
    
    # empirical coverage of one sided CI #
    
    plot(k.fac, C1.new.1sided[ ,1], type = "l", col = "blue", xlab = "k", ylab = "empirical coverage", lwd = 1.25, ylim = u.1sided)
    lines(k.fac, C1.new.1sided[ ,2], col = "red4", lwd = 1.25)
    title(main = "empirical coverage of 90% one-sided lower and upper PB based CIs': heteroskedastic case", cex.main = 1)
    abline(h = 0.9, col = 1, lty = 3)
    legend(x="topleft", lty = c(1,1), col = c("blue","red4"), legend = c("PB.Lower", "PB.Upper"), bty = "n", cex = 1, lwd = rep(1.25, 2))
   
     dev.off()
}
plot.code(n = 300, p = 500, coef.index = 10)
