rm(list = ls())
# plotting code for homoskedastic case
plot.code <- function(n, p, coef.index, err.type)
{
    # plotting for the homoskedastic case # 
    
    file.name.1 = paste("short-n-",n,"-p-",p,"-coef-",coef.index,"-et-",err.type,".Rdata", sep = "")
    load(file.name.1)

    
    C1.new <- C1[[8]][-(1:16), ]  # removes the first 16 rows from C1, as they correspond to k.fac \in (3,5], too high penalty
    
    # for our benefit we store the column names # 
    # > colnames(C1.new)
    # [1] "k.fac"     "count.vec" "or.ec"     "or.alen"   "rbL"       "rbU"      
    # [7] "rb2"       "rb2s"      "rb2.alen"  "rb2s.alen" "pbL"       "pbU"      
    #[13] "pb2"       "pb2s"      "pb2.alen"  "pb2s.alen"
    # --------------------------------------------------------------------------- #
    
    k.fac <- C1.new[ ,1]
    count.vec <- C1.new[ ,2]
    C1.new.2sided <- C1.new[ , c(3, 7, 8, 13, 14)] # or.ec, rb2, rb2s, pb2, pb2s
    C1.new.2sided.len <- C1.new[ , c(4, 9, 10, 15, 16)] # or.alen, rb2.alen, rb2s.alen, pb2.alen, pb2s.alen
    
    # empirical coverage and avg. length of two sided CI #
    my.file.name.new.2 = paste("plots-all-homoskedastic-n-",n,"-p-",p,"-coef-beta-",coef.index,"-errtype-",err.type,".pdf", sep = "")
    pdf(my.file.name.new.2, width = 14, height = 12)

    par(mfrow = c(2,2))
    min.2s <- min(C1.new.2sided)
    max.2s <- max(C1.new.2sided)
    
    plot(k.fac, C1.new[ ,3], type = "l", xlab = "k", ylab = "empirical coverage", col = 1, lwd = 1.25, ylim = c(min.2s, max.2s))
    lines(k.fac, C1.new[ ,7], col = 4, lwd = 1.25)
    title(main = "empirical coverage of 90% two-sided CIs: homoskedastic case", cex.main = 1)
    lines(k.fac, C1.new[ ,8], col = "yellow4", lwd = 1.25)
    lines(k.fac, C1.new[ ,13], col = "red4", lty = 2, lwd = 1.25)
    lines(k.fac, C1.new[ ,14], col = "seagreen4", lty = 2, lwd = 1.25)
    #lines(k.fac, rep(0.9, length(k.fac)), lty = 3, lwd = 1.25, col = "red")
    abline(h = 0.9, col = "red", lty = 3)
    legend(x="topleft", lty = c(1,1,1,2,2), col = c("black","blue","yellow4","red4", "seagreen4"), legend = c("Orac.2", "RB.2", "RB.2Sym", "PB.2", "PB.2Sym"), bty = "n", cex = 1, lwd = rep(1.25, 5))
    
 
    min.2slen <- min(C1.new.2sided.len)
    max.2slen <- max(C1.new.2sided.len)
    
           
    plot(k.fac, C1.new[ ,4], type = "l", xlab = "k", ylab = "average length", col = 1, lwd = 1.25, ylim = c(min.2slen, max.2slen))
    lines(k.fac, C1.new[ ,9], col = 4, lwd = 1.25)
    title(main = "average length of 90% two-sided CIs: homoskedastic case", cex.main = 1)
    lines(k.fac, C1.new[ ,10], col = "yellow4", lwd = 1.25)
    lines(k.fac, C1.new[ ,15], col = "red4", lty = 2, lwd = 1.25)
    lines(k.fac, C1.new[ ,16], col = "seagreen4", lty = 2, lwd = 1.25)
    #lines(k.fac, rep(0.9, length(k.fac)), lty = 3, lwd = 1.25, col = "red")
    #abline(h = 0.9, col = "red", lty = 3)
    legend(x="topleft", lty = c(1,1,1,2,2), col = c("black","blue","yellow4","red4", "seagreen4"), legend = c("Orac.2", "RB.2", "RB.2Sym", "PB.2", "PB.2Sym"), bty = "n", cex = 1, lwd = rep(1.25, 5))
    
    #dev.off()
    
    # empirical coverage of one sided CI #
    
    C1.new.1L <- C1.new[ , c(5, 11)] # rbL, pbL
    C1.new.1U <- C1.new[ , c(6, 12)] # rbU, pbU

    #my.file.name.new.1 = paste("plot1s-n-",n,"-p-",p,"-coef-beta-",coef.index,"-errtype-",err.type,".pdf", sep = "")
    #pdf(my.file.name.new.1, width = 10, height = 10)

    #par(mfrow = c(2,1))
    min.1L <- min(C1.new.1L)
    max.1L <- max(C1.new.1L)

    min.1U <- min(C1.new.1U)
    max.1U <- max(C1.new.1U)
    
    plot(k.fac, C1.new[ ,5], type = "l", col = "blue", xlab = "k", ylab = "empirical coverage", lwd = 1.25, ylim = c(min.1L, max(0.9,max.1L)))
    lines(k.fac, C1.new[ ,11], col = "red4", lwd = 1.25)
    title(main = "empirical coverage of 90% one-sided lower CIs: homoskedastic case", cex.main = 1)
    abline(h = 0.9, col = 1, lty = 3)
    legend(x="topleft", lty = c(1,1), col = c("blue","red4"), legend = c("RB.Lower", "PB.Lower"), bty = "n", cex = 1, lwd = rep(1.25, 2))
   
    plot(k.fac, C1.new[ ,6], type = "l", col = "blue", xlab = "k", ylab = "empirical coverage", lwd = 1.25, ylim = c(min.1U, max.1U))
    lines(k.fac, C1.new[ ,12], col = "red4", lwd = 1.25)
    title(main = "empirical coverage of 90% one-sided upper CIs: homoskedastic case", cex.main = 1)
    abline(h = 0.9, col = 1, lty = 3)
    legend(x="topleft", lty = c(1,1), col = c("blue","red4"), legend = c("RB.Upper", "PB.Upper"), bty = "n", cex = 1, lwd = rep(1.25, 2))

   dev.off()
}
plot.code(n = 150, p = 500, coef.index = 10, err.type = 2)
