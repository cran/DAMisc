mnlSig <- function(obj, pval=.05, two.sided=TRUE, flag.sig=TRUE, insig.blank=FALSE){
	smulti <- summary(obj)
	multi.t <- smulti$coefficients/smulti$standard.errors
	multi.p <- (2^as.numeric(two.sided))*pnorm(abs(multi.t), lower.tail=F)
	b <- matrix(sprintf("%.3f", smulti$coefficients), ncol=ncol(multi.t))
	sig.vec <- c(" ", "*")
	sig.obs <- as.numeric(multi.p < pval) + 1
	if(flag.sig){
		b <- matrix(paste(b, sig.vec[sig.obs], sep=""), ncol=ncol(multi.t))
	}
	if(insig.blank){
	b[which(multi.p > pval, arr.ind=T)] <- ""
	}
	rownames(b) <- rownames(multi.t)
	colnames(b) <- colnames(multi.t)
	b <- as.data.frame(b)
	return(b)
}
