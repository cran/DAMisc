poisGOF <- function(obj){
	if(!("glm" %in% class(obj))){
		stop("poisGOF only works on objects of class glm (with a poisson family)\n")
	} 
	if(!("y" %in% names(obj))){
		obj <- update(obj, y=TRUE)
	}
	ind.chisq <- (((obj$y - obj$fitted)^2)/obj$fitted)
	dev <- obj$deviance
	df <- obj$df.residual
	p.chisq <- pchisq(sum(ind.chisq), df=df, lower.tail=FALSE)
	p.dev <- pchisq(dev, df=df, lower.tail=FALSE)
	vec <- sprintf("%.3f", c(sum(ind.chisq), dev, p.chisq, p.dev))
	mat <- matrix(vec, ncol=2, nrow=2)
	rownames(mat) <- c("Chi-squared", "Deviance")
	colnames(mat) <- c("Stat", "p-value")
	mat <- as.data.frame(mat)
	return(mat)
}
