simPredpolr <- function (object, coefs, n.coef) 
{
X <- model.matrix(object)
xint <- match("(Intercept)", colnames(X), nomatch = 0L)
if (xint > 0L) 
    X <- X[, -xint, drop = FALSE]
n <- nrow(X)
q <- length(coefs) - n.coef
eta <- {if(n.coef > 0){
		drop(X %*% coefs[1:n.coef])
	} else{
		rep(0, n)
	}
}
pfun <- switch(object$method, logistic = plogis, probit = pnorm, 
    cauchit = pcauchy)
cumpr <- matrix(pfun(matrix(coefs[(n.coef+1):length(coefs)], n, q, byrow = TRUE) - 
    eta), , q)
Y <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
return(Y)
}

