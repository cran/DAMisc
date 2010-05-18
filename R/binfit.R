binfit <- function(mod){
	mod <- update(mod, x=TRUE, y=TRUE)
	y <- mod$y
	null.mod <- update(mod, ".~1")
	b <- mod$coef[-1]
	var.ystar <- t(b) %*% var(model.matrix(mod)[,-1]) %*% b
	G <- -2*(logLik(null.mod) - logLik(mod))
	res.col1 <- c(logLik(null.mod), deviance(mod), NA, 
		1-(logLik(mod)/logLik(null.mod)), 
		1-(exp(logLik(null.mod))/exp(logLik(mod)))^(2/length(mod$residuals)), 
		var.ystar / (var.ystar + switch(mod$family[[2]], logit = pi^2/3, probit=1)), 
		mean(mod$y == as.numeric(fitted(mod) > .5)), stats:::BIC(mod))
	res.col2 <- c(logLik(mod), G, pchisq(G, 8, lower.tail=F), 
		1-((logLik(mod)-mod$rank)/logLik(null.mod)), 
		res.col1[5]/(1-(exp(logLik(null.mod))^(2/length(mod[["residuals"]])))), 
		1-(sum((y-fitted(mod))^2)/sum((y-mean(y))^2)), 
		(sum(mod$y == as.numeric(fitted(mod) > .5)) - 
			max(table(y)))/(length(mod$residuals) - max(table(y))), 
		stats:::AIC(mod))
	res.vec <- c(res.col1, res.col2)[-3]
	res.col1 <- sprintf("%3.3f", res.col1)
	res.col1[3] <- ""
	res.df <- data.frame(
	Names1 = c("Log-Lik Intercept Only:", paste("D(", mod$df.residual, "):", sep=""), 
		" ", "McFadden's R2:", "ML (Cox-Snell) R2:", "McKelvey & Zavoina R2:", 
		"Count R2:", "BIC:"), 
	vals1 = res.col1, 
	Names2 = c("Log-Lik Full Model:", paste("LR(", mod$rank-1, "):", sep=""), 
		"Prob > LR:", "McFadden's Adk R2:", "Cragg-Uhler (Nagelkerke) R2:", 
		"Efron's R2:", "Adj Count R2:", "AIC:"), 
	vals2 = sprintf("%3.3f", res.col2))
	return(res.df)
}
