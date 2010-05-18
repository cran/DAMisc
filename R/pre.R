pre <- function(mod1, mod2 = NULL, sim = FALSE, R = 2500){
require(MASS)
if(!is.null(mod2)){
 	if(mean(class(mod1) == class(mod2)) != 1){
	stop("Model 2 must be either NULL or of the same class as Model 1\n")
	}
}
if(!any(class(mod1) %in% c("polr", "multinom", "glm"))){
	stop("pre only works on models of class glm (with binomial family), polr or multinom\n")
}
if("glm" %in% class(mod1)){
	if(!(family(mod1)$link %in% c("logit", "probit", "cloglog", "cauchit"))){
		stop("PRE only calculated for models with logit, probit, cloglog or cauchit links\n")
	}
  if(is.null(mod2)){
  	y <- mod1[['y']]
  	mod2 <- update(mod1, ". ~ 1")
  }
pred.mod2 <- as.numeric(
  predict(mod2, type="response") >= .5)
pmc <- mean(mod2$y == pred.mod2)
pred.y <- as.numeric(
    predict(mod1, type="response") >= .5)
pcp <- mean(pred.y == mod1$y)
pre <- (pcp-pmc)/(1-pmc)

pred.prob1 <- predict(mod1, type="response")
pred.prob2 <- predict(mod2, type="response")
epcp <- (1/length(pred.prob1))*(
    sum(pred.prob1[which(mod1$y == 1)]) + 
    sum(1-pred.prob1[which(mod1$y == 0)]))
epmc <- (1/length(pred.prob2))*(
    sum(pred.prob2[which(mod2$y == 1)]) + 
    sum(1-pred.prob2[which(mod2$y == 0)]))
epre <- (epcp - epmc)/(1-epmc)
if(sim){
b1.sim <- mvrnorm(R, 
    coef(mod1), vcov(mod1))
b2.sim <- mvrnorm(R, 
    coef(mod2), vcov(mod2))
mod1.probs <- family(mod1)$linkinv(
  model.matrix(mod1) %*% t(b1.sim))
mod2.probs <- family(mod2)$linkinv(
  model.matrix(mod2) %*% t(b2.sim))
pmcs <- apply(mod2.probs, 2, 
  function(x)mean(as.numeric(x > .5) == mod2$y))
pcps <- apply(mod1.probs, 2, function(x)mean(
    as.numeric(x>.5) == mod1$y))
pre.sim <- (pcps-pmcs)/(1-pmcs)

epmc.sim <- apply(mod2.probs, 2, function(x)
    (1/length(x))*
    (sum(x[which(mod2$y == 1)]) + 
    sum(1-x[which(mod2$y == 0)])))

epcp.sim <- apply(mod1.probs, 2, 
  function(x)(1/length(x))*
  (sum(x[which(mod1$y == 1)]) + 
  sum(1-x[which(mod1$y == 0)])))

epre.sim <- (epcp.sim - epmc.sim)/(1-epmc.sim)
}  
}
if("multinom" %in% class(mod1)){
    mod1 <- update(mod1, ".~.", model=TRUE, trace=FALSE)
	if(is.null(mod2)){
	    mod2 <- update(mod1, ". ~ 1", data=mod1$model, trace=FALSE)
	  }

pred.prob1 <- predict(mod1, type="prob")
pred.prob2 <- predict(mod2, type="prob")
pred.cat1 <- apply(pred.prob1, 1, which.max)
pred.cat2 <- apply(pred.prob2, 1, which.max)
pcp <- mean(as.numeric(mod1$model[,1]) == pred.cat1)
pmc <- max(table(as.numeric(mod2$model[,1]))/sum(table(mod2$model[,1])))
pre <- (pcp-pmc)/(1-pmc)

epcp <- mean(pred.prob1[cbind(1:nrow(pred.prob1), 
    as.numeric(mod1$model[,1]))])
tab <- table(mod1$model[,1])/
    sum(table(mod1$model[,1]))
epmc <- mean(tab[as.numeric(mod2$model[,1])])
epre <- (epcp - epmc)/(1-epmc)
if(sim){
	b1.sim <- mvrnorm(R, 
	    c(t(coef(mod1))), vcov(mod1))
	b2.sim <- mvrnorm(R, 
	    c(t(coef(mod2))), vcov(mod2))
	mod.levs <- rownames(coef(mod1))
	var.levs <- rownames(contrasts(mod1$model[,1]))
	tmp.reord <- c(match(mod.levs, var.levs), (1:length(var.levs))[-match(mod.levs, var.levs)])
	reord <- match(1:length(var.levs), tmp.reord)
	tmp <- lapply(1:nrow(b1.sim), function(x)matrix(b1.sim[x,], ncol=ncol(coef(mod1)), byrow=T))
	tmp2 <- lapply(tmp, function(x)cbind(model.matrix(mod1)%*%t(x), 0))
	tmp2 <- lapply(tmp2, function(x)x[, reord])
	mod1.probs <- lapply(tmp2, function(x)exp(x)/apply(exp(x), 1, sum))
	tmp <- lapply(1:nrow(b2.sim), function(x)matrix(b2.sim[x,], ncol=ncol(coef(mod2)), byrow=T))
	tmp2 <- lapply(tmp, function(x)cbind(model.matrix(mod2)%*%t(x), 0))
	tmp2 <- lapply(tmp2, function(x)x[, reord])
	mod2.probs <- lapply(tmp2, function(x)exp(x)/apply(exp(x), 1, sum))

	pred.cat1 <- lapply(mod1.probs, function(x)apply(x, 1, which.max))
	pred.cat2 <- lapply(mod2.probs, function(x)apply(x, 1, which.max))
	pcp.sim <- sapply(pred.cat1, function(x)mean(x == as.numeric(mod1$model[,1])))
	pmc.sim <- sapply(pred.cat2, function(x)mean(x == as.numeric(mod1$model[,1])))
	pre.sim <- (pcp.sim-pmc.sim)/(1-pmc.sim)
	epcp.sim <- sapply(mod1.probs, function(z)mean(z[cbind(1:nrow(mod1$model), as.numeric(mod1$model[,1]))]))
	epmc.sim <- sapply(mod2.probs, function(z)mean(z[cbind(1:nrow(mod1$model), as.numeric(mod1$model[,1]))]))
	epre.sim <- (epcp.sim - epmc.sim)/(1-epmc.sim)
}
}

if("polr" %in% class(mod1)){
if(is.null(mod1$Hess)){
	mod1 <- update(mod1, Hess=TRUE)
}	
if(is.null(mod2)){
    mod2 <- update(mod1, ". ~ 1", data=mod1$model, model=TRUE, Hess=TRUE)
  }
pred.prob1 <- predict(mod1, type="prob")
pred.prob2 <- predict(mod2, type="prob")
pred.cat1 <- apply(pred.prob1, 1, which.max)
pred.cat2 <- apply(pred.prob2, 1, which.max)
pcp <- mean(as.numeric(mod1$model[,1]) == pred.cat1)
pmc <- max(table(as.numeric(mod2$model[,1]))/sum(table(mod2$model[,1])))
pre <- (pcp-pmc)/(1-pmc)

epcp <- mean(pred.prob1[cbind(1:nrow(pred.prob1), 
    as.numeric(mod1$model[,1]))])
tab <- table(mod1$model[,1])/
    sum(table(mod1$model[,1]))
epmc <- mean(tab[as.numeric(mod2$model[,1])])
epre <- (epcp - epmc)/(1-epmc)
if(sim){
	b1 <- c(coef(mod1), mod1$zeta)
	b2 <- c(coef(mod2), mod2$zeta)
	v1 <- invisible(vcov(mod1))
	v2 <- invisible(vcov(mod2))
	b1.sim <- mvrnorm(R, b1, v1)
	b2.sim <- mvrnorm(R, b2, v2)
	mod1.probs <- lapply(1:nrow(b1.sim), function(x)simPredpolr(mod1, b1.sim[x,], n.coef = length(coef(mod1))))
	mod2.probs <- lapply(1:nrow(b2.sim), function(x)simPredpolr(mod2, b2.sim[x,], n.coef = length(coef(mod2))))

	pred.cat1 <- lapply(mod1.probs, function(x)apply(x, 1, which.max))
	pred.cat2 <- lapply(mod2.probs, function(x)apply(x, 1, which.max))
	pcp.sim <- sapply(pred.cat1, function(x)mean(x == as.numeric(mod1$model[,1])))
	pmc.sim <- sapply(pred.cat2, function(x)mean(x == as.numeric(mod1$model[,1])))
	pre.sim <- (pcp.sim-pmc.sim)/(1-pmc.sim)
	epcp.sim <- sapply(mod1.probs, function(z)mean(z[cbind(1:nrow(mod1$model), as.numeric(mod1$model[,1]))]))
	epmc.sim <- sapply(mod2.probs, function(z)mean(z[cbind(1:nrow(mod1$model), as.numeric(mod1$model[,1]))]))
	epre.sim <- (epcp.sim - epmc.sim)/(1-epmc.sim)
}
}
ret <- list()
ret$pre <- pre
ret$epre <- epre
form1 <- formula(mod1)
form2 <- formula(mod2)
ret$m1form <- paste(form1[2], form1[1], form1[3], sep=" ")
ret$m2form <- paste(form2[2], form2[1], form2[3], sep=" ")
ret$pcp <- pcp
ret$pmc <- pmc
ret$epmc <- epmc
ret$epcp <- epcp
if(sim){
  ret$pre.sim <- pre.sim
  ret$epre.sim <- epre.sim
}
class(ret) <- "pre"
return(ret)
}  


print.pre <- function(x, ..., sim.ci = .95){
  cat("mod1: ", as.character(x$m1form), "\n")
  cat("mod2: ", as.character(x$m2form), "\n\n")
  cat("Analytical Results\n")
  cat(" PMC = ", sprintf("%2.3f", x$pmc), "\n")
  cat(" PCP = ", sprintf("%2.3f", x$pcp), "\n")
  cat(" PRE = ", sprintf("%2.3f", x$pre), "\n")
  cat("ePMC = ", sprintf("%2.3f", x$epmc), "\n")
  cat("ePCP = ", sprintf("%2.3f", x$epcp), "\n")
  cat("ePRE = ", sprintf("%2.3f", x$epre), "\n\n")
  low <- (1-sim.ci)/2
  up <- 1-low
  if(exists("pre.sim", x)){
    pre.ci <- sprintf("%2.3f", quantile(x$pre.sim, c(.5,low, up)))
    epre.ci <- sprintf("%2.3f", quantile(x$epre.sim, c(.5,low, up)))
    tmp <- rbind(pre.ci, epre.ci)
    rownames(tmp) <- c(" PRE", "ePRE")
    colnames(tmp) <- c("median", "lower", "upper")
    cat("Simulated Results\n")
    print(tmp, quote=F)
    }  
}