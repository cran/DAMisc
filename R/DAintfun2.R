DAintfun2 <- function(obj, varnames, name.stem = "cond_eff", plot.type = "screen"){
	rseq <- function(x){rx <- range(x, na.rm=TRUE); seq(rx[1], rx[2], length=25)}
	v1 <- varnames[1]
	v2 <- varnames[2]
	ind1 <- grep(v1, names(obj$coef))
	ind2 <- grep(v2, names(obj$coef))
	s1 <- rseq(model.matrix(obj)[,v1])
	s2 <- rseq(model.matrix(obj)[,v2])
	a1 <- a2 <- matrix(0, nrow=25, ncol=ncol(model.matrix(obj)))
	a1[,ind1[1]] <- 1
	a1[, ind1[2]] <- s2
	a2[,ind2[1]] <- 1
	a2[,ind2[2]] <- s1
	eff1 <- a1 %*% obj$coef
	se.eff1 <- sqrt(diag(a1 %*% vcov(obj) %*% t(a1)))
	low1 <- eff1 - qt(.975,obj$df.residual)*se.eff1
	up1 <- eff1 + qt(.975,obj$df.residual)*se.eff1
	eff2 <- a2 %*% obj$coef
	se.eff2 <- sqrt(diag(a2 %*% vcov(obj) %*% t(a2)))
	low2 <- eff2 - qt(.975,obj$df.residual)*se.eff2
	up2 <- eff2 + qt(.975,obj$df.residual)*se.eff2
	
	if(!plot.type %in% c("pdf", "png", "eps", "screen")){print("plot type must be one of - pdf, png or eps")}
	else{
	if(plot.type == "pdf"){pdf(paste(name.stem, "_", v1, ".pdf", sep=""), height=6, width=6)}
	if(plot.type == "png"){png(paste(name.stem, "_", v1, ".png", sep=""))}
	if(plot.type == "eps"){old.psopts <- ps.options(); setEPS(); postscript(paste(name.stem, "_", v1, ".eps", sep=""))}
	if(plot.type == "screen"){oldpar <- par(); par(mfrow=c(1,2))}
	plot(s2, eff1, type="n", ylim=range(c(low1, up1)), xlab=toupper(v2), ylab=paste("Conditional Effect of ", toupper(v1), " | ", toupper(v2), sep=""))
	if(par()$usr[3] < 0 & par()$usr[4] > 0){abline(h=0, col="gray50")}
	lines(s2, eff1)
	lines(s2, low1, lty=2)
	lines(s2, up1, lty=2)
	if(plot.type != "screen"){dev.off()}

	if(plot.type == "pdf"){pdf(paste(name.stem, "_", v2, ".pdf", sep=""), height=6, width=6)}
	if(plot.type == "png"){png(paste(name.stem, "_", v2, ".png", sep=""))}
	if(plot.type == "eps"){postscript(paste(name.stem, "_", v2, ".eps", sep=""))}
	plot(s1, eff2, type="n", ylim=range(c(low2, up2)), xlab=toupper(v1), ylab=paste("Conditional Effect of ", toupper(v2), " | ", toupper(v1), sep=""))
	if(par()$usr[3] < 0 & par()$usr[4] > 0){abline(h=0, col="gray50")}
	lines(s1, eff2)
	lines(s1, low2, lty=2)
	lines(s1, up2, lty=2)
	if(plot.type != "screen"){dev.off()}

	if(plot.type == "eps"){ps.options <- old.psopts}
	if(plot.type == "screen"){par <- oldpar}
	}
}

