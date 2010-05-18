DAintfun2 <- function(obj, varnames, rug=TRUE, ticksize=-.03, hist=FALSE, hist.col="gray75", 
	nclass=c(10, 10), scale.hist=.5, border=NA, name.stem = "cond_eff", plot.type = "screen"){
	rseq <- function(x){rx <- range(x, na.rm=TRUE); seq(rx[1], rx[2], length=25)}
	if(!("model" %in% names(obj))){obj <- update(obj, model=T)}
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
	if(hist == TRUE){
		rng <- diff(par()$usr[3:4])
		h2 <- hist(obj$model[[v2]], nclass=nclass[1], plot=FALSE)
		prop2 <- h2$counts/sum(h2$counts)
		plot.prop2 <- (prop2/max(prop2))*rng*scale.hist + par()$usr[3]
		av2 <- pretty(prop2, n=3)
		axis(4, at=(av2/max(prop2))*rng*scale.hist + par()$usr[3], labels=av2)
		br2 <- h2$breaks
		for(i in 1:(length(br2)-1)){
			polygon(
				x=c(br2[i], br2[(i+1)], br2[(i+1)], br2[i], br2[i]),
				y=c(par()$usr[3], par()$usr[3], plot.prop2[i], plot.prop2[i], par()$usr[3]), 
				col=hist.col, border=border)
		}
	}
	if(rug == TRUE){rug(obj$model[[v2]], ticksize=ticksize)}
	if(par()$usr[3] < 0 & par()$usr[4] > 0){abline(h=0, col="gray50")}
	lines(s2, eff1)
	lines(s2, low1, lty=2)
	lines(s2, up1, lty=2)
	if(plot.type != "screen"){dev.off()}
	if(plot.type == "pdf"){pdf(paste(name.stem, "_", v2, ".pdf", sep=""), height=6, width=6)}
	if(plot.type == "png"){png(paste(name.stem, "_", v2, ".png", sep=""))}
	if(plot.type == "eps"){postscript(paste(name.stem, "_", v2, ".eps", sep=""))}
	plot(s1, eff2, type="n", ylim=range(c(low2, up2)), xlab=toupper(v1), ylab=paste("Conditional Effect of ", toupper(v2), " | ", toupper(v1), sep=""))
	if(hist == TRUE){
		rng <- diff(par()$usr[3:4])
		h1 <- hist(obj$model[[v1]], nclass=nclass[2], plot=FALSE)
		prop1 <- h1$counts/sum(h1$counts)
		plot.prop1 <- (prop1/max(prop1))*rng*scale.hist + par()$usr[3]
		av1 <- pretty(prop1, n=3)
		axis(4, at=(av1/max(prop1))*rng*scale.hist + par()$usr[3], labels=av1)
		br1 <- h1$breaks
		for(i in 1:(length(br1)-1)){
			polygon(
				x=c(br1[i], br1[(i+1)], br1[(i+1)], br1[i], br1[i]),
				y=c(par()$usr[3], par()$usr[3], plot.prop1[i], plot.prop1[i], par()$usr[3]), 
				col=hist.col, border=border)
		}
	}
	if(rug == TRUE){rug(obj$model[[v1]], ticksize=ticksize)}
	if(par()$usr[3] < 0 & par()$usr[4] > 0){abline(h=0, col="gray50")}
	lines(s1, eff2)
	lines(s1, low2, lty=2)
	lines(s1, up2, lty=2)
	if(plot.type != "screen"){dev.off()}

	if(plot.type == "eps"){ps.options <- old.psopts}
	if(plot.type == "screen"){par <- oldpar}
	}
}

