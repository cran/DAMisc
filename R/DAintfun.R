DAintfun <- function(obj, varnames, theta=45, phi=10){
require(sm)
require(effects)
    if(length(varnames) !=2){
        stop("varnames must be a vector of 2 variable names")
    }
    else{
        v1 <- varnames[1]
        v2 <- varnames[2]
                
        ind <- unique(c(grep(v1, names(obj$coef)), grep(v2, names(obj$coef))))
        ind <- ind[order(ind)]        
        b <- obj$coef[ind]
        mod.x <- model.matrix(obj)                
        not.ind <- c(1:ncol(mod.x))[!(c(1:ncol(mod.x))%in% ind)]
        mod.x[, not.ind] <- 0
        dens <- sm.density(mod.x[,varnames], display="none")
        b <- obj$coef[ind]
        v1.seq <- dens$eval.points[,1]
        v2.seq <- dens$eval.points[,2]
        eff.fun <- function(x1,x2){
            b[1]*x1 + b[2]*x2 + b[3]*x1*x2
        }
        hcols <- paste("gray", seq(from=20, to=80, length=4),sep="")
predsurf <- outer(v1.seq, v2.seq, eff.fun)
cutoff <- quantile(dens$estimate, prob=c(0.25, 0.5, 0.75))
pred1 <- predsurf
pred1[dens$estimate < cutoff[1]] <- NA
pred2 <- predsurf
pred2[dens$estimate < cutoff[2]] <- NA
pred3 <- predsurf
pred3[dens$estimate < cutoff[3]] <- NA
persp(v1.seq, v2.seq, predsurf, xlab = toupper(v1), ylab = toupper(v2), zlab=toupper("Predictions"), col=hcols[1], theta=theta, phi=phi)
par(new=TRUE)
persp(v1.seq, v2.seq, pred1, col=hcols[2], axes=FALSE, xlab="", ylab="", zlab="", theta=theta, phi=phi,
    zlim=c(min(c(predsurf)), max(c(predsurf))), ylim=c(min(v2.seq), max(v2.seq)), xlim=c(min(v1.seq), max(v1.seq)))
par(new=TRUE)
persp(v1.seq, v2.seq, pred2, col=hcols[3], axes=FALSE, xlab="", ylab="", zlab="", theta=theta,phi=phi,
    zlim=c(min(c(predsurf)), max(c(predsurf))), ylim=c(min(v2.seq), max(v2.seq)), xlim=c(min(v1.seq), max(v1.seq)))
par(new=TRUE)
persp(v1.seq, v2.seq, pred3, col=hcols[4], axes=FALSE, xlab="", ylab="", zlab="", theta=theta, phi=phi,
    zlim=c(min(c(predsurf)), max(c(predsurf))), ylim=c(min(v2.seq), max(v2.seq)), xlim=c(min(v1.seq), max(v1.seq)))
    }
}

