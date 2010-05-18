factorplot <-function(obj, factor.variable=NULL, pval=0.05, two.sided=TRUE, bonferroni=FALSE){
require(effects)
tmp.classes <- attr(terms(obj), "dataClasses")
tmp.classes <- tmp.classes[tmp.classes == "factor"]
tmp.levs <- NULL
for(i in 1:length(tmp.classes)){
    tmp.levs <- c(tmp.levs, length(levels(obj$model[[names(tmp.classes)[i]]])))
}
tmp.class <- tmp.classes[tmp.levs > 2]
if(is.null(factor.variable)){
{if(length(tmp.classes) > 1){
    myvar <- names(tmp.classes)[menu(names(tmp.classes))]
}
else{ 
    myvar <- names(tmp.classes[1])
}}
}
else{
	myvar <- factor.variable
}
bcols <- grep(myvar, names(coef(obj)))
b <- coef(obj)[bcols]
v <- vcov(obj)[bcols, bcols]

b.diff <- b.sd <- matrix(NA, ncol=length(bcols), nrow=length(bcols))
b.diff[1, ] <- (-1)*b
b.sd[1, ] <- sqrt(diag(v))
for(i in seq(from=1, to=(length(bcols)-1))){
    for(j in seq(from=i+1, to=length(bcols))){
        b.diff[(i+1),j] <- b[i]-b[j]
        b.sd[(i+1),j] <- sqrt(v[i,i] + v[j,j] + 2*(-1)*v[j,i])
        }
}
b.t <- b.diff/b.sd
r.bdiff <- b.diff[rev(1:nrow(b.diff)), ]
r.bt <- b.t[rev(1:nrow(b.t)), ]
r.bsd <- b.sd[rev(1:nrow(b.sd)), ]
rns <- abbreviate(rownames(contrasts(obj$model[[myvar]]))[1:length(bcols)], minlength=8)
cns <- abbreviate(colnames(contrasts(obj$model[[myvar]])), minlength=8)
rns.split <- strsplit(rns, split=" ")
rns.out <- NULL
for(i in 1:length(rns.split)){
    rns.out <- c(rns.out, paste(rns.split[[i]], collapse="\n "))
}
cns.split <- strsplit(cns, split=" ")
cns.out <- NULL
for(i in 1:length(cns.split)){
    cns.out <- c(cns.out, paste(cns.split[[i]], collapse="\n "))
}
rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- abbreviate(rns, minlength=8)
colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- abbreviate(cns, minlength=8)
pval.print <- pval
if(bonferroni){
	pval <- pval/((ncol(b.t)*(nrow(b.t)-1)) + nrow(b.t))
	pval.print <- formatC(pval/((ncol(b.t)*(nrow(b.t)-1)) + nrow(b.t)), digits=3, format="E")
}
par(mar=c(0,0,0,0))
plot(c((1-length(bcols)*.15),length(bcols)+1+(length(bcols)*.01)), 
    c(1, length(bcols)+1+(.15*length(bcols))), type="n", main="", xlab="", ylab="", axes=FALSE)
axis(3, at=seq(from=1.5, to=length(bcols)+.5, by=1), labels=cns.out, tick=FALSE, lwd=0)
axis(2, at=seq(from=1.5, to=length(bcols)+.5, by=1), labels=rns.out[rev(rns.out)], tick=FALSE, lwd=0)
rseq <- rev(1:length(bcols))
m <- 1
for(i in rseq){ 
    for(j in m:length(bcols)){
        colvec <- c("gray90", "gray70")
        col.ind <- as.numeric((2^as.numeric(two.sided))*pt(abs(b.t[m,j]), obj$df.residual, 
            lower.tail=FALSE) < pval) + 1
        squares(c(j,i), col=colvec[col.ind])
    text(j+.5, i+(.5+min(.02*length(bcols))), round(r.bdiff[i,j], 2), font=2, 
        cex=1-(.025*(length(bcols)-2)))
    text(j+.5, i+(.5-min(.02*length(bcols))), round(r.bsd[i,j], 2), font=3, 
        cex=1-(.025*(length(bcols)-2)))
    if(i == length(bcols)){
        text(j+.5, length(bcols)+1+(.09*length(bcols)), cns[j], srt=90)    
        }
    }
    text(1, i+.5, rns[m], pos=2)
m <- m+1
}
squares(c(1+(.1*length(bcols)), .95), 0.25, col="gray90")
squares(c(1+(.1*length(bcols)), 1.25), 0.25, col="gray70")
text(1.25+(.1*length(bcols)), 1.375, paste("p < ", pval.print, sep=""), pos=4, 
    cex=1-(.025*(length(bcols)-2)))
text(1.25+(.1*length(bcols)), 1.075, paste("p > ", pval.print, sep=""), pos=4, 
    cex=1-(.025*(length(bcols)-2)))
if(bonferroni){
text(1.25+(.3*length(bcols)), 1.3, "Bonferroni", pos=4, 
    cex=1-(.025*(length(bcols)-2)))
text(1.25+(.3*length(bcols)), 1.1, "Corrected", pos=4, 
    cex=1-(.025*(length(bcols)-2)))

}
text(1+(.15*length(bcols)), 1.5+.1*length(bcols), "bold", font=2, pos=2, 
    cex=1-(.025*(length(bcols)-2)))
text(1+(.15*length(bcols)), 1.5+.1*length(bcols), expression(b[row]-b[col]), font=2, pos=4, 
    cex=1-(.025*(length(bcols)-2)))
text(1+(.15*length(bcols)), 1.5+.05*length(bcols), "ital", font=3, pos=2, 
    cex=1-(.025*(length(bcols)-2)))
text(1+(.15*length(bcols)), 1.5+.05*length(bcols), expression(S.E.(b[row]-b[col])), 
    font=2, pos=4, cex=1-(.025*(length(bcols)-2)))
}
