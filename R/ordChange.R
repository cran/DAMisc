ordChange <- function(obj, data, typical.dat=NULL){
vars <- names(attr(terms(obj), "dataClasses"))[-1]
pols <- grep("poly", vars)
if(length(pols) > 0){
	poly.split <- strsplit(vars[pols], split="")
	start <- lapply(poly.split, function(x)grep("(", x, fixed=T)+1)
	stop <- lapply(poly.split, function(x)grep(",", x, fixed=T)-1)
	pol.vars <- sapply(1:length(poly.split), function(x)paste(poly.split[[x]][start[[x]]:stop[[x]]], collapse=""))
	vars[pols] <- pol.vars
}
rn <- vars
var.classes <- sapply(vars, function(x)class(data[[x]]))
minmax <- lapply(vars, function(x)c(NA, NA))
meds <- lapply(vars, function(x)NA)
names(minmax) <- names(meds) <- vars
levs <- obj$xlevels
if(length(levs) > 0){
	for(i in 1:length(levs)){
		tmp.levs <- paste(names(levs)[i], unlist(levs[i]), sep="")
		col.inds <- match(tmp.levs, names(obj$coef))
		if(length(grep("1$", names(obj$coef)[col.inds])) > 0){
			col.inds <- c(col.inds[which(is.na(col.inds))], col.inds[grep("1$", names(col.inds))])
			names(col.inds) <- gsub("1$", "", names(col.inds))
			col.inds <- col.inds[match(tmp.levs, names(col.inds))]
		}
		tmp.coefs <- obj$coef[col.inds]
		tmp.coefs[which(is.na(tmp.coefs))] <- 0
		mm <- c(which.min(tmp.coefs), which.max(tmp.coefs))
		minmax[[names(levs)[i]]] <- factor(levs[[i]][mm], levels=levs[[i]])
		tmp.tab <- table(data[[names(levs)[i]]])
		meds[[names(levs)[i]]] <- factor(names(tmp.tab)[which.max(tmp.tab)], levels=levs[[i]])
	}
}

vars <- vars[sapply(minmax, function(x)is.na(x[1]))]
for(i in 1:length(vars)){
	minmax[[vars[i]]] <- range(data[[vars[i]]], na.rm=T)
	meds[[vars[i]]] <- median(data[[vars[i]]], na.rm=T)
}
tmp.df <- do.call(data.frame, lapply(meds, function(x)rep(x, length(meds)*2)))
if(!is.null(typical.dat)){
	notin <- which(!(names(typical.dat) %in% names(tmp.df)))
	if(length(notin) > 0){
		cat("The following variables in typical.dat were not found in the prediction data: ", names(typical.dat)[notin], "\n\n", sep="")
		typical.dat <- typical.dat[,-notin]
	}
for(j in 1:ncol(typical.dat)){
	tmp.df[[names(typical.dat)[j]]] <-  typical.dat[1,j]
	meds[names(typical.dat)[j]] <- as.numeric(typical.dat[1,j])
}
}
inds <- seq(1, nrow(tmp.df), by=2)
for(j in 1:length(minmax)){
	tmp.df[inds[j]:(inds[j]+1), j] <- minmax[[j]]
}
preds <- predict(obj, newdata=tmp.df, type="probs")
preds.min <- preds[seq(1, nrow(preds), by=2), ]
preds.max <- preds[seq(2, nrow(preds), by=2), ]
diffs <- preds.max-preds.min
rownames(preds.min) <- rownames(preds.max) <- rownames(diffs) <- rn
minmax.mat <- do.call(data.frame, minmax)
minmax.mat <- rbind(do.call(data.frame, meds), minmax.mat)
rownames(minmax.mat) <- c("typical", "min", "max")
ret <- list(diffs=diffs, minmax = minmax.mat, minPred = preds.min, maxPred=preds.max)
class(ret) <- "change"
return(ret)
}
