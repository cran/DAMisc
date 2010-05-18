ziChange <- function(obj, data, typical.dat=NULL, type="count"){
require(gdata)
vars <- as.character(formula(obj))[3]
vars <- gsub("*", "+", vars, fixed=T)
vars <- strsplit(vars, split="|", fixed=T)[[1]]
vars <- c(unlist(strsplit(vars, "+", fixed=T)))
vars <- unique(trim(vars))
pols <- grep("poly", vars)
if(length(pols) > 0){
	poly.split <- strsplit(vars[pols], split="")
	start <- lapply(poly.split, function(x)grep("(", x, fixed=T)+1)
	stop <- lapply(poly.split, function(x)grep(",", x, fixed=T)-1)
	pol.vars <- sapply(1:length(poly.split), function(x)paste(poly.split[[x]][start[[x]]:stop[[x]]], collapse=""))
	vars[pols] <- pol.vars
}
rn <- vars
vars.type <- as.character(formula(obj))[3]
vars.type <- gsub("*", "+", vars.type, fixed=T)
vars.type <- strsplit(vars.type, split="|", fixed=T)[[1]][ifelse(type == "count", 1, 2)]
vars.type <- c(unlist(strsplit(vars.type, "+", fixed=T)))
vars.type <- unique(trim(vars.type))
pols <- grep("poly", vars.type)
if(length(pols) > 0){
	poly.split <- strsplit(vars.type[pols], split="")
	start <- lapply(poly.split, function(x)grep("(", x, fixed=T)+1)
	stop <- lapply(poly.split, function(x)grep(",", x, fixed=T)-1)
	pol.vars.type <- sapply(1:length(poly.split), function(x)paste(poly.split[[x]][start[[x]]:stop[[x]]], collapse=""))
	vars.type[pols] <- pol.vars.type
}
var.classes <- sapply(vars, function(x)class(data[[x]]))
minmax <- lapply(vars, function(x)c(NA, NA))
meds <- lapply(vars, function(x)NA)
names(minmax) <- names(meds) <- vars
levs <- obj$levels
if(length(levs) > 0){
	for(i in 1:length(levs)){
		tmp.levs <- paste(names(levs)[i], unlist(levs[i]), sep="")
		tmp.coefs <- obj$coef[[type]][match(tmp.levs, names(obj$coef[[type]]))]
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
preds <- matrix(predict(obj, newdata=tmp.df, type=type), ncol=2, byrow=T)
diffs <- cbind(preds, apply(preds, 1, diff))
colnames(diffs) <- c("min", "max", "diff")
rownames(diffs) <- rn
diffs <- diffs[which(rownames(diffs) %in% vars.type),]
minmax.mat <- do.call(cbind, minmax)
minmax.mat <- rbind(c(unlist(meds)), minmax.mat)
rownames(minmax.mat) <- c("typical", "min", "max")
ret <- list(diffs=diffs, minmax = minmax.mat)
class(ret) <- "change"
return(ret)
}

