searchVarLabels <- function(dat, str){
	if("var.labels" %in% names(attributes(dat))){
		vlat <- "var.labels"
	}
	if("variable.labels" %in% names(attributes(dat))){
		vlat <- "variable.labels"
	}
	ind <- grep(str, attr(dat, vlat), ignore.case=T)
	vldf <- data.frame(ind=ind, label = attr(dat, vlat)[ind])
	rownames(vldf) <- names(dat)[ind]
	vldf
}