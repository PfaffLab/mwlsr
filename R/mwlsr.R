#' mwlsr
#' 
#' Multiple Weighted Least Squares Regression (mwlsr). Used to fit gaussian
#' glm against multiple responses simultaneously.
#' 
#' @param data Input response matrix with responses in columns
#' @param design Design matrix. See \link{model.matrix}
#' @param weights Weights matrix
#' @param scale.weights If TRUE then weights are scaled (default behavior)
#' @param data.err Additional per-response-value uncertainty that should be
#' considered in the final sum of squared residual. Useful if your response
#' values have some knowm measurement uncertainty that you'd like to 
#' have considered in the models.
#' @param coef.method Method used to compute coefficients. This setting is
#' passed to \link{mols.coefs} or \link{wls.coefs}
#' @param coef.tol Tolerance setting for svd based coefficient calculation.
#' Passed to \link{mols.coefs} or \link{wls.coefs}
#' @param coefs.only Stop at the coefficient calculation and return only 
#' the coefficients of the models.
#' 
#' @return List with the following elements:
#' \item{coefficients}{Model coefficients}
#' \item{residuals}{Residuals of the fit}
#' \item{fitted.values}{Fitted values. Same dimension as the input response matrix.}
#' \item{deviance}{Sum of squared residuals}
#' \item{dispersion}{deviance / df.residual}
#' \item{null.deviance}{Sum of squared residuals for the NULL model (intercept only)}
#' \item{weights}{Weights matrix}
#' \item{prior.weights}{Weights matrix pre-scaling}
#' \item{weighted}{TRUE if fit was a weighted fit}
#' \item{df.residual}{Degrees of freedom of the model. \code{nrows(data) - ncol(design)}}
#' \item{df.null}{Degrees of freedom of the null model. \code{nrows(data) - 1}}
#' \item{y}{Input data matrix}
#' \item{y.err}{Input \code{data.err} matrix}
#' \item{X}{Design matrix}
#' \item{x}{If design matrix was based on factor levels then this will be a 
#' factor vector that matches the original grouping vector}
#' \item{intercept}{TRUE if the fit has an Intercept}
#' \item{coef.hat}{If the fit has an Intercept then this is a matrix of 
#' modified coefficients that represent the per-group averages. This is
#' calculated by adding the Intercept coefficients to each of the other 
#' coefficients. This only makes sense if your design was based on a single
#' multi-level factor}
#' 
#' @import parallel doParallel foreach doRNG 
#' @export
#' @examples
#' # Using the iris data.
#' design <- model.matrix(~Species, data=iris)
#' fit <- mwlsr(iris[, 1:4], design)
#' # test data association with the Species factor
#' result <- mwlsr.Ftest(fit)
#' print(table(result$F.padj < 0.05))
mwlsr <- function(data, design, weights=NULL, scale.weights=TRUE, data.err=NULL, 
	coef.method=c("chol", "ginv", "svd", "qr"), coef.tol=1e-7, coefs.only=FALSE) {
	
	# check parameters
	if(!inherits(data, "matrix")) {
		if(inherits(data, "data.frame")) {
			data <- as.matrix(data)
		} else {
			data <- matrix(data)
			colnames(data) <- "response"
		}
	}
	if(!missing(weights)) {
		if(!inherits(weights, "matrix")) {
			weights <- matrix(weights)
		}
	}
	if(!missing(data.err)) {
		if(!inherits(data.err, "matrix")) {
			data.err <- matrix(data.err)
		}
	}
	
	if(nrow(data) != nrow(design)) {
		stop("Design and data do not match")
	}
	
	# initalize some variables
	coef.method <- match.arg(coef.method)
	n <- nrow(design)
	p <- ncol(design)
	num.fits <- ncol(data)
	df.null <- n
	df.residual <- n-p
	use.weights <- FALSE
	intercept <- FALSE
	coef.names <- colnames(design)
	
	# check for intercept in the design
	if(grepl("intercept", coef.names[1], ignore.case=TRUE)) {
		intercept <- TRUE
		df.null <- n-1
	}

	# user provided prior variances per response value to propagate into 
	# the model's residuals
	if(missing(data.err)) {
		data.err <- data*0
	}
	
	# weights?
	if(is.null(weights)) {
		# no weights - make a weights matrix of 1's
		weights0 <- weights <- matrix(1, ncol=ncol(data), nrow=nrow(data))
	} else {
		if(!all.equal(dim(data), dim(weights))) {
			stop("Weights matrix doesn't match data dimension")
		}
		weights0 <- weights
		if(scale.weights) {
			# normalize by mean
			weights <- sweep(weights, 2, colMeans(weights), "/")
		}
		use.weights <- TRUE
	}

	if(use.weights) {
		# calculate weighted coefficients. this must run once per 
		# model or, if I was interested in streamlining it, finding groups 
		# of models that may share weights and calculating them in 
		# groups. 
		coefficients <- matrix(0, ncol=num.fits, nrow=ncol(design))
		wcoefficients <- coefficients
#		for(i in 1:num.fits) {
#			y <- data[, i]
#			w <- weights[, i]
#			# call wls.coefs
#			b <- drop(wls.coefs(design, y, weights=w, method=coef.method, tol=coef.tol))
#			coefficients[, i] <- b
#		}
		
		if(!require(parallel)) stop("Missing package 'parallel'")
		if(!require(doParallel)) stop("Missing package 'doParallel'")
		if(!require(foreach)) stop("Missing package 'foreach'")
		if(!require(doRNG)) stop("Missing package 'doRNG'")
		
		ncores <- detectCores() / 2 - 1
		cat("creating", ncores, "core cluster for fitting models\n")
		cl <- makeCluster(ncores)
		doParallel::registerDoParallel(cl)

		rres <- tryCatch({
			foreach(i=1:num.fits, .export=c("wls.coefs")) %dorng% { 
				y <- data[, i]
				w <- weights[, i]
				# call wls.coefs
				b <- drop(wls.coefs(design, y, weights=w, method=coef.method, tol=coef.tol))
				b
			}
		}, error=function(e) {
			# stop cluster
			# parallel::stopCluster(cl)
			print(paste("MY_ERROR: ", e))
			NULL
		}, finally = {
			stopCluster(cl)
		})
		
		if(is.null(rres)) stop("Something went wrong calculating coefficients in parallel.")
		
		coefficients <- do.call(cbind, rres)

	} else {
		# no weights so we can calc the coefficients in one shot
		coefficients <- mols.coefs(design, data, method=coef.method, tol=coef.tol)
	}

	if(coefs.only) {
		return(coefficients)
	}
	
	coefsHat <- NULL
	if(intercept && ncol(design) > 1) {	
		# if we have an intercept we can create a version of the coefficients that 
		# represents the condition means rather than the intercept + relative means
		coefsHat <- coefficients	
		for(i in 2:nrow(coefsHat)) {
			# add intercept
			coefsHat[i, ] <- coefsHat[i, ]+coefsHat[1, ]
		}
	}
	
	
	# calculate "fit" and residuals
	fitted.values <- design %*% coefficients
	residuals <- data-fitted.values

	# calculate null deviance with prior variances
	if(intercept) {
		# null deviance is relative to weighted mean
		null.deviance <- colSums(weights * sweep(data, 2, colSums(weights*data)/colSums(weights), "-")^2)
	} else {
		# no intercept so the null deviance is relative to 0
		null.deviance <- colSums(weights * data^2)
	}
	# not entirely sure if the prior errors need to be weighted...but it kinda makes sense
	null.deviance.err <- colSums(weights * data.err)
	null.deviance <- null.deviance + null.deviance.err
	
	# calculate residual deviance with prior variances - prior variances can't be
	# explained by the model so they get added in on top of the residuals
	deviance <- colSums(weights * residuals^2)
	deviance.err <- null.deviance.err 
	deviance <- deviance + deviance.err
	
	# final dispersion per fit
	wfactor <- df.residual*colSums(weights)/n
	dispersion <- deviance/wfactor
	
	# annotate all of the matrices and vectors
	dimnames(residuals) <- dimnames(fitted.values) <- dimnames(data)
	rownames(coefficients) <- colnames(design)
	colnames(coefficients) <- colnames(data)
	names(dispersion) <- names(deviance) <- names(null.deviance) <- colnames(data)
	
	# extract factor from design matrix if it's that kinda model
	x <- NULL
	if(ncol(design) > 1) {
		try(x <- mwlsr.design2factor(design))
	}
	
	# build the output list so it kinda resembles the output list of glm
	lout <- list(
			coefficients=coefficients, 
			residuals=residuals, 
			fitted.values=fitted.values, 
			deviance=deviance, 
			dispersion=dispersion,
			null.deviance=null.deviance,
			weights=weights, 
			prior.weights=weights0,
			weighted=use.weights,
			df.residual=df.residual, 
			df.null=df.null, 
			y=data, y.err=data.err, X=design, x=x, 
			intercept=intercept)
	
	if(intercept) {
		lout$coef.hat <- coefsHat
	}
	
	class(lout) <- c("mwlsr", class(lout))
	
	return(lout)
	
}

#' print.mwlsr
#' 
#' Override of generic \link{print} method for mwlsr objects.
#' 
#' @export
print.mwlsr <- function(x, ...) {
	
	cat("\n")
	cat("Multiple LS regression result (list)\n")
	
	cat("\n")
	cat("Members:\n")
	print(names(x))
	
}


#' mwlsr.rSquared
#' 
#' Calculate r-squared for each model in fit
#' 
#' @param fit Result of mwlsr fit
#' @return Input mwlsr object (list) with \code{rquared} element attached
#' @export
#' 
mwlsr.rSquared <- function(fit) {
	
	# calculate r^2 for the fit and attach it to the object
	rsq <- 1-fit$deviance/fit$null.deviance
	fit$rsquared <- rsq
	return(fit)
	
}


#' mwlsr.Fstatistic
#' 
#' Calculate F-statistic for each model.
#' 
#' @param fit mwlsr fit object
#' @return Input mwlsr object with \code{F} and \code{F.pval} elements
#' attached
#' @export
#' 
mwlsr.Fstatistic <- function(fit) {
	
	if(is.null(fit$rsquared)) {
		fit <- mwlsr.rSquared(fit)
	}
	
	dff <- nrow(fit$X)-1
	
	fit$F <- with(fit, (rsquared*df.residual)/((1-rsquared)*(dff-df.residual)))
	fit$F.pval <- with(fit, pf(F, dff-df.residual, df.residual, lower.tail=FALSE))
	
	return(fit)
	
}


#' mwlsr.Ftest
#' 
#' Calculates F-statistic and p-values for all models in the fit. Returns 
#' a table of the results.  This is only sensible if your design included
#' an intercept.
#' 
#' @param fit mwlsr fit object
#' @return data.frame with F-test results
#' @export
mwlsr.Ftest <- function(fit) {

	if(is.null(fit$F)) {
		fit <- mwlsr.Fstatistic(fit)
	}

	vid <- colnames(fit$coefficients)
	if(is.null(vid)) {
		vid <- 1:ncol(fit$coefficients)
	}

	df <- data.frame(varid=vid, df.null=fit$df.null, df=fit$df.residual,
		null.deviance=fit$null.deviance, deviance=fit$deviance, change=fit$null.deviance-fit$deviance, 
		r2=fit$rsquared, F=fit$F, F.pval=fit$F.pval, F.padj=p.adjust(fit$F.pval, method="BH"))

	rownames(df) <- vid

	return(df)

}


#' mwlsr.overallFstatistic
#' 
#' Calculates F-statistic and p-value for all models. In this test we 
#' sum all of the residual deviance and compare it to the total sum of 
#' null deviance.
#' 
#' @param fit mwlsr fit object
#' @return Vector containing the results
#' @export
mwlsr.overallFstatistic <- function(fit) {

	dd <- sum(fit$deviance)
	dd0 <- sum(fit$null.deviance)

	F <- ((dd0-dd)/(fit$df.null-fit$df.residual))/(dd/fit$df.residual)
	pval <- pf(F, fit$df.null-fit$df.residual, fit$df.residual, lower.tail=FALSE)

	cout <- c(fit$df.null, fit$df.residual, dd0, dd, dd0-dd, F, pval)
	names(cout) <- c("df.null", "df.residual", "null.deviance", "deviance", "change", "F", "pval")

	return(cout)

}


#' mwlsr.coefStats
#' 
#' Calculate coefficient standard errors, t-values and p-values. Results
#' are appended to the input mwlsr object.
#' 
#' @param fit mwlsr fit object
#' @return mwlsr fit object with coefficient statistics results appended
#' @export
mwlsr.coefStats <- function(fit) {

	# if weights then we have to calculate all of the weighted pseudo-inverse vectors
	if(fit$weighted) {

		n <- ncol(fit$coefficients)
		sccm <- sapply(1:n, function(i) {
			xhat <- t(fit$X) %*% diag(fit$weights[, i]) %*% fit$X
			return(diag(chol2inv(chol(xhat))))
		})
		fit$sccm <- sccm
		fit$coef.stderr <- sqrt(sccm %*% diag(fit$dispersion))
		fit$coef.tvals <- fit$coefficients/fit$coef.stderr
		fit$coef.pvals <- pt(abs(fit$coef.tvals), fit$df.residual, lower.tail=FALSE)*2

	} else {
		sccm <- diag(chol2inv(chol(crossprod(fit$X))))
		fit$sccm <- sccm
		fit$coef.stderr <- sqrt(sccm %*% matrix(fit$dispersion, 1))
		fit$coef.tvals <- fit$coefficients / fit$coef.stderr
		fit$coef.pvals <- pt(abs(fit$coef.tvals), fit$df.residual, lower.tail=FALSE)*2
	}

	dimnames(fit$coef.stderr) <- dimnames(fit$coef.tvals) <- dimnames(fit$coef.pvals) <- dimnames(fit$coefficients)
	
	return(fit)

}

#' mwlsr.groupStats
#' 
#' If your design was based on a single multi-level factor then you can
#' use this function to calculate per-group deviance, variance (dispersion), 
#' and standard error. Can be handy if you need to calculate group 
#' level variances.
#' 
#' @param fit mwlsr fit object
#' @return mwlsr fit object with group-level statistics appended
#' @importFrom MASS ginv
#' @export
mwlsr.groupStats <- function(fit) {
	
	wresid <- fit$weights * fit$residuals^2
	werr <- fit$weights * fit$y.err
	
	# sum of weighted squared residuals per group plus the sum of the weighted 
	# errors per group
	group.deviance <- (t(fit$X) %*% wresid) + (t(fit$X) %*% werr)
	group.rel <- (t(fit$X) %*% wresid)/group.deviance
	# make group dispersions
	n <- colSums(fit$X)-1
	if(any(n==0)) {
		n[n==0] <- 1
	}
	group.dispersion <- ginv(diag(n)) %*% group.deviance
	
	n <- colSums(fit$X)
	group.stderr <- ginv(diag(n)) %*% group.dispersion
	rownames(group.rel) <- rownames(group.stderr) <- rownames(group.deviance) <- rownames(group.dispersion) <- colnames(fit$X)
	
	fit$group.deviance <- group.deviance
	fit$group.dispersion <- group.dispersion
	fit$group.stderr <- group.stderr
	fit$group.rel <- group.rel
	
	return(fit)
	
}

#' mwlsr.contrastModelMatrix
#' 
#' Transforms a design matrix to a contrast matrix. Instead of calculating
#' the contrast result directly, as in \link{mwlsr.contrastTest}, 
#' this method would be used to create "reduced" model with a contrast
#' matrix and then you can evaluate significant associations with and LRT.
#' The code for this function is based on code within edgeR::glmLRT.
#' 
#' @param design Full design matrix
#' @param contrast Single column contrast matrix indicating which levels
#' of the full design to contrast against one other.
#' @return New, reduced, design matrix
#' 
#' @export
mwlsr.contrastModelMatrix <- function(design, contrast) {
	##
	# NOTE: the bulk, if not all, of this code is from edgeR::glmLRT
	
	contrast0 <- contrast
	contrast <- as.matrix(contrast)
	if(nrow(contrast) != ncol(design)) stop("contrast does not match design matrix dimension")
	
	coef.names <- colnames(design)
	nlibs <- ncol(design)
	
	qrc <- qr(contrast)
	ncontrasts <- qrc$rank
	if(ncontrasts==0) stop("contrasts are all zero")
	
	coef <- 1:ncontrasts
	
	if(ncontrasts > 1) {
		coef.name <- paste("LR test on", ncontrasts, "degrees of freedom")
	} else {
		contrast <- drop(contrast)
		i <- contrast != 0
		coef.name <- paste(paste(contrast[i], coef.names[i], sep="*"), collapse=" ")
	}
	
	Dvec <- rep.int(1, nlibs)
	Dvec[coef] <- diag(qrc$qr)[coef]
	Q <- qr.Q(qrc, complete=TRUE, Dvec=Dvec)
	design <- design %*% Q
	
	design0 <- design[, -coef]
	
	colnames(design0) <- paste("coef", 1:ncol(design0), sep="")
	attr(design0, "contrast") <- contrast0
	attr(design0, "coef.name") <- coef.name
	
	return(design0)
	
}

#' mwlsr.contrastTest
#' 
#' Contrast one or more levels of the design factors against one or more
#' other levels. This kind of test uses the full model's dispersions
#' as a basis for the comparison of the means of two conditions. Useful for
#' comparing levels of a single multi-level factor, such as different 
#' groups in an RNA-Seq experiment, to one another to check for statistical
#' difference in means.
#' 
#' @param fit mwlsr fit object
#' @param contrast Single-column contrast matrix
#' @param coef Coefficient to test (for intercept models)
#' @param ncomps Sidak post-hoc correction factor. Default behavior is
#' no correction.
#' @param squeeze.var Employ limma's 'squeeze.var' method which not only
#' adjusts the model's dispersion but also the residual degrees of freedom
#' @return data.frame with results of the test
#' @importFrom limma squeezeVar
#' 
#' @export
mwlsr.contrastTest <- function(fit, contrast=NULL, coef=NULL, ncomps=NULL, 
		squeeze.var=FALSE) {
	
	if(is.null(coef) & is.null(contrast)) {
		stop("Must specify either a coefficient or a contrast")
	}

	# figure out number of comparisons within model for Sidak correction
	if(is.null(ncomps)) {
		nconds <- length(levels(fit$x))
		ncomps <- nconds*(nconds-1)/2
	}
	
	if(squeeze.var) {
		if(!require(limma)) stop("Missing package 'limma'. Cannot perform 'squeeze.var' without it.")
		out <- limma::squeezeVar(fit$dispersion, fit$df.residual)
		fit$df.prior <- fit$df.residual
		fit$df.residual <- fit$df.residual + out$df.prior
		fit$dispersion.prior <- fit$dispersion
		fit$dispersion <- out$var.pos
	}

	if(!missing(coef)) {
		# return single coefficient statistics

		if(coef > ncol(fit$X)) {
			stop("Coefficient is beyond the design dimension")
		}
		if(is.null(fit$coef.stderr) || is.null(fit$coef.tvals) || is.null(fit$coef.pvals)) {
			message("Calculating coefficient statistics...")
			fit <- mwlsr.coefStats(fit)
		}

		mref <- fit$coefficients[1, ]
		mtarget <- mref+fit$coefficients[coef, ]
		tnum <- fit$coefficients[coef, ]
		tstat <- fit$coef.tvals[coef, ]
		tdenom <- fit$coef.stderr[coef, ]
		pval <- mwlsr.p.adjust(fit$coef.pvals[coef, ], n.comps=ncomps, method="sidak")
		baseMean <- (mref+mtarget)/2

	} else if(!missing(contrast)) {
		# return contrast test statistics
		
		coefs <- fit$coefficients
		design <- fit$X
		if(fit$intercept) {
			# use adjusted coefficients if we had an intercept
			coefs <- fit$coef.hat
			# if it is a factor based design..
			if(all(design==0 | design==1)) {
				idx_fix <- which(rowSums(design) > 1)
				design[idx_fix, 1] <- 0
			}
		}
		
		contrast <- matrix(drop(contrast))
		n <- ncol(fit$coefficients)

		# numerator for t stat also the change between the 
		# conditions being contrasted
		tnum <- drop(t(contrast) %*% coefs)

		if(fit$weighted) {	
			# weighted - one calculation per model		
			tdenom <- sapply(1:n, function(i) {
				wt <- fit$weights[, i]
				sscm <- chol2inv(chol((t(design) %*% diag(wt) %*% design)))
				rres <- t(contrast) %*% sscm %*% contrast
				return(rres)
			})
		} else {
			# no weights so we can do this more fasterer
			tdenom <- drop(t(contrast) %*% chol2inv(chol(crossprod(design))) %*% contrast)
		}
		# finish the standard error calculation
		tdenom <- sqrt(fit$dispersion * tdenom)

		# tstaistic and pvalue
		tstat <- tnum/tdenom
		pval <- pt(abs(tstat), fit$df.residual, lower.tail=FALSE)*2
		pval <- mwlsr.p.adjust(pval, n.comps=ncomps, method="sidak")

		# make the group means
		tmp <- contrast > 0
		ctmp <- contrast
		ctmp[tmp] <- 0
		mref <- drop(t(abs(ctmp)) %*% coefs)

		tmp <- contrast < 0
		ctmp <- contrast
		ctmp[tmp] <- 0
		mtarget <- drop(t(ctmp) %*% coefs)

		baseMean <- drop(t(abs(contrast)) %*% coefs)/2
				
	}

	# output table
	dres <- data.frame(row.names=colnames(fit$coefficients), id=colnames(fit$coefficients), 
		baseMean=baseMean, condA=mref, condB=mtarget, change=tnum, 
		stderr=tdenom, tstat=tstat, pval=pval, padj=p.adjust(pval, method="BH"), stringsAsFactors=FALSE)

	dres$status <- "n.s."
	tmp <- sapply(dres$padj, pval2stars)
	mm <- tmp != "N.S."
	if(any(mm)) {
		mhat <- mm & dres$change < 0
		if(any(mhat)) {
			dres$status[mhat] <- paste("sig.neg", tmp[mhat], sep="")
		}
		mhat <- mm & dres$change > 0
		if(any(mhat)) {
			dres$status[mhat] <- paste("sig.pos", tmp[mhat], sep="")
		}
		
	}

	if(!missing(coef)) {
		names(dres)[3:4] <- colnames(fit$X)[c(1, coef)]
	}
	
#	lout <- list(result=dres)
#	class(lout) <- c("mwlsrContrastResult", "list")

	# coming soon!
#	return(newDEResult(lout))
	
	return(dres)	
}


#' mwlsr.p.adjust
#' 
#' Performs multi-contrast within model p-value adjustment. This adjustment
#' is performed at each p-value and is based on the number of contrasts 
#' being tested in the model.
#' 
#' @param p Vector of p-values
#' @param n.comps Number of contrasts being tested within the model. Defaults
#' to the total number of possible pairwise contrasts (probably excessive)
#' @param n.levels Number of levels in the factor design (or columns of the 
#' design matrix.
#' @param n.samples Required for "scheffe" method.
#' @param method Adjustment method. Sidak is the default.
#' 
#' @return Adjusted p-values
#' 
#' @export
mwlsr.p.adjust <- function(p, n.comps=NA, n.levels=NA, n.samples=NA, method=c("sidak", "scheffe", "bonferroni")) {
	
	method <- match.arg(method)
	
	if(is.na(n.comps) & is.na(n.levels)) {
		stop("You must specifiy either the number of comparisons (n.comp) or the number of factor levels in the model (n.levels)")
	}
	
	if(method=="sidak" | method=="bonferroni") {
		if(is.na(n.comps)) {
			n.comp <- n.levels*(n.levels-1)/2
		}
	} else if(method=="scheffe") {
		if(is.na(n.levels)) {
			stop("Cannot calculate scheffe correction without total number of levels (n.levels)")
		}
		if(is.na(n.samples)) {
			stop("Cannot calculate scheffe correction without total number of samples (n.samples)")
		}
	}
	
	if(method=="scheffe") {
		if(!(any(p > 1))) {
			message("WARNING: Input for scheffe correction is supposed to be the LSD (t) statistic")
		}
	}
	
	p0 <- switch(method, 
			sidak={
				# this is just a transform by scaling
				1 - (1 - p)^n.comps
			}, 
			scheffe={
				F <- p^2/(n.levels-1)
				df1 <- n.levels-1
				df2 <- n.samples-n.levels
				# return p-values for the F-statistic
				pf(F, df1, df2, lower.tail=FALSE)
			}, 
			bonferroni={
				p*n.comps
			})
	
	if(any(p0==0)) {
		idx <- p0==0
		# set to smallest double such that 1-x != 1
		p0[idx] <- .Machine$double.neg.eps
	}
	
	return(p0)
	
}


#' mwlsr.tukeyHSD
#' 
#' Implementation of Tukey HSD per model. This only works for intercept 
#' designs.
#' 
#' @param fit mwlsr fit object
#' @return list with everything in it (TODO: explain results)
#' 
#' @export
mwlsr.tukeyHSD <- function(fit) {

	X <- fit$X
	if(!all(X==1 | X==0)) {
		stop("Unsure how to apply Tukey to this design")
	}

	terms <- colnames(X)
	if(fit$intercept) {
		terms[1] <- "Intercept"
	}

	f <- fit$x
	means <- fit$coefficients
	if(fit$intercept) {
		# if intercept then add the intercept to all of the other 
		# rows so that each row is now a term mean
		for(i in 2:nrow(means)) {
			means[i, ] <- means[i, ] + means[1, ]
		}
	}
	flevels <- levels(f)
	nn <- table(f)
	df.residual <- fit$df.residual
	MSE <- fit$dispersion

	pares <- combn(1:nrow(means), 2)
	center <- t(apply(pares, 2, function(x) means[x[2], ]-means[x[1], ]))

	onn <- apply(pares, 2, function(x) sum(1/nn[x]))
	SE <- t(sqrt(matrix(MSE/2) %*% matrix(onn, 1)))

	width <- qtukey(0.95, nrow(means), df.residual) * SE
	est <- center/SE
	pval <- ptukey(abs(est), nrow(means), df.residual, lower.tail=FALSE)

	# setup condition comparison labels
	lab0 <- apply(combn(1:length(flevels), 2), 2, function(x) flevels[rev(x)])
	lab <- apply(lab0, 2, function(x) paste(x, collapse="-"))
	
	# setup variable labels
	vid <- colnames(fit$coefficients)
	if(is.null(vid)) {
		vid <- as.character(1:ncol(fit$coefficients))
	}

	# setup 95% ci boundaries
	lower <- center - width
	upper <- center + width

	# build a list of tables - one for each comparison
	lres <- vector(mode="list", length=length(lab))
	names(lres) <- lab

	for(i in 1:length(lab)) {
		df <- data.frame(id=vid, change=center[i, ], lower=lower[i, ], 
			upper=upper[i, ], pval=pval[i, ], stringsAsFactors=FALSE)
		names(df)[2] <- paste(lab[i], "change", sep=".")
		lres[[i]] <- df
	}

	# build a list with everything
	lout <- list(results=lres, change=center, lower=lower, upper=upper, pval=pval)

	return(lout)
}

#' mols.coefs
#' 
#' Multiple ordinary least squares coefficients. Used interally by
#' \link{mwlsr} to compute coefficients without weights.
#' 
#' @param x Design matrix
#' @param y Response matrix
#' @param method Coefficient calculation method. \code{chol} is the fastest
#' and \code{svd} is said to be the most reliable but maybe the slowest.
#' @param tol Tolerance setting for the \code{svd} method.
#' @return Matrix of fit coefficients.
#' @importFrom MASS ginv
#' @export
mols.coefs <- function(x, y, method=c("chol", "ginv", "svd", "qr"), tol=1e-7) {

	method <- match.arg(method)

	if(!inherits(x, "matrix")) {
		stop("Expected x to be a design matrix such as the output of model.matrix")
	}

	if(!inherits(y, "matrix")) {
		y <- matrix(y)
	}
		
	ydim <- nrow(y)
	
	# check dimensions of things
	if(nrow(x) != ydim) {
		stop("response dimension doesn't match design")
	}
		
	if(method=="qr") {
		coefs <- qr.solve(x, y)
	} else if(method=="svd") {
		# use the SVD method - the most robust! this one doesn't care
		# if your design matrix isn't invertable because all that has to be
		# inverted are all eigenvalues > tol. it will, however, drop
		# coefficients out if they correspond to a eigenvector with 
		# < tol variance
		s <- svd(x)
		r <- max(which(s$d > tol))
		v1 <- s$v[, 1:r]
		sr <- diag(s$d[1:r])
		u1 <- s$u[, 1:r]
		coefs <- v1 %*% (ginv(sr) %*% t(u1) %*% y)
	} else {
		# use either the chol method (fastest of all) or the 
		# classic pseudoinverse with the ginv function. 
		# ginv failes less frequently than 'solve'
		sccm <- crossprod(x)
		if(method=="chol") {
			tmp <- chol2inv(chol(sccm))
		} else if(method=="ginv") {
			tmp <- ginv(sccm)
		}
		coefs <- tmp %*% t(x) %*% y
	}

	# annotate the rows and columns
	colnames(coefs) <- colnames(y)
	rownames(coefs) <- colnames(x)

	return(coefs)	

}

#' wls.coefs
#' 
#' Calculates coefficients for a single response with weights.
#' 
#' @param x Design matrix
#' @param y Response vector
#' @param weights Weights vector.
#' @param method Coefficient calculation method. See \link{mols.coefs}.
#' @param tol Tolerance for \code{svd} coefficient method.
#' @return Vector of fit coefficients
#' @importFrom MASS ginv
#' @export
wls.coefs <- function(x, y, weights=NULL, method=c("chol", "ginv", "svd", "qr"), tol=1e-7) {
	
	method <- match.arg(method)

	if(!inherits(x, "matrix")) {
		stop("Expected x to be a design matrix such as the output of model.matrix")
	}

	if(!inherits(y, "matrix")) {
		y <- matrix(y)
	}
	
	if(ncol(y) > 1) {
		stop("y matrix must only be a single response")
	}
	
	ydim <- nrow(y)
	
	# check dimensions of things
	if(nrow(x) != ydim) {
		stop("response dimension doesn't match design")
	}
	# check the weights out - don't normalize them leave that up to the 
	# calling codes
	if(!missing(weights)) {
		weights <- diag(drop(weights))
		if(ncol(weights) != ydim) {
			stop("weight dimension doesn't match response")
		}
	} else {
		# no weights - just make a diagonal of 1's so the below calculations
		# don't have to be changed
		weights <- diag(rep(1, nrow(y)))
	}
	
	# we can pre-weight x and y with the square root of the weights. this 
	# makes it so we can use the same x and y in all of the below
	# calculations
	wt <- weights^0.5
	xhat <- wt %*% x
	yhat <- wt %*% y
	
	if(method=="qr") {
		# use the QR method, second fastest
		wt <- weights^0.5
		coefs <- qr.solve(xhat, yhat)
	} else if(method=="svd") {
		# use the SVD method - the most robust! this one doesn't care
		# if your design matrix isn't invertable because all that has to be
		# inverted are all eigenvalues > tol. it will, however, drop
		# coefficients out if they correspond to a eigenvector with 
		# < tol variance
		s <- svd(xhat)
		r <- max(which(s$d > tol))
		v1 <- s$v[, 1:r]
		sr <- diag(s$d[1:r])
		u1 <- s$u[, 1:r]
		coefs <- v1 %*% (ginv(sr) %*% t(u1) %*% yhat)
	} else {
		# use either the chol method (fastest of all) or the 
		# classic pseudoinverse with the ginv function. 
		# ginv failes less frequently than 'solve'
		sccm <- crossprod(xhat)
		if(method=="chol") {
			tmp <- chol2inv(chol(sccm))
		} else if(method=="ginv") {
			tmp <- ginv(sccm)
		}
		coefs <- tmp %*% t(xhat) %*% yhat
	}

	coefs <- drop(coefs)
	names(coefs) <- colnames(x)
	return(coefs)
}

#
# turns a design matrix into a factor vector if the design matrix was
# based on a factor type model. it just fails otherwise.

#' mwlsr.design2factor
#' 
#' Derives a single multi-level factor vector from a design matrix. 
#' This would produce senseless results for regression models.
#' 
#' @param X Design matrix
#' @export
mwlsr.design2factor <- function(X) {
	
	if(!all(X==0 | X==1)) {
		return(NULL)
	}
	
	n <- ncol(X)
	m <- nrow(X)
	fac.levels <- colnames(X)
	fac <- character(m)
	if(sum(X[, 1])==m) {
		# intercept!
		fac <- rep("level0", m)
		for(i in 2:n) {
			idx <- which(X[, i]==1)
			fac[idx] <- fac.levels[i]
		}
	} else {
		for(i in 1:n) {
			idx <- which(X[, i]==1)
			fac[idx] <- fac.levels[i]
		}
	}
	
	return(factor(fac))
}


#' mwlsr.LRT
#' 
#' Perform likelihood ratio test (LRT) between two models (typically a 
#' full model and a reduced model). Using an F based LRT on a full 
#' model compared to an intercept-only model should give the same
#' results as \link{mwlsr.Ftest}.
#' 
#' @param full.m Full model mwlsr object
#' @param reduced.m Reduced model mwlsr object
#' @param test Type of test to perform. For gaussian models, which is
#' all that mwlsr can do, you should use the F-test.
#' 
#' @return data.frame with results of the test for all models.
#'  
#' @export
mwlsr.LRT <- function(full.m, reduced.m, test=c("F", "LRT")) {
	# f is the original model, f0 is the reduced model
	
	# make sure these are sorted
	ltests <- list(full.m, reduced.m)
	# resid
	dfresid <- sapply(ltests, function(x) x$df.residual)
	o <- order(dfresid)
	full.m <- ltests[[o[1]]]
	reduced.m <- ltests[[o[2]]]

	varnames <- colnames(full.m$coefficients)
	if(is.null(varnames)) {
		varnames <- as.character(1:ncol(full.m$coefficients))
	}

	test <- match.arg(test)
	deviance <- reduced.m$deviance - full.m$deviance
	df.test <- reduced.m$df.residual - full.m$df.residual
	dispersion <- full.m$deviance/full.m$df.residual
	
	if(test=="F") {
		LR <- (deviance/df.test)/dispersion
		pval <- pf(LR, df.test, full.m$df.residual, lower.tail=FALSE)
	} else if(test=="LRT") {
		LR <- deviance/dispersion
		pval <- pchisq(LR, df.test, lower.tail=FALSE)
		
		# edgeR defines the LR as what I have as deviance in this code. their
		# p-value looks like this:
		# pval <- pchisq(deviance, df.test, lower.tail=FALSE)
		
	}
	
	padj <- p.adjust(pval, method="BH")
	
	dout <- data.frame(variable=varnames, deviance=deviance, 
		df=rep(df.test, length(deviance)), 
		LR=LR, pval=pval, padj=padj, stringsAsFactors=FALSE)
	
	return(dout)
}


#' mwlsr.contrastCoefficients
#' 
#' Calculate the coefficients of a contrast fit.
#' 
#' @param fit mwlsr fit object
#' @param contrast Single-column contrast matrix
#' @return mwlsr fit object with results appended
#' 
#' @export
mwlsr.contrastCoefficients <- function(fit, contrast) {

	contrast0 <- contrast
	contrast <- drop(contrast)

	# split contrast vector into positive and negative
	pmask <- ifelse(contrast > 0, 1, 0)
	nmask <- ifelse(contrast < 0, -1, 0)

	coef.pos <- drop(matrix(contrast*pmask, 1) %*% fit$coefficients)
	coef.neg <- drop(matrix(contrast*nmask, 1) %*% fit$coefficients)

	cout <- rbind(coef.neg, coef.pos)
	colnames(cout) <- colnames(fit$coefficients)

	# combine names of levels
	names.pos <- paste(rownames(fit$coefficients)[which(pmask != 0)], collapse=":")
	names.neg <- paste(rownames(fit$coefficients)[which(nmask != 0)], collapse=":")
	rownames(cout) <- c(names.neg, names.pos)

	fit$contrast <- contrast0
	fit$contrasts.coefs <- cout
	
	return(fit)

}

#' mwlsr.makeContrast
#' 
#' Create a contrast matrix for computing statistical difference 
#' in means between levels of a design. For example if your model
#' is based on a single multi-level factor and you want to compare 
#' the average values of one level verses another within the context
#' of the full model then you'd specify those levels in this function
#' and then call on \link{mwlsr.contrastTest} to perform the test.
#' 
#' @param y First level or levels to contrast.
#' @param lvls Either the factor levels or a full length factor vector or
#' the design matrix from the mwlsr object or the mwlsr object itself.
#' @param x Factor level or levels to contrast against those specified in 
#' \code{y}. If this is omitted then the level or levels specified in 
#' \code{y} are contrasted against all other levels.
#' @return Single-column matrix specifying the contrast
#' @details Levels specified for \code{y} or \code{x} must match levels
#' in the design matrix associated with the model you plan to perform 
#' the contrast test within.
#' 
#' @export
mwlsr.makeContrast <- function(y, lvls, x=NULL) {
	
	if(inherits(lvls, "mwlsr")) {
		lvls <- colnames(lvls$X)
	} else if(inherits(lvls, "factor")) {
		lvls <- levels(factor)
	} else if(inherits(lvls, "matrix")) {
		# assuming this is a design matrix
		lvls <- colnames(lvls)
	} else if(is.null(dim(lvls)) & length(lvls) > 0) {
		lvls <- levels(factor(lvls))
	}
		
#	nvalid <- lvls != make.names(lvls)
#	if(any(nvalid)) {
#		stop("The levels must be valid names (use make.names)")
#	}
	
	tmp <- c(y, x)
	if(!all(tmp %in% lvls)) 
		stop("one or more of the specified levels for your contrast are not in the design")
	
	ny <- length(y)
	nyidx <- sapply(y, function(a) { match(a, lvls) })
	idx <- 1:length(lvls)
	
	if(!is.null(x)) {
		nx <- length(x)
		nxidx <- sapply(x, function(a) { match(a, lvls) })
	} else {
		nx <- length(lvls)-ny
		nxidx <- idx[-nyidx]
	}
	
	ci <- rep(0, length(lvls))
	ci[nyidx] <- 1/ny
	ci[nxidx] <- -1/nx
	
	# build contrast string
	o <- order(ci, decreasing=TRUE)
	ci.tmp <- ci[o]
	lvls.tmp <- lvls[o]
	
	i <- ci.tmp != 0
	sz <- paste(paste(round(ci.tmp[i], 4), lvls.tmp[i], sep="*"), collapse=" ")
	
	ci <- matrix(ci, dimnames=list(Levels=lvls, Contrast="coefs"))
	attr(ci, "contrast") <- sz
	attr(ci, "levels") <- lvls
	
	return(ci)
}


pval2stars <- function(p, max.stars=4) {
	
	psig <- p < 0.05
	rres <- ifelse(psig, "*", "N.S.")
	if(any(p==0)) {
		rres[p==0] <- paste(rep("*", max.stars+1), collapse="")
	}
	
	mm <- p <= 0.1
	plog <- -log10(p)
	if(any(mm)) {
		for(i in which(mm)) {
			rres[i] <- paste(rep("*", min(c(max.stars, plog[i]))), collapse="")
		}
	}
	
	return(rres)
	
}
