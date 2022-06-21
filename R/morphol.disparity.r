#' Morphological disparity for one or more groups of specimens 
#'
#' Function estimates morphological disparity and performs pairwise comparisons among groups 
#'
#' The function estimates morphological disparity and performs pairwise comparisons to identify differences
#' among groups. Morphological disparity is estimated as the Procrustes variance, overall or for groups, 
#' using residuals of a linear model fit.  Procrustes variance is the same sum of the diagonal elements 
#' of the group covariance matrix divided by the number of observations in the group (e.g., Zelditch et al. 2012).
#' The function takes as input a formula to describe the linear model fit,
#' plus a formulaic indication of groups (e.g., ~ groups).  It is assumed that the formula describes shape data that 
#' have been GPA-aligned [e.g., \code{\link{gpagen}}], although the function can work with any multivariate data.
#' 
#' Partial disparities (Foote 1993) can also be calculated, but only for model formulas containing only an intercept 
#' (e.g., coords ~ 1).  Partial disparity has the same numerator as Procrustes variance (with respect to an overall mean)
#' but the denominator is N - 1 for all N observations, rather than n, the group size.  (The sum of all group n equals N.)
#' Partial disparities have the appeal that the sum of group partial disparities it the total disparity.
#' 
#' Absolute differences in Procrustes variances are test statistics that can be used to test differences
#' in morphological disparity among groups.  These differences are  statistically evaluated through permutation, 
#' where the vectors of residuals are randomized among groups. The function can be used to obtain disparity for the whole
#' dataset by using "a dummy group factor "~ 1" as the right-hand portion of the formula, in which case only Procrustes 
#' variance is returned.  Additionally, if the right-hand portion of the formula only contains (continuous) covariates, 
#' e.g., "~ Csize", Procrustes variance will be calculated for the whole data set or groups, after accounting for the
#' linear regression described.  Finally, different factors can be indicated in the formula and for groups, if one wishes to 
#' compare morphological disparities for groups comprising only a portion of or collapsing of the groups in a more complex model 
#' (see examples).
#' 
#' This function can be used with an object of class "procD.lm" or "lm.rrpp", if such analyses have already been performed.  
#' This is especially useful for analyses performed with  \code{\link{procD.pgls}}.  In this case, residuals obtained from PGLS estimation
#' of coefficients, rather than OLS estimation, will be used in the analysis.  Thus, one can account for phylogeny when comparing
#' morphological disparity among groups.  However, one should be aware that PGLS finds the PGLS residuals and transforms them via the 
#' phylogenetic covariance matrix in order to calculate variance.  Thus, disparity (dispersion in tangent space) and variance (dispersion
#' in transformed space, after accounting for phylogeny) can be considered different things.  This function uses an argument, 
#' \code{transform}, to allow users to use either GLS residuals in tangent space or transformed residuals in the transformed space to 
#' calculate disparity.  The former characterizes the dispersion of actual shapes around GLS means; the latter estimates GLS variances, 
#' by group.
#'
#' \subsection{Notes for geomorph 3.1.0 and subsequent versions}{ 
#'  The function \code{\link{pairwise}} in the \code{RRPP} package can also be used to evaluate morphological 
#'  disparity, and will yield results identical to those of the current function. A simple example is shown below
#'  }
#'
#' @param f1 A formula describing the linear model used.  The left-hand portion of the formula should be
#'  a 3D array (p x k x n) containing Procrustes shape variables for a set of specimens, or a matrix (n x variables). 
#'  The right-hand portion of the formula should be " ~1" to use the overall mean, or "~ x1 + x2 + x3 +...", where each x is a 
#'  covariate or factor.  (Interactions and nested terms also work.)  Alternatively, one can use an object of class "procD.lm" or
#'  "lm.rrpp", which already has a formula defined.  This is especially helpful for analyses performed with \code{\link{procD.pgls}},
#'  as residuals from PGLS will be used for analysis (see examples).
#' @param groups Either a formula designating groups, e.g., groups = ~ groups, or a factor, or 
#' a vector coercible to factor.  If NULL and if f1 is a procD.lm or lm.rrpp model fit, morphol.disparity
#' will attempt to define groups based on the terms of the model.  If there are no groups inherently
#' indicated in f1 and groups is NULL, a single Procrustes variance will be returned for the entire data set.
#' @param partial A logical value to indicate whether partial disparities should be calculated, sensu
#' Foote (1993).  If TRUE, the model formula should have only an intercept (e.g., coords ~ 1); otherwise an error 
#' will be returned.
#' @param transform A logical value to indicate whether disparity should be measured in the transformed
#' space rather than the tangent space.  This is only applicable if a PGLS model is used, in which case
#' the transformation would be necessary for a variance that is independent of phylogeny.  If transform = FALSE,
#' disparity tracks the dispersion one would observe in a phylo-PC plot.  If transform = TRUE, disparity tracks
#' the dispersion in a PC plot of the transformed space (lacking phylogenetic signal).  
#' See Collyer and Adams 2021, and examples below.
#' @param iter Number of iterations for permutation test
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}}
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @param ... Arguments passed on to procD.lm (typically associated with the lm function,
#' such as weights or offset, plus SS.type). 
#' @keywords analysis
#' @export
#' @author Emma Sherratt and Michael Collyer
#' @return Objects of class "morphol.disparity" return a list with the following components 
#' (if groups are specified): 
#'   \item{Procrustes.var}{Observed Procrustes variances.}
#'   \item{PV.dist}{Observed pairwise absolute differences (as distances) 
#'   among group Procrustes variances.}
#'   \item{PV.dist.Pval}{P-values associated with pairwise differences.}
#'   \item{random.PV.dist}{Pairwise distance matrices produced in the resampling procedure.}
#'   \item{permutations}{Number of random permutations in resampling procedure.}
#'   \item{partial}{Logical value to indicate if partial disparities were calculated.}
#'   \item{call}{The match call}
#'   
#' @references Zelditch, M. L., D. L. Swiderski, H. D. Sheets, and W. L. Fink. 2012. Geometric morphometrics 
#'   for biologists: a primer. 2nd edition. Elsevier/Academic Press, Amsterdam.
#' @references Foote, M. 1993. Contributions of individual taxa to overall morphological disparity. 
#' Paleobiology, 19: 403-419.
#' @references Collyer, M.L and D.C. Adams, 2021. Phylogenetically Aligned component analysis. 
#' Methods in Ecology and Evolution, 12: 369-372.
#' 
#' @examples
#' data(plethodon)
#' Y.gpa<-gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment
#' gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
#' site = plethodon$site)
#' 
#' # Morphological disparity for entire data set
#' morphol.disparity(coords ~ 1, groups = NULL, data = gdf, 
#' iter = 299, print.progress = FALSE)
#' 
#' # Morphological disparity for entire data set, accounting for allometry
#' morphol.disparity(coords ~ Csize, groups= NULL, data = gdf, 
#' iter = 299, print.progress = FALSE)
#' 
#' # Morphological disparity without covariates, using overall mean
#' morphol.disparity(coords ~ 1, groups= ~ species*site, data = gdf, 
#' iter = 299, print.progress = FALSE)
#' 
#' # Morphological partial disparities for overal mean
#' morphol.disparity(coords ~ 1, groups= ~ species*site, partial = TRUE, 
#' data = gdf, iter = 299, print.progress = FALSE)
#' 
#' # Morphological disparity without covariates, using group means
#' morphol.disparity(coords ~ species*site, groups= ~species*site, 
#' data = gdf, iter = 299, print.progress = FALSE)
#' 
#' # Morphological disparity of different groups than those 
#' # described by the linear model
#' morphol.disparity(coords ~ Csize + species*site, groups= ~ species, 
#' data = gdf, iter = 299, print.progress = FALSE)
#' 
#' # Extracting components
#' MD <- morphol.disparity(coords ~ Csize + species*site, groups= ~ species, 
#' data = gdf, iter = 299, print.progress = FALSE)
#' MD$Procrustes.var # just the Procrustes variances
#' 
#' 
#' ### Morphol.disparity can be used with previously-defined 
#' ### procD.lm or lm.rrpp class objects
#' 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land, print.progress = FALSE)    #GPA-alignment
#' gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
#' names(gp.end)<-plethspecies$phy$tip
#' 
#' gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy, 
#' gp.end = gp.end)
#' 
#' pleth.ols <- procD.lm(coords ~ Csize + gp.end, 
#' data = gdf, iter = 299) # ordinary least squares
#' pleth.pgls <- procD.pgls(coords ~ Csize + gp.end, phy = phy, 
#' data = gdf, iter = 299) # phylogenetic generalized least squares
#' 
#' summary(pleth.ols)
#' summary(pleth.pgls)
#' 
#' morphol.disparity(f1 = pleth.ols, groups = ~ gp.end, data = gdf, 
#' iter = 299, print.progress = FALSE)
#' morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end, 
#' transform = FALSE, data = gdf, 
#' iter = 299, print.progress = FALSE) # disparity in tangent space
#' morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end,
#' transform = TRUE, data = gdf, 
#' iter = 299, print.progress = FALSE) # disparity in transformed space 
#' 
#' # Three plots that correspond to the three tests
#' PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#' pPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
#' GLS = TRUE, transform = FALSE)
#' tpPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
#' GLS = TRUE, transform = TRUE)
#' 
#' par(mfrow = c(1,3))
#' 
#' # Phylomorphospace
#' PC.plot <- plot(PCA, pch = 19, phylo = TRUE, main = "PCA-OLS")
#' shapeHulls(PC.plot, groups = gp.end)
#' 
#' # Phylo-PCA
#' pPC.plot <- plot(pPCA, pch = 19, phylo = TRUE, main = "pPCA - GLS, not transformed")
#' shapeHulls(pPC.plot, groups = gp.end)
#' 
#' # Transformed phylo-PCA
#' tpPC.plot <- plot(tpPCA, pch = 19, phylo = TRUE, main = "tpPCA - GLS, transformed")
#' shapeHulls(tpPC.plot, groups = gp.end)
#' 
#' par(mfrow = c(1,1))
#' 
#'  ### Variance using RRPP (not necessarily the same as disparity)
#'  PW <- pairwise(pleth.ols, groups = gp.end)
#'  summary(PW, test.type = 'var')
#'  PW2 <- pairwise(pleth.pgls, groups = gp.end)
#'  summary(PW2, test.type = 'var')
#' 
morphol.disparity <- function(f1, groups = NULL, partial = FALSE, transform = FALSE,
                              iter = 999, seed = NULL, 
                              data = NULL, print.progress = TRUE, ...){
  
 
  if(inherits(f1, "lm")) {
    R <- as.matrix(lm$residuals)
    form <- formula(f1)
    Terms <- terms(form)
    df <- f1$model
    k <- length(attr(Terms, "term.labels"))
    if(k == 0) int.model <- TRUE else
      int.model = FALSE
  }
  
  if(inherits(f1, "lm.rrpp")) {
    if(f1$LM$gls) R <- as.matrix(f1$LM$gls.residuals) else
      R <- as.matrix(f1$LM$residuals)
    df <- f1$LM$data
    form <- f1$LM$form
    Terms <- f1$LM$Terms
    k <- length(attr(Terms, "term.labels"))
    if(k == 0) int.model <- TRUE else
      int.model = FALSE
    Pcov <- if(f1$LM$gls) f1$LM$Pcov else NULL
  }
  
  if(inherits(f1, "formula")) {
    form <- f1
    Terms <- terms(form)
    k <- length(attr(Terms, "term.labels"))
    if(k == 0) int.model <- TRUE else
      int.model = FALSE
    fit <- procD.lm(form, data = data, seed = seed, print.progress = print.progress, iter = 0, ...)
    if(fit$LM$gls) R <- as.matrix(fit$LM$gls.residuals) else
      R <- as.matrix(fit$LM$residuals)
    df <- fit$LM$data
  }
  
  if(!inherits(f1, c("lm", "lm.rrpp", "formula")))
    stop("f1 is neither a formula, lm object, nor lm.rrpp object.\n", call. = FALSE)
  
  d <- rowSums(R^2)
  N <- length(d)
  
  if(inherits(f1, "lm") || inherits(f1, "lm.rrpp")){
    
    if(is.null(groups) && !int.model) {
      if(length(attr(Terms, "term.labels") > 0)) {
        cat("\n *** Attempting to define groups from terms in the model fit.",  
            "If results are peculiar, define groups precisely and check model formula.)\n\n")
        fac.check <- which(sapply(df, is.factor))
        if(length(fac.check) == 0) groups <- NULL else {
          fac.dat <- as.matrix(df[fac.check])
          groups <- factor(apply(fac.dat, 1, function(x) paste(x, collapse = ".")))
        }
      }
    }
  }
  
  if(!is.null(groups) && inherits(groups, "formula")){
    if(length(groups) > 2) groups <- groups[-2]
    gp.var <- all.vars(groups)
    gps <- try(mget(gp.var, envir = as.environment(df)), silent = TRUE)
    if(inherits(gps, "try-error"))
      gps <- try(mget(gp.var, envir = as.environment(data)), silent = TRUE) 
    if(inherits(gps, "try-error"))
      gps <- try(mget(gp.var, envir = parent.frame()), silent = TRUE)
    if(inherits(gps, "try-error")) stop("Cannot find groups either in data frame or global environment.\n", 
                                        call. = FALSE)
    groups <- factor(apply(as.data.frame(gps), 1, function(x) paste(x, collapse = ".")))
  }
  

  if(!is.null(groups)) {
    if(!is.factor(groups) && !is.vector(groups)) 
      stop("Groups must be a factor, a vector coercible to factor, or a formula.\n", call. = FALSE)
  }

  if(!is.null(groups)) groups <- as.factor(groups)
  
  if(!int.model && partial) 
    stop("\n It is not possible to calculate partial disparities unless the model formula contains only an intercept; e.g., coords ~ 1\n")

  if(is.null(groups)) {
    pv <- sum(d)/N
    if(partial) pv <- pv * N / (N - 1)
  } else {
    if(print.progress) {
      cat("\n\nPerformimg pairwise comparisons of disparity\n")
      pb <- txtProgressBar(min = 0, max = iter+1, initial = 0, style=3)
    }
    
    newDf <- data.frame(d = d, groups = groups)
    fit <- lm.rrpp(d ~ groups + 0, iter = 0, data = newDf, print.progress = FALSE,
                   seed = seed, ...)
    Q <- fit$LM$QR
    if(transform){
      if(!is.null(Pcov)) {
        Q <- qr(Pcov %*% fit$LM$X)
      }
    }
    H <- tcrossprod(solve(qr.R(Q)), qr.Q(Q))
    ind <- perm.index(N, iter, seed)
    pv <- sapply(1:(iter + 1), function(j){
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      H %*% d[ind[[j]]]
    })
    
    if(print.progress) close(pb)
    
    if(partial){
      gpn <- as.vector(by(groups, groups, length))
      part.mat <- matrix(gpn, length(gpn), NCOL(pv))/(N - 1)
      pv <- pv * part.mat
    }
  }
  
  if(!is.null(groups)) {
    pvd <- lapply(1:(iter+1), function(j) as.matrix(dist(matrix(pv[, j])))) 
    for(i in 1:(iter + 1)) dimnames(pvd[[i]]) <- list(rownames(pv), rownames(pv))
    p.val <- Reduce("+", lapply(1:(iter+1), function(j){
      x <- pvd[[1]]
      y <- pvd[[j]]
      ifelse(y >= x, 1, 0)
    }))/(iter + 1)
  } else {
    pvd <- p.val <- permutations <- NULL
  }
  
  if(is.null(groups)) {
    pv.obs <- pv
    cat("No factor in formula or model terms from which to define groups.\n")
    cat("Procrustes variance")
    if(partial) cat("(Foote's disparity)")
    cat(":\n")
    cat(pv, "\n")
    return(pv.obs)
  } else {
    PV.dist = pvd[[1]]
    dimnames(PV.dist) <- dimnames(p.val) <- list(levels(groups), levels(groups))
    pv.obs <- pv[,1]
    names(pv.obs) <- levels(groups)
    out <- list(Procrustes.var = pv.obs, PV.dist = PV.dist, PV.dist.Pval = p.val,
                random.PV.dist = pvd, permutations = iter+1, 
                partial = partial, call = match.call())
    
    class(out) <-"morphol.disparity"
    out
    
  }
}

