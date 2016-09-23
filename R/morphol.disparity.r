#' Morphological disparity for one or more groups of specimens 
#'
#' Function estimates morphological disparity and performs pairwise comparisons among groups. 
#'
#' The function estimates morphological disparity and performs pairwise comparisons to identify differences
#' between groups. Morphological disparity is estimated as the Procrustes variance, overall or for groups, 
#' using residuals of a linear model fit.  Procrustes variance is the same sum of the diagonal elements 
#' of the group covariance matrix divided by the number of observations in the group (e.g., Zelditch et al. 2012).
#' The function takes as input a formula to describe the linear model fit,
#' plus a formulaic indication of groups (e.g., ~ groups).  It is assumed that the formula describes shape data that 
#' have been GPA-aligned [e.g., \code{\link{gpagen}}], although the function can work with any multivariate data.
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
#' This function can be used with an object of class "procD.lm" or "advanced.procD.lm", if such analyses have already been performed.  This is especially helpful for analyses performed with \code{\link{gpagen}},
#' This is specially useful for analyses performed with  \code{\link{procD.pgls}}.  In this case, residuals obtained from PGLS estimation
#' of coefficients, rather than OLS estimation, will be used in the analysis.  Thus, one can account for phylogeny when comparing
#' morphological disparity among groups.  However, one should be aware that this approach only adjusts expected values because of phylogeny
#' and does not assert a null hypothesis of equal variances based on phylogenetic relatedness.  (For example, specials means can be adjusted
#' because of phylogenetic relatedness, but the null hypothesis of equal variances is conditioned on the estimation of means.) 
#'
#' @param f1 A formula describing the linear model used.  The left-hand portion of the formula should be
#'  a 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens, or a matrix (n x variables). 
#'  The right-hand portion of the formula should be " ~1" to use the overall mean, or "~ x1 + x2 + x3 +...", where each x is a 
#'  covariate or factor.  (Interactions and nested terms also work.)  Alternatively, one can use an object of class "procD.lm" or
#'  "advanced.procD.lm", which already has a formula defined.  This is especially helpful for analyses performed with \code{\link{procD.pgls}},
#'  as residuals from PGLS will be used for analysis (see examples).
#' @param groups A formula designating groups, e.g., groups = ~ groups.  If NULL, morphol.disparity
#' will attempt to define groups based on the linear model formula, f1.  If there are no groups inherently
#' indicated in f1 and groups is NULL, a single Procrustes variance will be returned for the entire data set.
#' @param iter Number of iterations for permutation test
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}}
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @keywords analysis
#' @export
#' @author Emma Sherratt and Michael Collyer
#' @return Objects of class "morphol.disparity" return a list with the following components 
#' (if groups are specified): 
#'   \item{Procrustes.var}{Observed Procrustes variances.}
#'   \item{PV.dist}{Observed pairwise absolute differences (distances) 
#'   among group Procrustes variances.}
#'   \item{PV.dist.Pval}{P-values associated with pairwise differences.}
#'   \item{random.PV.dist}{Pairwise distance matrices produced in the resampling procedure.}
#'   \item{permutations}{Number of random permutations in resampling procedure.}
#'   \item{call}{The match call}
#'   
#' @references Zelditch, M. L., D. L. Swiderski, H. D. Sheets, and W. L. Fink. 2012. Geometric morphometrics 
#'   for biologists: a primer. 2nd edition. Elsevier/Academic Press, Amsterdam.
#' @examples
#' data(plethodon)
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)
#' 
#' # Morphological disparity for entire data set
#' morphol.disparity(coords ~ 1, groups= NULL, data = gdf, iter=999)
#' 
#' # Morphological disparity for entire data set, accounting for allometry
#' morphol.disparity(coords ~ Csize, groups= NULL, data = gdf, iter=999)
#' 
#' # Morphological disparity without covariates, using overall mean
#' morphol.disparity(coords ~ 1, groups= ~ species*site, data = gdf, iter=999)
#' 
#' # Morphological disparity without covariates, using group means
#' morphol.disparity(coords ~ species*site, groups= ~species*site, data = gdf, iter=999)
#' 
#' # Morphological disparity of different groups than those described by the linear model
#' morphol.disparity(coords ~ Csize + species*site, groups= ~ species, data = gdf, iter=999)
#' 
#' # Extracting components
#' MD <- morphol.disparity(coords ~ Csize + species*site, groups= ~ species, data = gdf, iter=999)
#' MD$Procrustes.var # just the Procrustes variances
#' 
#' 
#' ### Morphol.disparity can be used with procD.lm or advanced.procd.lm class objects
#' 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
#' names(gp.end)<-plethspecies$phy$tip
#' 
#' gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy, gp.end = gp.end)
#' 
#' pleth.ols <- procD.lm(coords ~ Csize + gp.end, 
#' data = gdf, iter = 999) # ordinary least squares
#' pleth.pgls <- procD.pgls(coords ~ Csize + gp.end, phy = phy, 
#' data = gdf, iter = 999) # phylogenetic generalized least squares
#' 
#' summary(pleth.ols)
#' summary(pleth.pgls)
#' 
#' morphol.disparity(f1 = pleth.ols, groups = ~ gp.end, data = gdf, iter = 999)
#' morphol.disparity(f1 = pleth.pgls, groups = ~ gp.end, data = gdf, iter = 999)
#' 
morphol.disparity <- function(f1, groups = NULL, iter = 999, seed = NULL, 
                              data = NULL, print.progress = TRUE, ...){
  if(!is.null(data)) data <- droplevels(data)
  if(!is.null(groups) & class(groups) != "formula") stop("groups must be a formula; e.g., groups = ~ X")
  if(class(f1) == "formula") {
    pfit <- procD.fit(f1, data=data)
    if(is.null(groups)) gps <- single.factor(pfit) else {
      data.types <- lapply(data, class)
      keep = sapply(data.types, function(x) x != "array" & x != "phylo" & x != "dist")
      dat2 <- as.data.frame(data[keep])
      gps <- model.frame(groups, data= dat2)
      if(ncol(gps) > 1) gps <- factor(apply(gps, 1,function(x) paste(x, collapse=":"))) else 
        gps <- as.factor(unlist(gps))
    }
    k <- length(pfit$wResiduals)
    R <- as.matrix(pfit$wResiduals[[k]])
    if(length(gps) == 0) pv = sum(R^2)/nrow(R) else 
      pv = sapply(1:nlevels(gps), function(j){
        x <- R[gps==levels(gps)[j],]
        sum(x^2)/nrow(x)
      })
  }
  
  if(class(f1) == "procD.lm" || class(f1) == "advanced.procD.lm"){
    if(is.null(f1$pgls.residuals)) R <- as.matrix(f1$residuals) else 
      R <- as.matrix(f1$pgls.residuals)
    k <- length(f1$term.labels) + 1
    if(is.null(groups) || is.null(data)){
      X <- f1$X
      dummy <- function(x){
        check <- x == 0 | x == 1 | x == -1
        ifelse(all(check), TRUE, FALSE)
      }
      keep <- apply(X,2,dummy)
      X <- X[,keep]
      if(NCOL(X) > 1 && sum(X[,1]) == NROW(X)) X <- X[,-1]
      unX <- as.matrix(unique(X)); unX <- unX[order(unX[,1]),]
      gps <- as.factor(match(as.data.frame(t(X)), as.data.frame(t(unX))))
      cat("\n\n
          Because either no groups were defined or no geomorph data frame was provided, 
          an attempt was made to define groups from the model design matrix.  Group levels are thus 
          defined by design matrix levels in numrical order.  To be precise, rerun the anlaysis with 
          groups defined and a geomorph data frame provided (see example in help file).
          \n\n\n.")
      
      if(length(levels(gps)) == 1) pv = sum(R^2)/nrow(R) else 
        pv = sapply(1:nlevels(gps), function(j){
          x <- R[gps==levels(gps)[j],]
          sum(x^2)/nrow(x)
        })
    } else {
      data.types <- lapply(data, class)
      keep = sapply(data.types, function(x) x != "array" & x != "phylo" & x != "dist")
      dat2 <- as.data.frame(data[keep])
      gps <- model.frame(groups, data= dat2)
      if(ncol(gps) > 1) gps <- factor(apply(gps, 1,function(x) paste(x, collapse=":"))) else 
        gps <- as.factor(unlist(gps))
      if(length(gps) == 0) pv = sum(R^2)/nrow(R) else 
        pv = sapply(1:nlevels(gps), function(j){
          x <- R[gps==levels(gps)[j],]
          sum(x^2)/nrow(x)
        })
    }
  }
  
  names(pv) <- levels(gps)
  if(length(gps) == 0) {
    cat("No factor in formula from which to define groups.\n")
    cat("Procrustes variance:\n")
    out <- as.numeric(pv)
  } else{
    if(print.progress){
      ind <- perm.index(nrow(R),iter, seed=seed)
      pb <- txtProgressBar(min = 0, max = iter, initial = 0, style=3) 
      P <- lapply(1:(iter+1), function(j){
        r <- R[ind[[j]],]
        setTxtProgressBar(pb,j)
        pvr <- sapply(1:nlevels(gps), function(j){
          x <- r[gps==levels(gps)[j],]
          sum(x^2)/nrow(x)
        })
        names(pvr)<-levels(gps)
        as.matrix(dist(pvr))
      } )
      close(pb)
      if(iter > 0) names(P) <- c("obs", 1:iter)
      p.val <- Pval.matrix(simplify2array(P))
      pvd <- P[[1]]
      out <- list(Procrustes.var = pv, PV.dist = pvd, PV.dist.Pval = p.val,
                  random.PV.dist = P, permutations = iter+1, call = match.call())
    } else {
      ind <- perm.index(nrow(R),iter, seed=seed)
      P <- lapply(1:(iter+1), function(j){
        r <- R[ind[[j]],]
        pvr <- sapply(1:nlevels(gps), function(j){
          x <- r[gps==levels(gps)[j],]
          sum(x^2)/nrow(x)
        })
        names(pvr)<-levels(gps)
        as.matrix(dist(pvr))
      } )
      if(iter > 0) names(P) <- c("obs", 1:iter)
      p.val <- Pval.matrix(simplify2array(P))
      pvd <- P[[1]]
      out <- list(Procrustes.var = pv, PV.dist = pvd, PV.dist.Pval = p.val,
                  random.PV.dist = P, permutations = iter+1, call = match.call())
    }
      
      class(out) <-"morphol.disparity"
  }
  out
}
