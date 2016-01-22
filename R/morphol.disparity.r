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
#' Absoluted differences in Procrustes variances are test statistics that can be used to test differences
#' in morphological disparity among groups.  These differences are  statistically evaluated through permutation, 
#' where the vectors of residuals are randomized among groups. The function can be used to obtain disparity for the whole
#' dataset by using "a dummy group factor "~ 1" as the right-hand portion of the formula, in which case only Procrustes 
#' variance is returned.  Additionally, if the right-hand portion of the formula only contains (continuous) covariates, 
#' e.g., "~ Csize", Procrustes variance will be calculated for the whole data set or groups, after accounting for the
#' linear regression described.  Finally, different factors can be indicated in the formula and for groups, if one wishes to 
#' compare morphological disparities for groups comprising only a portion of or collapsing of the groups in a more complex model 
#' (see examples).
#'
#' @param f1 A formula describing the linear model used.  The left-hand portion of the formula should be
#'  a matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens.  However,
#'  any matrix of multivariate data will work.  The right-hand portion of the formula should be " ~1" to use the overall mean,
#'  or "~ x1 + x2 + x3 +...", where each x is a covariate or factor.  (Interactions and nested terms also work.)
#' @param groups A formula designating groups, e.g., groups = ~ groups.  If NULL, morphol.disparity
#' will attempt to define groups based on the linear model formula, f1.  If there are no groups inherently
#' indicated in f1 and groups is NULL, a single Procrustes variance will be returned for the entire data set.
#' @param iter Number of iterations for permutation test
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}}
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
morphol.disparity <- function(f1, groups = NULL, iter = 999, seed = NULL, data = NULL, ...){
  pfit <- procD.fit(f1, data=data)
  if(!is.null(groups) & class(groups) != "formula") stop("groups must be a formula; e.g., groups = ~ X")
  if(is.null(groups)) gps <- single.factor(pfit) else {
    data.types <- lapply(data, class)
    keep = sapply(data.types, function(x) x != "array" & x != "phylo")
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
  names(pv) <- levels(gps)
  if(length(gps) == 0) {
    warning("No factor in formula from which to define groups.")
    out = (noquote(paste("Procrustes variance =",round(pv, 8))))
  } else{
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
      class(out) <-"morphol.disparity"
  }
  out
}
