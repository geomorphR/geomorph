#' Phylogenetic ANOVA/regression for shape data
#'
#' Function performs Procrustes ANOVA in a phylogenetic framework and uses  permutation procedures to assess 
#' statistical hypotheses describing patterns of shape variation and covariation for a set of Procrustes-aligned coordinates
#'
#' The function performs ANOVA and regression models in a phylogenetic context under a Brownian motion model of evolution, 
#' in a manner that can accommodate 
#' high-dimensional datasets. The approach is derived from the statistical equivalency between parametric methods 
#' utilizing covariance matrices and methods based on distance matrices (Adams 2014). Data input is specified by 
#' a formula (e.g., y~X), where 'y' specifies the response variables (shape data), and 'X' contains one or more 
#' independent variables (discrete or continuous). The response matrix 'y' can be either in the form of a two-dimensional data 
#' matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that the landmarks have previously 
#' been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].
#' Linear model fits (using the  \code{\link{lm}} function)
#' can also be input in place of a formula.  Arguments for \code{\link{lm}} can also be passed on via this function.
#' The user must also specify a phylogeny describing the evolutionary relationships among species (of class phylo).
#' Note that the specimen labels for both X and y must match the labels on the tips of the phylogeny.
#'
#'   The function \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates; however this step is no longer necessary, as procD.lm can receive 3D arrays as dependent variables.  It is also 
#'   recommended that \code{\link{geomorph.data.frame}} is used to create and input a data frame.  This will reduce problems caused
#'   by conflicts between the global and function environments.  In the absence of a specified data frame, procD.pgls will attempt to 
#'   coerce input data into a data frame, but success is not guaranteed.
#'   
#'   From the phylogeny, a phylogenetic transformation matrix is obtained under a Brownian motion model, and used to 
#'   transform the X and y variables. Next, the Gower-centered distance matrix is obtained from predicted values from the
#'   model (y~X), from which sums-of-squares, F-ratios, and R^2 are estimated for each factor in the model (see Adams, 2014). 
#'   Data are then permuted across the tips of the phylogeny, and all estimates of statistical values are obtained for the permuted data,
#'   which are compared to the observed value to assess significance. This approach has been shown to have appropriate type I error
#'   rates, whereas an alternative procedure for phylogenetic regression of morphometric shape data displays elevated type I error rates
#'   (see Adams and Collyer 2015). 
#'   
#'   Two possible resampling procedures are provided. First, if RRPP=FALSE, 
#'   the rows of the matrix of shape variables 
#'   are randomized relative to the design matrix. This is analogous to a 'full' randomization. Second, if RRPP=TRUE,
#'   a residual randomization permutation procedure is utilized (Collyer et al. 2015). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). Effect-sizes (Z-scores) are computed as standard deviates of the sampling 
#'   distributions (of F values) generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2015).  In the case  
#'   that multiple factor or factor-covariate interactions are used in the model formula, one can specify whether all main effects should be  
#'    added to the model first, or interactions should precede subsequent main effects 
#'   (i.e., Y ~ a + b + c + a:b + ..., or Y ~ a + b + a:b + c + ..., respectively.)
#'
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{procD.pgls}}.
#'   The generic function, \code{\link{plot}}, produces diagnostic plots for Procrustes residuals of the linear fit.
#'   
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}} 
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return procD.lm.pgls returns an object of class "procD.lm".  
#' See \code{\link{procD.lm}} for a description of the list of results generated.  Additionally, procD.pgls provides
#' the phylogenetic correction matrix, Pcor.
#' @references Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional 
#' multivariate data. Evolution. 68:2675-2688. 
#' @references Adams, D.C., and M.L. Collyer. 2015. Permutation tests for phylogenetic comparative analyses of high-dimensional 
#' shape data: what you shuffle matters. Evolution. 69:823-829.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @examples
#' ### Example of D-PGLS for high-dimensional data 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' gdf <- geomorph.data.frame(Y.gpa, phy = plethspecies$phy)
#' procD.pgls(coords ~ Csize, phy = phy, data = gdf, iter = 999, RRPP = FALSE) # randomize raw values
#' procD.pgls(coords ~ Csize, phy = phy, data = gdf, iter = 999, RRPP = TRUE) # randomize residuals
#' 
#' ### Extracting objects
#' pleth.pgls <- procD.pgls(coords ~ Csize, phy = phy, data = gdf, iter = 999, RRPP = TRUE)
#' summary(pleth.pgls)
#' plot(pleth.pgls)
#' pleth.pgls$Pcor # the phylogenetic transformation (correction) matrix
procD.pgls<-function(f1, phy, iter=999, seed=NULL, int.first = FALSE, 
                         RRPP=TRUE, data=NULL, ...){
  if(int.first==TRUE) ko = TRUE else ko = FALSE
  pfit <- procD.fit(f1, data=data, keep.order=ko)
  Terms <- pfit$Terms
  k <- length(pfit$term.labels) 
  Y <- as.matrix(pfit$wY)
  phy.name <- deparse(substitute(phy))
  phy.match <- match(phy.name, names(data))
  if(length(phy.match) > 1) stop("More than one object of class phylo in data frame")
  if(all(is.na(phy.match))) phy <- phy else phy <- data[[phy.match]]
  N<-length(phy$tip.label)
  if(length(match(rownames(Y), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(Y)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  C <- vcv.phylo(phy)
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[,1:(length(lambda))]
  Pcor <- fast.solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect)) 
  dimnames(Pcor) <- dimnames(C)
  if(RRPP == TRUE) SSr <- Fpgls.iter(pfit, Yalt="RRPP", Pcor, iter=iter, seed=seed) else 
    SSr <- Fpgls.iter(pfit, Yalt="resample", Pcor, iter=iter, seed=seed)
  anova.parts.obs <- anova.parts.pgls(pfit, SSr)
  anova.tab <-anova.parts.obs$anova.table 
  P <- SSr$Fs
  if(is.matrix(P)){
    P.val <- apply(P,1,pval)
    Z <- apply(P,1,effect.size) 
    rownames(P) <- pfit$term.labels
    colnames(P) <- c("obs", paste("iter", 1:iter, sep=":"))
  } else {
    P.val <- pval(P)
    Z <- effect.size(P) 
    names(P) <-c("obs", paste("iter", 1:iter, sep=":"))
  }
  tab <- data.frame(anova.tab, Z = c(Z, NA, NA), Pr = c(P.val, NA, NA))
  colnames(tab)[1] <- "Df"
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  class(tab) <- c("anova", class(tab))
  pfit <- procD.fit(f1, data=data, keep.order=ko, pca=FALSE)
  out = list(aov.table = tab, call = match.call(),
             coefficients=pfit$coefficients, 
             Y=pfit$Y,  X=pfit$X, 
             Pcor=Pcor, 
             QR = pfit$QRs[[k+1]], fitted=pfit$fitted[[k+1]], 
             residuals = pfit$residuals[[k+1]], 
             weights = pfit$w, Terms = pfit$Terms, term.labels = pfit$term.labels,
             SS = anova.parts.obs$SS, df = anova.parts.obs$df, R2 = anova.parts.obs$R2[1:k], 
             F = anova.parts.obs$Fs[1:k], permutations = iter+1,
             random.SS = P, perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"))
  class(out) <- "procD.lm"
  out
}
