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
#'   matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].
#'   The user must also specify a phylogeny describing the evolutionary relationships among species (of class phylo).
#'   Note that the specimen labels for both x and y must match the labels on the tips of the phylogeny.
#'
#'   From the phylogeny, a phylogenetic transformation matrix is obtained under a Brownian motion model, and used to 
#'   transform the x and y variables. Next, the Gower-centered distance matrix is obtained from predicted values from the
#'   model (y~x), from which sums-of-squares, F-ratios, and R^2 are estimated for each factor in the model (see Adams, 2014). 
#'   Data are then permuted across the tips of the phylogeny, and estimates of statistical values are obtained for the permuted data,
#'   which are  compared to the observed value to assess significance. 
#'   
#'   Two possible resampling procedures are provided. First, if RRPP=FALSE, 
#'   the rows of the matrix of shape variables 
#'   are randomized relative to the design matrix. This is analogous to a 'full' randomization. Second, if RRPP=TRUE,
#'   a residual randomization permutation procedure is utilized (Collyer et al. 2014). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). Effect-sizes (Z-scores) are computed as standard deviates of the sampling 
#'   distributions (of F values) generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2014).  In the case  
#'   that multiple factor or factor-covariate interactions are used in the model formula, one can specify whether all main effects should be  
#'    added to the model first, or interactions should precede subsequent main effects 
#'   (i.e., Y ~ a + b + c + a:b + ..., or Y ~ a + b + a:b + c + ..., respectively.)
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param iter Number of iterations for significance testing
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param verbose A logical value specifying whether additional output should be displayed
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS,
#' F ratio, Prand, and Rsquare.
#' @references Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional 
#' multivariate data. Evolution. 68:2675-2688. 
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
#' @examples
#' ### Example of D-PGLS for high-dimensional data 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' procD.pgls(Y.gpa$coords ~ Y.gpa$Csize,plethspecies$phy,iter=49)
#'
#' ### Example of D-PGLS for high-dimensional data, using RRPP
#' procD.pgls(Y.gpa$coords ~ Y.gpa$Csize,plethspecies$phy,iter=49,RRPP=TRUE)
procD.pgls<-function(f1, phy, iter=999, int.first = FALSE, RRPP=FALSE, verbose=FALSE){
  data=NULL
  form.in <- formula(f1)
  Y <- eval(form.in[[2]], parent.frame())
  if(length(dim(Y)) == 3)  Y <- two.d.array(Y) else Y <- as.matrix(Y)
  form.in <- as.formula(paste(c("Y",form.in[[3]]),collapse="~"))
  if(int.first == TRUE) ko = TRUE else ko = FALSE
  Terms <- terms(form.in, keep.order = ko)
  k <- length(attr(Terms, "term.labels"))
  mf <- model.frame(form.in)
  if (any(is.na(Y)) == T) {
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (is.null(dimnames(Y)[[1]])) {
    stop("No species names with Y-data")
  }
  N<-length(phy$tip.label)
  if(length(match(rownames(Y), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(Y)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  C<-vcv.phylo(phy); C<-C[rownames(Y),rownames(Y)]  
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[,1:(length(lambda))]
  Pcor <- solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect)) 
  PY <- Pcor%*%Y   
  Xs = mod.mats(form.in, mf)
  
  anova.parts.obs <- anova.pgls.parts(form.in, X=NULL, Pcor, Yalt = "observed", keep.order=ko)
  anova.tab <-anova.parts.obs$table  
  df <- anova.parts.obs$df[1:k]
  dfE <-anova.parts.obs$df[k+1]
  
  P <-array(0, c(k, 1, iter+1))
  P[,,1] <- SS.obs <- anova.parts.obs$F[1:k]
  for(i in 1:iter){
    if(RRPP == TRUE) {    
      SS.ran <- SS.pgls.random(Y, Xs, SS=SS.obs, Pcor,Yalt = "RRPP")
    } else SS.ran <- SS.pgls.random(Y, Xs, Pcor, SS=SS.obs, Yalt = "resample")
    SS.r <- SS.ran$SS
    Yr <- SS.ran$Y
    SSE.r <- SS.ran$SSE
    Fs.r <- (SS.r/df)/(SSE.r/dfE)
    P[,,i+1] <- Fs.r
  }  
  P.val <- Pval.matrix(P)
  Z <- Effect.size.matrix(P)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val, NA, NA))
  anova.title = "\nRandomization of Raw Values used\n"
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  if(verbose==TRUE)  {
    list(anova.table = anova.tab, call=match.call(), SS.rand = P)
  } else anova.tab
}
