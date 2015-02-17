#' Pairwise Group Comparisons
#'
#' Function performs pairwise comparisons among groups using the Euclidean distances among group means.
#'
#' The function performs pairwise comparisons to identify shape differences among groups. The function is designed as a post-hoc
#'  test to Procrustes ANOVA, where the latter has identified significant shape variation explained by a grouping factor. 
#'  
#'  The Function can handle, single-factor ANOVA, factorial ANOVA, and ANCOVA designs.
#'  
#'  As input the user provides a formula describing the linear model of how shape (y) varies as a function of a factor (a) 
#'  or factorial interaction (a*b). A single covariate, matrix of covariates, or data frame of covariates can also be added.
#'  E.g., covariates = x1, covariates = cbind(x1, x2, x3,...), or covariates = data.frame(x1, x2, x3,...).
#'  The shape data (y) must be in the form of a two-dimensional data matrix of dimension (n x [p x k]), rather than a 3D array.  
#'  It is assumed that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) 
#'  [e.g., with \code{\link{gpagen}}]. The function \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix 
#'  from a 3D array of landmark coordinates. From the data, the Euclidean distances among group means are estimated, and used 
#'  as test values.
#'   
#'   To evaluate significance of group differences, two possible resampling procedures are provided. First, if 
#'   RRPP=FALSE, the rows of the matrix of shape variables are randomized relative to the design matrix. This is 
#'   analogous to a 'full' randomization. Second, if RRPP=TRUE, a residual randomization permutation procedure 
#'   is utilized (Collyer et al. 2014). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). Effect-sizes (Z-scores) are computed as standard deviates of the SS 
#'   sampling distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2014).
#'   
#'   This test is essentially the same as procD.lm with post-hoc comparisons among least squares (LS) means for appropriate
#'   models.  However, differences in means are calculated simultaneously with the same random permutations peformed for ANOVA,
#'   making it less so a post-hoc test and more so a simultaneous test of pairwise contrasts (see Collyer et al. 2014).
#'
#' @param f1 A formula for the linear model from which groups are to be compared (e.g., y~x1*x2)
#' @param f2 A right-side formula for one or more covariates (e.g., ~ CS + altitude)
#' @param iter Number of iterations for permutation test
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @keywords analysis
#' @export
#' @author Michael Collyer and Dean Adams
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Computation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
#' @return Function returns a list with the following components: 
#'   \item{Obs.dist}{A matrix of Euclidean distances among group means or Least Squares group means}
#'   \item{Prob.Dist}{A matrix of pairwise significance levels based on permutation}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y<-two.d.array(Y.gpa$coords)
#' ### Procrustes ANOVA
#' procD.lm(y~plethodon$species,iter=99)
#' 
#' ### Pairwise comparisons: ANOVA design with full randomization
#' pairwiseD.test(y~plethodon$species*plethodon$site,iter=9)
#' 
#' ### Pairwise comparisons: ANOVA design with residual randomization
#' pairwiseD.test(y~plethodon$species*plethodon$site,iter=9,RRPP=TRUE)
#' 
#' ### Pairwise comparisons: ANCOVA design with full randomization
#' pairwiseD.test(y~plethodon$species*plethodon$site, ~ Y.gpa$Csize, iter=9)
#' 
#' ### Pairwise comparisons: ANCOVA design with residual randomization
#' pairwiseD.test(y~plethodon$species*plethodon$site, ~ Y.gpa$Csize, iter=9, RRPP = TRUE)
#' 
pairwiseD.test <- function(f1, f2 = NULL, RRPP = FALSE, int.first = FALSE, iter= 999){
  f1 <- as.formula(f1)
  if(int.first == TRUE) ko = TRUE else ko = FALSE
  fTerms <- terms(as.formula(f1), keep.order = ko)
  fac.mf <- model.frame(f1)
  Y <- as.matrix(fac.mf[1])
  if (length(dim(Y)) != 2) {
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")
  }
  if (any(is.na(Y)) == T) {
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  }
  Xfacs <- model.matrix(fTerms)
  newfac <- single.factor(f1, keep.order = ko)
  if(is.null(f2)){
    cov.mf = NULL
  } else {cov.mf = model.frame(as.formula(f2))}
  if(any(Xfacs != 1 & Xfacs != 0)) stop("Only factors allowed as independent variables in model formula. 
                                        \nMake sure covariates are input separately. \ne.g., f2 = ~ X1 + X2 +...")
  if(ncol(model.matrix(~newfac)) != ncol(Xfacs)) stop("Model formula must be for either a single factor or full-factorial model\n  e.g., Shape ~ Factor.A  -or-  Shape ~ Factor.A * Factor.B * ...")
  if(is.null(f2)) {Xs <-mod.mats(fac.mf, keep.order=ko) 
  } else {Xs <- mod.mats.w.cov(fac.mf, cov.mf, keep.order=ko)}
  k <- length(Xs$Xs) - 1
  X <- Xs$Xs[[k+1]]
  anova.parts.obs <- anova.parts(f1, X = Xs,Yalt = "observed", keep.order=ko)
  anova.tab <-anova.parts.obs$table
  SS.obs <- anova.parts.obs$SS[1:k]
  ls.means.obs <- ls.means(newfac, cov.mf, Y)
  dm <- as.matrix(dist(ls.means.obs))
  P <- array(0, c(k, 1, iter+1))
  P[,,1] <- SS.obs
  P.dist <- array(0,c(dim(dm), iter+1))
  P.dist[,,1] <- dm
  
  for(i in 1:iter){
    if(RRPP == TRUE) {
      SSr <- SS.random(Y, Xs, SS.obs, Yalt = "RRPP")
    } else {
      SSr <- SS.random(Y, Xs, SS.obs, Yalt = "resample")
    }
    P[,,i+1] <- SSr$SS
    P.dist[,,i+1] <- as.matrix(dist(ls.means(newfac, cov.mf, SSr$Y)))
  }
  P.val <- Pval.matrix(P)
  Z <- Effect.size.matrix(P)
  dimnames(P.dist)[1:2] = dimnames(dm)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA,NA), P.value = c(P.val, NA, NA))
  if(RRPP == TRUE) {anova.title = "\nRandomized Residual Permutation Procedure used\n"
  } else {anova.title = "\nRandomization of Raw Values used\n"}
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  Prob.dist <- Pval.matrix(P.dist)
  list(anova.table=anova.tab, Obs.dist = dm, Prob.dist = Prob.dist)	
}