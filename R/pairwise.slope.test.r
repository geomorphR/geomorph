#' Pairwise Comparisons of Slopes
#'
#' Function performs pairwise comparisons among slopes for groups as specified by a linear model.
#'
#' The function performs pairwise comparisons to identify differences in slopes between groups. The function is 
#' designed as a post-hoc test to MANCOVA, where the latter has identified significant shape variation explained by a 
#' covariate*group interaction term. 
#' 
#'  As input the user provides a formula describing the linear model of how shape (y) varies as a function of a factor (a) 
#'  or factorial interaction (a*b). A single covariate must also be added.  Shape data (y) can be
#'  in the form of a two-dimensional data matrix of dimension (n x [p x k]) or a 3D array (p x k x n).  
#'  It is assumed that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA)
#'  [e.g., with \code{\link{gpagen}}]. From the data, the slopes for each group are estimated, and pairwise differences 
#'  in slopes determined.
#'   
#'   It is assumed that one has verified a significant group*covariate interaction [e.g., with \code{\link{procD.lm}}].
#'   To evaluate significance of the pairwise differences, two possible resampling procedures are provided. First, if 
#'   RRPP=FALSE, the rows of the matrix of shape variables are randomized relative to the design matrix. This is 
#'   analogous to a 'full' randomization. Second, if RRPP=TRUE, a residual randomization permutation procedure 
#'   is utilized (Collyer et al. 2015). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). Effect-sizes (Z-scores) are computed as standard deviates of the sampling 
#'   distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2014).
#'   
#'   Slopes can differ in two ways: the amount of shape change per covariate unit change and the direction of shape change 
#'   associated with covariate change.  Tests statistics to compare these attributes between groups are the differences in length 
#'   and direction between slope vectors, respectively.  These statistics are calculated with the exact same random permutations used
#'   to calculate random SS for ANOVA.
#'
#'   This test is essentially the same as procD.lm with post-hoc comparisons among slopes for appropriate
#'   models.  However, differences in slopes are calculated simultaneously with the same random permutations peformed for ANOVA,
#'   making it less so a post-hoc test and more so a simultaneous test of pairwise contrasts (see Collyer et al. 2014).
#'   
#' @param f1 A formula for the linear model from which groups are to be compared (e.g., y~x1*x2)
#' @param f2 A right side formula for the covariate (e.g., ~ CS).
#' @param iter Number of iterations for permutation test
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param angle.type A value specifying whether differences between slopes should be represented by vector
#' correlations (r), radians (rad) or degrees (deg)
#' @param RRPP a logical value indicating whether residual randomization should be used for significance testing
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Computation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
#' @return Function returns a list with the following components: 
#'   \item{ANOVA.table}{An ANOVA table assessing the linear model}
#'   \item{Slope.Dist}{A matrix of pairwise differences between slope magnitudes}
#'   \item{Prob.Dist}{Associated matrix of pairwise significance levels based on permutations}
#'   \item{Slope.cor}{A matrix of pairwise slope vector correlations (if vector correlation is chosen)}
#'   \item{Prob.cor}{Associated matrix of pairwise significance levels based on permutations}
#'   \item{Slope.angle}{A matrix of pairwise angular differences in slope (if "rad" or "deg" chosen)}
#'   \item{Prob.angle}{Associated matrix of pairwise significance levels based on permutations}
#'   @examples
#' ### MANCOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' y <- Y.gpa$coords
#' 
#' ## Pairwise slope vector correlations
#' pairwise.slope.test(y~plethodon$site, ~ Y.gpa$Csize, iter=24, angle.type="r")
#' 
#' ## Pairwise angular difference between slopes
#' pairwise.slope.test(y~plethodon$site, ~ Y.gpa$Csize, iter=24, angle.type="rad")
#' 
#' ## Using RRPP
#' pairwise.slope.test(y~plethodon$site, ~Y.gpa$Csize, iter=24, angle.type="rad", RRPP=TRUE)
pairwise.slope.test <- function (f1, f2, iter = 999, int.first = FALSE, angle.type = c("r", "deg", "rad"), RRPP = FALSE){
  f1 <- as.formula(f1)
  Y <- eval(f1[[2]], parent.frame())
  if(length(dim(Y)) == 3)  Y <- two.d.array(Y) else Y <- as.matrix(Y)
  f1 <- as.formula(paste(c("Y",f1[[3]]),collapse="~"))
  fac.mf <- model.frame(f1)
  angle.type = match.arg(angle.type)
  if(int.first == TRUE) ko = TRUE else ko = FALSE
  Terms <- terms(f1, keep.order = ko)
  if (any(is.na(Y)) == T) stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  Xfacs <- model.matrix(Terms)
  newfac <- single.factor(f1, keep.order = ko)
  if(any(Xfacs != 1 & Xfacs != 0)) stop("Only factors allowed as independent variables in model formula. 
                                        \nMake sure the covariate is input separately.")
  if(is.null(f2)) stop("A covariate formula must be included, separate from the model formula")
  if(ncol(model.matrix(~newfac)) != ncol(Xfacs)) stop("Model formula must be for either a single factor or full-factorial model\n  e.g., Shape ~ Factor.A  -or-  Shape ~ Factor.A * Factor.B * ...")
  
  cov.mf <- model.frame(as.formula(f2))
  if(length(attr(terms(cov.mf),"term.labels")) > 1) stop("Pairwise comparisons of slopes can only be considered for a single covariate. \n Use advanced.procD.lm for multiple covariates or complex model designs.")
  Xs <- mod.mats.w.cov(f1=f1, f2=f2, dat1=fac.mf, dat2=cov.mf, keep.order=ko, interaction = TRUE)
  k <- length(Xs$Xs) - 1
  X <- Xs$Xs[[k+1]]
  anova.parts.obs <- anova.parts(f1, X = Xs,Yalt = "observed", keep.order=ko)
  anova.tab <-anova.parts.obs$table 
  SS.obs <- anova.parts.obs$SS[1:k]
  Bslopes <- slopes(newfac, cov.mf, Y)
  slope.lengths <- sqrt(diag(Bslopes%*%t(Bslopes)))
  db <- as.matrix(dist(Bslopes))
  cb <- vec.ang.matrix(Bslopes, type = angle.type)
  dimnames(cb)=dimnames(db)
  P <- array(0, c(k, 1, iter+1))
  P[,,1] <- SS.obs
  P.sl <- array(,c(length(slope.lengths),1,iter+1))
  P.sl[,,1] <- slope.lengths
  P.dist <- P.cor <- array(0,c(dim(db), iter+1))
  P.dist[,,1] <- db
  P.cor[,,1] <- 1 - vec.cor.matrix(Bslopes)
  
  for(i in 1: iter){
    if(RRPP == TRUE) {
      SSr <- SS.random(Y, Xs, SS.obs, Yalt = "RRPP")
      Bslopes.r <- slopes(newfac, cov.mf, SSr$Y)
    } else {
      SSr <- SS.random(Y, Xs, SS.obs, Yalt = "resample")
      Bslopes.r <- slopes(newfac, cov.mf, SSr$Y)
    }
    P[,,i+1] <- SSr$SS
    P.dist[,,i+1] <- as.matrix(dist(Bslopes.r))	
    P.cor[,, i+1] <- 1 - vec.cor.matrix(Bslopes.r)
    P.sl[,,i+1] <- sqrt(diag(Bslopes.r%*%t(Bslopes.r)))
  }
  P.val <- Pval.matrix(P)
  Z <- Effect.size.matrix(P)
  dimnames(P.dist)[1:2] = dimnames(P.cor) = dimnames(db)
  dimnames(P.sl)[[1]] <- names(slope.lengths)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val, NA, NA))
  if(RRPP == TRUE) {anova.title = "\nRandomized Residual Permutation Procedure used\n"
  } else {anova.title = "\nRandomization of Raw Values used\n"}
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  Prob.dist <- Pval.matrix(P.dist)
  Prob.cor <- Pval.matrix(P.cor)
  Prob.slope.length <- as.vector(Pval.matrix(P.sl))
  names(Prob.slope.length) <- names(slope.lengths)
  if(angle.type == "r") {
    list(anova.table=anova.tab, Slope.lengths = slope.lengths, Prob.slope.length = Prob.slope.length, 
         Slope.dist = db, Prob.dist = Prob.dist, Slope.cor = cb, Prob.cor = Prob.cor)
  } else list(anova.table=anova.tab, Slope.lengths = slope.lengths, Prob.slope.length = Prob.slope.length, Slope.dist = db, Prob.dist = Prob.dist, Slope.angle = cb, Prob.angle = Prob.cor)
}
