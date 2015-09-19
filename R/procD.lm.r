#' Procrustes ANOVA/regression for shape data
#'
#' Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing 
#'   patterns of shape variation and covariation for a set of Procrustes-aligned coordinates
#'
#' The function quantifies the relative amount of shape variation attributable to one or more factors in a 
#'   linear model and estimates teh probability of this variation ("signifiance") for a null model, via distributions generated 
#'   from resampling permutations. Data input is specified by a formula (e.g., 
#'   y~X), where 'y' specifies the response variables (shape data), and 'X' contains one or more independent 
#'   variables (discrete or continuous). The response matrix 'y' can be either in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. 
#'   The names specified for the independent (x) variables in the formula represent one or more 
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the 
#'   shape matrix matches the order of values in the independent variables.  Linear model fits (using the  \code{\link{lm}} function)
#'   can also be input in place of a formula.  Arguments for \code{\link{lm}} can also be passed on via this function.
#'   
#'   The function \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates; however this step is no longer necessary, as procD.lm can receive 3D arrays as depedendent variables.
#'
#'   The function performs statistical assessment of the terms in the model using Procrustes distances among 
#'   specimens, rather than explained covariance matrices among variables. With this approach, the sum-of-squared 
#'   Procrustes distances are used as a measure of SS (see Goodall 1991). The observed SS are evaluated through 
#'   permutation. In morphometrics this approach is known as a Procrustes ANOVA (Goodall 1991), which is equivalent
#'   to distance-based anova designs (Anderson 2001). Two possible resampling procedures are provided. First, if RRPP=FALSE, 
#'   the rows of the matrix of shape variables 
#'   are randomized relative to the design matrix. This is analogous to a 'full' randomization. Second, if RRPP=TRUE,
#'   a residual randomization permutation procedure is utilized (Collyer et al. 2014). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). Effect-sizes (Z-scores) are computed as standard deviates of the SS sampling 
#'   distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2014).  In the case that multiple 
#'   factor or factor-covariate interactions are used in the model formula, one can specify whether all main effects should be added to the 
#'   model first, or interactions should precede subsequent main effects 
#'   (i.e., Y ~ a + b + c + a:b + ..., or Y ~ a + b + a:b + c + ..., respectively.)
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param iter Number of iterations for significance testing
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param verbose A logical value specifying whether additional output should be displayed
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return Function returns an ANOVA table of statistical results for all factors: df (for each factor), SS, MS,
#' Rsquare, F ratio, Z-score, and Prand.  If verbose=TRUE, random SS are provided.
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. 
#'    Austral Ecology 26: 32-46.
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Copmutation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357–365.
#' @references Goodall, C. R. 1991. Procrustes methods in the statistical analysis of shape. Journal of the 
#'    Royal Statistical Society B 53:285-339.
#' @examples
#' ### MANOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#'
#' procD.lm(Y.gpa$coords ~ plethodon$species*plethodon$site,iter=499)
#'
#' ### Regression example
#' data(ratland)
#' rat.gpa<-gpagen(ratland)         #GPA-alignment
#'
#' procD.lm(rat.gpa$coords ~ rat.gpa$Csize,iter=499)
#' 
#' ## using RRPP
#'  procD.lm(rat.gpa$coords ~ rat.gpa$Csize,iter=49,RRPP=TRUE)
procD.lm<- function(f1, iter = 999, RRPP = FALSE, int.first = FALSE, verbose=FALSE, ...){
  if(int.first==TRUE) ko = TRUE else ko = FALSE
  dat <- as.data.frame(model.frame(f1[-2]))
  if(any(class(f1)=="lm")) pf = procD.fit(f1,weights=f1$weights, contrasts=f1$contrasts, offset=f1$offset, data=dat) else 
    pf= procD.fit(f1, data=dat,...)
  anova.parts.obs <- anova.parts(pf, keep.order=ko)
  anova.tab <-anova.parts.obs$table  
  Xs <- pf$Xs
  if(RRPP == TRUE) P <- SS.random(pf,Yalt="RRPP", iter=iter) else P <- SS.random(pf, Yalt="resample", iter=iter)
  P.val <- Pval.matrix(P)
  Z <- Effect.size.matrix(P)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val, NA, NA))
  if(RRPP == TRUE) {
    anova.title = "\nRandomized Residual Permutation Procedure used\n"
  } else anova.title = "\nRandomization of Raw Values used\n"
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  if(verbose==TRUE)  {
    list(anova.table = anova.tab, call=match.call(), SS.rand = P)
  } else anova.tab
}
