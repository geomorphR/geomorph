#' Procrustes ANOVA/regression for shape data
#'
#' Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing 
#'   patterns of shape variation and covariation for a set of Procrustes-aligned coordinates
#'
#' The function quantifies the relative amount of shape variation attributable to one or more factors in a 
#'   linear model and estimates the probability of this variation ("significance") for a null model, via distributions generated 
#'   from resampling permutations. Data input is specified by a formula (e.g., 
#'   y~X), where 'y' specifies the response variables (shape data), and 'X' contains one or more independent 
#'   variables (discrete or continuous). The response matrix 'y' can be either in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that  -if the data based
#'   on landmark coordinates - the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) 
#'   [e.g., with \code{\link{gpagen}}]. 
#'   The names specified for the independent (x) variables in the formula represent one or more 
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the 
#'   shape matrix matches the order of values in the independent variables.  Linear model fits (using the  \code{\link{lm}} function)
#'   can also be input in place of a formula.  Arguments for \code{\link{lm}} can also be passed on via this function.
#'   
#'   The function \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates; however this step is no longer necessary, as procD.lm can receive 3D arrays as dependent variables.  It is also 
#'   recommended that \code{\link{geomorph.data.frame}} is used to create and input a data frame.  This will reduce problems caused
#'   by conflicts between the global and function environments.  In the absence of a specified data frame, procD.lm will attempt to 
#'   coerce input data into a data frame, but success is not guaranteed.
#'
#'   The function performs statistical assessment of the terms in the model using Procrustes distances among 
#'   specimens, rather than explained covariance matrices among variables. With this approach, the sum-of-squared 
#'   Procrustes distances are used as a measure of SS (see Goodall 1991). The observed SS are evaluated through 
#'   permutation. In morphometrics this approach is known as a Procrustes ANOVA (Goodall 1991), which is equivalent
#'   to distance-based anova designs (Anderson 2001). Two possible resampling procedures are provided. First, if RRPP=FALSE, 
#'   the rows of the matrix of shape variables are randomized relative to the design matrix. 
#'   This is analogous to a 'full' randomization. Second, if RRPP=TRUE, a residual randomization permutation procedure is utilized 
#'   (Collyer et al. 2015). Here, residual shape values from a reduced model are
#'   obtained, and are randomized with respect to the linear model under consideration. These are then added to 
#'   predicted values from the remaining effects to obtain pseudo-values from which SS are calculated. NOTE: for
#'   single-factor designs, the two approaches are identical.  However, when evaluating factorial models it has been
#'   shown that RRPP attains higher statistical power and thus has greater ability to identify patterns in data should
#'   they be present (see Anderson and terBraak 2003). Effect-sizes (Z-scores) are computed as standard deviates of the SS sampling 
#'   distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2015).  In the case that multiple 
#'   factor or factor-covariate interactions are used in the model formula, one can specify whether all main effects should be added to the 
#'   model first, or interactions should precede subsequent main effects 
#'   (i.e., Y ~ a + b + c + a:b + ..., or Y ~ a + b + a:b + c + ..., respectively.)
#'
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{procD.lm}}.
#'   The generic function, \code{\link{plot}}, produces diagnostic plots for Procrustes residuals of the linear fit.
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}} 
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return An object of class "procD.lm" is a list containing the following
#' \item{aov.table}{An analysis of variance table; the same as the summary.}
#' \item{call}{The matched call.}
#' \item{coefficients}{A vector or matrix of linear model coefficients.}
#' \item{Y}{The response data, in matrix form.}
#' \item{X}{The model matrix.}
#' \item{QR}{The QR decompositions of the model matrix.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{The residuals (observed responses - fitted responses).}
#' \item{weights}{The weights used in weighted least-squares fitting.  If no weights are used, 
#' NULL is returned.}
#' \item{Terms}{The results of the \code{\link{terms}} function applied to the model matrix}
#' \item{term.labels}{The terms used in constructing the aov.table.}
#' \item{SS}{The sums of squares for each term, model residuals, and the total.}
#' \item{df}{The degrees of freedom for each SS.}
#' \item{R2}{The coefficient of determination for each model term.}
#' \item{F}{The F values for each model term.}
#' \item{permutations}{The number of random permutations (including observed) used.}
#' \item{random.SS}{A matrix or vector of random SS found via the resampling procedure used.}
#' \item{perm.method}{A value indicating whether "Raw" values were shuffled or "RRPP" performed.}
#' 
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. 
#'    Austral Ecology 26: 32-46.
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Computation and Simulation 73: 85-113.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Goodall, C. R. 1991. Procrustes methods in the statistical analysis of shape. Journal of the 
#'    Royal Statistical Society B 53:285-339.
#' @seealso \code{\link{advanced.procD.lm}}, \code{\link{procD.pgls}}, and 
#' \code{\link{nested.update}} within geomorph; \code{\link[stats]{lm}} for more on linear model fits
#' @examples
#' ### MANOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)    #GPA-alignment  
#' gdf <- geomorph.data.frame(shape = Y.gpa$coords, 
#' site = plethodon$site, species = plethodon$species) # geomorph data frame
#'
#' procD.lm(shape ~ species * site, data = gdf, iter = 999, RRPP = FALSE) # randomize raw values
#' procD.lm(shape ~ species * site, data = gdf, iter = 999, RRPP = TRUE) # randomize residuals
#'
#' ### Regression example
#' data(ratland)
#' rat.gpa<-gpagen(ratland)         #GPA-alignment
#' gdf <- geomorph.data.frame(rat.gpa) # geomorph data frame is easy without additional input
#' 
#' procD.lm(coords ~ Csize, data = gdf, iter = 999, RRPP = FALSE) # randomize raw values
#' procD.lm(coords ~ Csize, data = gdf, iter = 999, RRPP = TRUE) # randomize raw values
#' # Outcomes should be exactly the same
#' 
#' ### Extracting objects and plotting residuals
#' rat.anova <- procD.lm(coords ~ Csize, data = gdf, iter = 999, RRPP = TRUE)
#' summary(rat.anova)
#' plot(rat.anova) # diagnostic plots
#' plot(rat.anova, outliers = TRUE) # diagnostic plots, including plotOutliers
#' attributes(rat.anova)
#' rat.anova$fitted # just the fitted values
procD.lm<- function(f1, iter = 999, seed=NULL, RRPP = TRUE, 
                        int.first = FALSE,  data=NULL, ...){
  if(int.first==TRUE) ko = TRUE else ko = FALSE
  pfit <- procD.fit(f1, data=data, keep.order=ko)
  k <- length(pfit$term.labels)
  if(RRPP == TRUE) P <- SS.iter(pfit,Yalt="RRPP", iter=iter, seed=seed) else 
    P <- SS.iter(pfit, Yalt="resample", iter=iter, seed=seed)
  anova.parts.obs <- anova.parts(pfit, P)
  anova.tab <-anova.parts.obs$anova.table 
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
  out <- list(aov.table = tab, call = match.call(),
              coefficients=pfit$coefficients[[k+1]], 
              Y=pfit$Y,  X=pfit$X, 
              QR = pfit$QRs[[k+1]], fitted=pfit$fitted[[k+1]],
              residuals = pfit$residuals[[k+1]], 
              weights = pfit$w, Terms = pfit$Terms, term.labels = pfit$term.labels,
              SS = anova.parts.obs$SS, df = anova.parts.obs$df, 
              R2 = anova.parts.obs$R2[1:k], F = anova.parts.obs$Fs[1:k], permutations = iter+1,
              random.SS = P, perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"))
  class(out) <- "procD.lm"
  out
}