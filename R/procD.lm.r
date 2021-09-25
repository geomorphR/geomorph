#' Procrustes ANOVA/regression for Procrustes shape variables
#'
#' Function performs Procrustes ANOVA with permutation procedures to assess statistical hypotheses describing 
#'   patterns of shape variation and covariation for a set of Procrustes shape variables
#'
#' The function quantifies the relative amount of shape variation attributable to one or more factors in a 
#'   linear model and estimates the probability of this variation ("significance") for a null model, via distributions generated 
#'   from resampling permutations. Data input is specified by a formula (e.g., 
#'   y~X), where 'y' specifies the response variables (Procrustes shape variables), and 'X' contains one or more independent 
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
#'   they be present (see Anderson and terBraak 2003). 
#'   
#'   Effect-sizes (Z scores) are computed as standard deviates of either the SS, 
#'   F, or Cohen's f-squared sampling distributions generated, which might be more intuitive for P-values than F-values 
#'   (see Collyer et al. 2015).  Values from these distributions are log-transformed prior to effect size estimation,
#'   to assure normally distributed data.  The SS type will influence how Cohen's f-squared values are calculated.  
#'   Cohen's f-squared values are based on partial eta-squared values that can be calculated sequentially or marginally, as with SS.
#'   
#'   In the case that multiple factor or factor-covariate interactions are used in the model 
#'   formula, one can specify whether all main effects should be added to the 
#'   model first, or interactions should precede subsequent main effects 
#'   (i.e., Y ~ a + b + c + a:b + ..., or Y ~ a + b + a:b + c + ..., respectively.)
#'
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{procD.lm}}.
#'   The generic function, \code{\link{plot}} has several options for plotting, using \code{\link{plot.procD.lm}}.  Diagnostics plots, 
#'   principal component plots (rotated to first PC of covariance matrix of fitted values), and regression plots can be performed.  
#'   One must provide a linear predictor, and
#'   can choose among predicted values (PredLine) or regression scores (RegScore). 
#'   In these plotting options, the predictor does not need to be size, and fitted values and residuals from the procD.lm fit are used rather 
#'   than mean-centered values. 
#'   
#'  \subsection{Notes for geomorph 3.1.0 and subsequent versions}{ 
#'  The procD.lm function is now a wrapper for the \code{\link{lm.rrpp}} function
#'  in the \code{RRPP} package.  Examples below illustrate how to utilize
#'  \code{RRPP} functions along with \code{geomorph} functions for procD.lm objects,
#'  increasing the breadth of possible downstream analyses.  
#'  
#'  An important update in version 3.1.0 is that advanced.procD.lm and nested.update have been deprecated.  
#'  The examples emphasize how pairwise comparisons can now be accomplished with \code{\link{pairwise}} and
#'  ANOVA updates for nested factors can be made with the \code{\link{anova.lm.rrpp}}, utilizing the error argument.
#'  These functions work on procD.lm objects that have already been created.
#' }
#' 
#'  \subsection{Notes for geomorph 3.0.6 and subsequent versions}{ 
#'  Compared to previous versions, GLS computations in random permutations are now possible in procD.lm.  One should use
#'  RRPP = TRUE if a covariance matrix is provided as an argument.  The method of SS calculations follows
#'  Adams and Collyer 2018.  Additional output with a "gls." prefix is also available.
#' }
#'   
#'  \subsection{Notes for geomorph 3.0.4 and subsequent versions}{ 
#'  Compared to previous versions of geomorph, users might notice differences in effect sizes.  Previous versions used z-scores calculated with 
#'  expected values of statistics from null hypotheses (sensu Collyer et al. 2015); however Adams and Collyer (2016) showed that expected values 
#'  for some statistics can vary with sample size and variable number, and recommended finding the expected value, empirically, as the mean from the set 
#'  of random outcomes.  Geomorph 3.0.4 and subsequent versions now center z-scores on their empirically estimated expected values and where appropriate, 
#'  log-transform values to assure statistics are normally distributed.  This can result in negative effect sizes, when statistics are smaller than 
#'  expected compared to the average random outcome.  For ANOVA-based functions, the option to choose among different statistics to measure effect size 
#'  is now a function argument.
#' }
#' 
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param iter Number of iterations for significance testing
#' @param turbo A logical value that if TRUE, suppresses coefficient estimation 
#' in every random permutation.  This will affect subsequent analyses that 
#' require random coefficients (see \code{\link{coef.lm.rrpp}})
#' but might be useful for large data sets for which only ANOVA is needed.
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing
#' @param SS.type SS.type A choice between type I (sequential), type II (hierarchical), or type III (marginal)
#' sums of squares and cross-products computations.
#' @param effect.type One of "F", "SS", or "cohen", to choose from which random distribution to estimate effect size.
#' (The option, "cohen", is for Cohen's f-squared values.  The default is "F".  Values are log-transformed before z-score calculation to
#' assure normally distributed data.)
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param Cov An optional covariance matrix that can be used for generalized least squares estimates of
#' coefficients and sums of squares and cross-products (see Adams and Collyer 2018).
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}} 
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @param Parallel Either a logical value to indicate whether parallel processing 
#' should be used or a numeric value to indicate the number of cores to use in 
#' parallel processing via the \code{parallel} library. 
#' If TRUE, this argument invokes forking of all processor cores, except one.  If
#' FALSE, only one core is used. A numeric value directs the number of cores to use,
#' but one core will always be spared.
#' @param ... Arguments not listed above that can passed on to \code{\link{lm.rrpp}}, like weights, offset.
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
#' \item{data}{The data frame for the model.}
#' \item{SS}{The sums of squares for each term, model residuals, and the total.}
#' \item{SS.type}{The type of sums of squares.  One of type I or type III.}
#' \item{df}{The degrees of freedom for each SS.}
#' \item{R2}{The coefficient of determination for each model term.}
#' \item{F}{The F values for each model term.}
#' \item{permutations}{The number of random permutations (including observed) used.}
#' \item{random.SS}{A matrix of random SS found via the resampling procedure used.}
#' \item{random.F}{A matrix or vector of random F values found via the resampling procedure used.}
#' \item{random.cohenf}{A matrix or vector of random Cohen's f-squared values
#'  found via the resampling procedure used.}
#' \item{permutations}{The number of random permutations (including observed) used.}
#' \item{effect.type}{The distribution used to estimate effect-size.}
#' \item{perm.method}{A value indicating whether "Raw" values were shuffled or "RRPP" performed.}
#' \item{gls}{This prefix will be used if a covariance matrix is provided to indicate GLS computations.}
#' @references Anderson MJ. 2001. A new method for non-parametric multivariate analysis of variance. 
#'    Austral Ecology 26: 32-46.
#' @references Anderson MJ. and C.J.F. terBraak. 2003. Permutation tests for multi-factorial analysis of variance.
#'    Journal of Statistical Computation and Simulation 73: 85-113.
#' @references Goodall, C.R. 1991. Procrustes methods in the statistical analysis of shape. Journal of the 
#'    Royal Statistical Society B 53:285-339.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C. and M.L. Collyer. 2018. Multivariate comparative methods: evaluations, comparisons, and
#' recommendations. Systematic Biology. 67:14-31.
#' @seealso \code{\link{procD.pgls}} and 
#' \code{\link{lm.rrpp}} for more on linear model fits with RRPP. See \code{\link{RRPP-package}} for further
#' details on functions that can use procD.lm objects.
#' @examples
#' \dontrun{
#' ### ANOVA example for Goodall's F test (multivariate shape vs. factors)
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)    #GPA-alignment  
#' gdf <- geomorph.data.frame(Y.gpa, 
#' site = plethodon$site, 
#' species = plethodon$species) # geomorph data frame
#'
#' fit1 <- procD.lm(coords ~ species * site, 
#' data = gdf, iter = 999, turbo = TRUE,
#' RRPP = FALSE, print.progress = FALSE) # randomize raw values
#' fit2 <- procD.lm(coords ~ species * site, 
#' data = gdf, iter = 999, turbo = TRUE,
#' RRPP = TRUE, print.progress = FALSE) # randomize residuals
#' 
#' summary(fit1)
#' summary(fit2)
#' 
#' # RRPP example applications
#' 
#' coef(fit2)
#' coef(fit2, test = TRUE)
#' anova(fit2) # same as summary
#' anova(fit2, effect.type = "Rsq")
#' # if species and site were modeled as random effects ...
#' anova(fit2, error = c("species:site", "species:site", "Residuals"))  
#' # not run, because it is a large object to print 
#' # remove # to run
#' # predict(fit2) 
#' 
#' # TPS plots for fitted values of some specimens
#' 
#' M <- Y.gpa$consensus
#' plotRefToTarget(M, fit2$GM$fitted[,,1], mag = 3)
#' plotRefToTarget(M, fit2$GM$fitted[,,20], mag = 3)
#' 
#' ### THE FOLLOWING ILLUSTRATES SIMPLER SOLUTIONS FOR 
#' ### THE NOW DEPRECATED advanced.procD.lm FUNCTION AND
#' ### PERFORM ANALYSES ALSO FOUND VIA THE morphol.disparity FUNCTION
#' ### USING THE pairwise FUNCTION
#' 
#' # Comparison of LS means, with log(Csize) as a covariate
#' 
#' # Model fits
#' fit.null <- procD.lm(coords ~ log(Csize) + species + site, data = gdf, 
#' iter = 999, print.progress = FALSE, turbo = TRUE)
#' fit.full <- procD.lm(coords ~ log(Csize) + species * site, data = gdf, 
#' iter = 999, print.progress = FALSE, turbo = TRUE)
#' 
#' # ANOVA, using anova.lm.rrpp function from the RRPP package 
#' # (replaces advanced.procD.lm)
#' anova(fit.null, fit.full, print.progress = FALSE)
#' 
#' # Pairwise tests, using pairwise function from the RRPP package
#' gp <-  interaction(gdf$species, gdf$site)
#' 
#' PW <- pairwise(fit.full, groups = gp, covariate = NULL)
#' 
#' # Pairwise distances between means, summarized two ways 
#' # (replaces advanced.procD.lm):
#' summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE)
#' summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)
#' 
#' # Pairwise comaprisons of group variances, two ways 
#' # (same as morphol.disaprity):
#' summary(PW, test.type = "var", confidence = 0.95, stat.table = TRUE)
#' summary(PW, test.type = "var", confidence = 0.95, stat.table = FALSE)
#' morphol.disparity(fit.full, groups = gp, iter = 999)
#' 
#' ### Regression example
#' data(ratland)
#' rat.gpa<-gpagen(ratland)         #GPA-alignment
#' gdf <- geomorph.data.frame(rat.gpa) # geomorph data frame is easy 
#' # without additional input
#' 
#' fit <- procD.lm(coords ~ Csize, data = gdf, iter = 999, turbo = TRUE,
#' RRPP = TRUE, print.progress = FALSE) 
#' summary(fit)
#' 
#' ### Extracting objects and plotting options
#' # diagnostic plots
#' plot(fit, type = "diagnostics") 
#' # diagnostic plots, including plotOutliers
#' plot(fit, type = "diagnostics", outliers = TRUE) 
#' 
#' # PC plot rotated to major axis of fitted values
#' plot(fit, type = "PC", pch = 19, col = "blue") 

#' # Regression-type plots 
#' 
#' # Use fitted values from the model to make prediction lines
#' plot(fit, type = "regression", 
#' predictor = gdf$Csize, reg.type = "RegScore", 
#' pch = 19, col = "green")
#' 
#' # Uses coefficients from the model to find the projected regression 
#' # scores
#' rat.plot <- plot(fit, type = "regression", 
#' predictor = gdf$Csize, reg.type = "RegScore", 
#' pch = 21, bg = "yellow") 
#' 
#' # TPS grids for min and max scores in previous plot
#' preds <- shape.predictor(fit$GM$fitted, x = rat.plot$RegScore, 
#'                         predmin = min(rat.plot$RegScore), 
#'                         predmax = max(rat.plot$RegScore))
#' M <- rat.gpa$consensus
#' plotRefToTarget(M, preds$predmin, mag=2)
#' plotRefToTarget(M, preds$predmax, mag=2)
#'                         
#' attributes(fit)
#' fit$fitted[1:3, ] # the fitted values (first three specimens)
#' fit$GM$fitted[,, 1:3] # the fitted values as Procrustes 
#' # coordinates (same specimens)
#' 
#' ### THE FOLLOWING ILLUSTRATES SIMPLER SOLUTIONS FOR 
#' ### THE NOW DEFUNCT nested.update FUNCTION USING
#' ### THE anova GENERIC FUNCTION
#' 
#' data("larvalMorph")
#' Y.gpa <- gpagen(larvalMorph$tailcoords, 
#' curves = larvalMorph$tail.sliders,
#' ProcD = TRUE, print.progress = FALSE)
#' gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
#' family = larvalMorph$family)
#' 
#' fit <- procD.lm(coords ~ treatment/family, data = gdf, turbo = TRUE,
#' print.progress = FALSE, iter = 999)
#' anova(fit) # treatment effect not adjusted
#' anova(fit, error = c("treatment:family", "Residuals")) 
#' # treatment effect updated (adjusted)
#' }
#' 
#' 
#' 
procD.lm <- function(f1, iter = 999, seed=NULL, RRPP = TRUE, 
                     SS.type = c("I", "II", "III"),
                     effect.type = c("F", "cohenf", "SS", "MS", "Rsq"),
                     int.first = FALSE,  Cov = NULL, 
                     turbo = TRUE, Parallel = FALSE,
                     data=NULL, print.progress = FALSE, ...){
  
  if(is.null(data)) {
    vars <- rownames(attr(terms(f1), "factors"))
    if(!is.null(vars)) {
      data <- list()
      data <- data[vars]
      names(data) <- vars
      for(i in 1:length(vars)) {
        f <- as.formula(paste("~", vars[i]))
        temp <- try(eval(f[[2]]), silent = TRUE)
        if(inherits(temp, "try-error")) stop("Cannot find data in global environment.\n",
                                             call. = FALSE) else
                                               data[[i]] <- temp
      }
    }
  }
  
  if(inherits(f1, "formula")){
    Y <- try(eval(f1[[2]], envir = data , enclos = parent.frame()), silent = TRUE)
    if(inherits(Y, "try-error"))
      Y <- try(eval(f1[[2]], envir = parent.frame), silent = TRUE)
    if(inherits(Y, "try-error")) stop("Cannot find data in data frame or global environment.\n",
                                      call. = FALSE)
    
    nms <- if(is.vector(Y)) names(Y) else if(inherits(Y, "matrix")) attr(Y, "Labels") else
        if(is.matrix(Y)) rownames(Y) else dimnames(Y)[[3]]
    dims.Y <- dim(Y)
    f <- update(f1, Y ~ .)
    if(length(dims.Y) == 3) {
      GM <- TRUE
      Y <- two.d.array(Y) 
      rownames(Y) <- nms
      p <- dims.Y[[1]]
      k <- dims.Y[[2]]
      n <- dims.Y[[3]]
    } else {
      GM <- FALSE
      Y <- as.matrix(Y)
      rownames(Y) <- nms
      if(isSymmetric(Y)) colnames(Y) <- nms
    }
    data$Y <- Y
    
  } else {
    f <- f1
    GM <- FALSE
  }
  
  out <- lm.rrpp(f, data = data, turbo = turbo,
                 seed = seed, RRPP = RRPP,
                 SS.type = SS.type, 
                 int.first = int.first,  
                 Cov = Cov, iter = iter,
                 print.progress = print.progress, 
                 Parallel = Parallel, ...)
  
  n <- out$LM$n
  out$ANOVA$effect.type <- match.arg(effect.type)
  out$GM <- NULL
  
  if(!out$LM$gls) {
    out$fitted <- out$LM$fitted
    out$residuals <- out$LM$residuals
    out$coefficients <- out$LM$coefficients
    if(GM) {
      out$GM$p <- p
      out$GM$k <- k
      out$GM$n <- n
      kk <- NROW(out$LM$coefficients)
      out$GM$fitted <- arrayspecs(out$LM$fitted, p, k)
      out$GM$residuals <- arrayspecs(out$LM$residuals, p, k)
      if(kk > 1) out$GM$coefficients <- arrayspecs(out$LM$coefficients, p, k) else {
        out$GM$coefficients <- array(matrix(out$LM$coefficients, p, k, byrow = TRUE), c(p,k,1))
      }
    }
  }
  
  
  if(out$LM$gls) {
    out$gls.fitted <- out$LM$gls.fitted
    out$gls.residuals <- out$LM$gls.residuals
    out$gls.coefficients <- out$LM$gls.coefficients
    out$gls.mean <- out$LM$gls.mean
    if(GM) {
      out$GM$p <- p
      out$GM$k <- k
      out$GM$n <- n
      kk <- NROW(out$LM$gls.coefficients)
      out$GM$gls.fitted <- arrayspecs(out$LM$gls.fitted, p, k)
      out$GM$gls.residuals <- arrayspecs(out$LM$gls.residuals, p, k)
      out$GM$gls.mean <- matrix(out$LM$gls.mean, out$GM$p, out$GM$k, byrow = TRUE)
      if(kk > 1) out$GM$gls.coefficients <- arrayspecs(out$LM$gls.coefficients, p, k) else {
        out$GM$coefficients <- array(matrix(out$LM$gls.coefficients, p, k, byrow = TRUE), c(p,k,1))
      }
    }
  }
  
  o.class <- class(out)
  out2 <- list()
  out2$aov.table <- anova.lm.rrpp(out)$table
  out2$call <- match.call()
  out$call <- out2$call
  if(out$LM$gls) out2$gls.coefficients <- out$LM$gls.coefficients else
    out2$coefficients <- out$LM$coefficients
  out2$Y <- out$LM$Y
  out2$X <- out$LM$X
  out2$QR <- out$LM$QR
  if(out$LM$gls) out2$gls.fitted <- out$LM$gls.fitted else
    out2$fitted <- out$LM$fitted
  if(out$LM$gls) out2$gls.residuals <- out$LM$gls.residuals else
    out2$residuals <- out$LM$residuals

  out2$weights <- if(!is.null(out$LM$weights)) out$LM$weights else NULL
  out2$Terms <- out$LM$Terms
  out2$term.labels <- out$LM$term.labels
  out2$data <- out$LM$data
  out2$random.SS <- out2$ANOVA$SS
  out <- c(out2, out)
  class(out) <- c("procD.lm", o.class)
  out
  
}
