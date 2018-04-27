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
#'   principal component plots (rotated to first PC of covariance matrix of fitted values), and regression plots can be performed.  The
#'   latter is fundamentally similar to the plotting options for \code{\link{procD.allometry}}.  One must provide a linear predictor, and
#'   can choose among common regression component (CRC), predicted values (PredLine), or regression scores (RegScore).  See \code{\link{procD.allometry}} 
#'   for details. In these plotting options, the predictor does not need to be size, and fitted values and residuals from the procD.lm fit are used rather 
#'   than mean-centered values. 
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
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing
#' @param effect.type One of "F", "SS", or "cohen", to choose from which random distribution to estimate effect size.
#' (The option, "cohen", is for Cohen's f-squared values.  The default is "F".  Values are log-transformed before z-score calculation to
#' assure normally distributed data.)
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param Cov An optional covariance matrix that can be used for generalized least squares estimates of
#' coefficients and sums of squares and cross-products (see Adams and Collyer 2018).
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}} 
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function,
#' such as weights or offset).  The function procD.fit can also currently
#' handle either type I, type II, or type III sums of squares and cross-products (SSCP) calculations.  Choice of SSCP type can be made with the argument,
#' SS.type; i.e., SS.type = "I" or SS.type = "III".  Only advanced users should consider using these additional arguments, as such arguments
#' are experimental in nature. 
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
#' @references Adams, D.C. and M.L. Collyer. 2018. Multivariate phylogenetic comparative methods: evaluations, comparisons, and
#' recommendations. Systematic Biology. 67:14-31.
#' @seealso \code{\link{advanced.procD.lm}}, \code{\link{procD.pgls}}, and 
#' \code{\link{nested.update}} within geomorph; \code{\link[stats]{lm}} for more on linear model fits.
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
#' ### Extracting objects and plotting options
#' rat.anova <- procD.lm(coords ~ Csize, data = gdf, iter = 999, RRPP = TRUE)
#' summary(rat.anova)
#' # diagnostic plots
#' plot(rat.anova, type = "diagnostics") 
#' # diagnostic plots, including plotOutliers
#' plot(rat.anova, type = "diagnostics", outliers = TRUE) 
#' # PC plot rotated to major axis of fitted values
#' plot(rat.anova, type = "PC", pch = 19, col = "blue") 
#' # Uses residuals from model to find the commonom regression component 
#' # for a predictor from the model
#' plot(rat.anova, type = "regression", predictor = gdf$Csize, reg.type = "CRC", 
#' pch = 19, col = "green")
#' # Uses residuals from model to find the projected regression scores
#' rat.plot <- plot(rat.anova, type = "regression", predictor = gdf$Csize, reg.type = "RegScore", 
#' pch = 21, bg = "yellow") 
#' 
#' # TPS grids for min and max scores in previous plot
#' preds <- shape.predictor(gdf$coords, x = rat.plot$RegScore, 
#'                         predmin = min(rat.plot$RegScore), 
#'                         predmax = max(rat.plot$RegScore))
#' M <- rat.gpa$consensus
#' plotRefToTarget(M, preds$predmin, mag=3)
#' plotRefToTarget(M, preds$predmax, mag=3)
#'                         
#' attributes(rat.anova)
#' rat.anova$fitted # just the fitted values
procD.lm<- function(f1, iter = 999, seed=NULL, RRPP = TRUE, effect.type = c("F", "SS", "cohen"),
                    int.first = FALSE,  Cov = NULL, data=NULL, print.progress = TRUE, ...){
  if(int.first) ko = TRUE else ko = FALSE
  if(!is.null(data)) data <- droplevels(data)
  pfit <- procD.fit(f1, data=data, keep.order=ko,  pca=FALSE, ...)
  Y <- pfit$Y
  n <- dim(pfit$Y)[[1]]
  p <- dim(pfit$Y)[[2]]
  k <- length(pfit$term.labels)
  if(p > n) pfitr <- procD.fit(f1, data=data, keep.order=ko,  pca=TRUE, ...) else
    pfitr <- pfit
  id <- rownames(Y)
  ind <- perm.index(n, iter=iter, seed = seed)
  SS.args <- list(fit = pfitr, ind = ind, P = NULL,
                  RRPP = RRPP, print.progress = print.progress)
  if(!is.null(Cov)){
    Cov.name <- deparse(substitute(Cov))
    Cov.match <- match(Cov.name, names(data))
    if(length(Cov.match) > 1) stop("More than one object matches covariance matrix name")
    if(all(is.na(Cov.match))) Cov <- Cov else Cov <- data[[Cov.match]]
    if(!is.matrix(Cov)) stop("The covariance matrix must be a matrix.")
    dimsC <- dim(Cov)
    if(!all(dimsC == n))
      stop("Either one or both of the dimensions of the covariance matrix do not match the number of observations.")
    if(is.null(id) || is.null(rownames(Cov))) 
      cat("\nWarning: No names to match between data and covariance matrix; consistent order is assumed.\n")
    Pcov <- Cov.proj(Cov, id)
    SS.args$P <- Pcov
  } else {
    Pcov <- NULL
  }
  if(k > 0) {
    SS <- do.call(SS.iter, SS.args)
    anova.parts.obs <- anova.parts(pfitr, SS) 
    anova.tab <-anova.parts.obs$anova.table 
    df <- anova.parts.obs$df
    effect.type <- match.arg(effect.type)
    if(effect.type == "SS" && !is.null(Cov)) {
      cat("\nWarning: Measuring effect size on SS with GLS does not make sense; F used instead.\n")
      effect.type = "F"
    }
    SS.type <- pfit$SS.type
    SSE <- SS[k + 1, ]
    SSY <- SS[k + 2, ]
    MSE <- SSE/df[k + 1]
    SSE.mat <- matrix(SSE, k, length(SSE), byrow = TRUE)
    MSE.mat <- matrix(MSE, k, length(MSE), byrow = TRUE)
    Fs <- (SS[1:k,]/df[1:k])/MSE.mat
    if(SS.type == "III") {
      if(k == 1) etas <- SS/(SS + SSE.mat) else etas <- SS[1:k,]/(SS[1:k,]+SSE.mat)
      cohenf <- etas/(1-etas)
    } else {
      if(k == 1) etas <- SS/SSY else etas <- SS[1:k,]/SSY
      unexp <- 1 - apply(etas, 2, cumsum)
      cohenf <- etas/unexp
    }
    if(k == 1) P.val <- pval(Fs) else P.val <- apply(Fs, 1, pval)
    if(effect.type == "SS") Z <- apply(log(SS[1:k,]), 1, effect.size) else
      if(effect.type == "F") Z <- apply(log(Fs), 1, effect.size) else
        Z <- apply(log(cohenf), 1, effect.size) 
    colnames(SS) <- colnames(Fs) <- colnames(cohenf) <- c("obs", paste("iter", 1:iter, sep=":"))
    tab <- data.frame(anova.tab, Z = c(Z, NA, NA), Pr = c(P.val, NA, NA))
    colnames(tab)[1] <- "Df"
    colnames(tab)[ncol(tab)] <- "Pr(>F)"
    if(effect.type == "SS") colnames(tab)[ncol(tab)] <- "Pr(>SS)"
    if(effect.type == "cohen") colnames(tab)[ncol(tab)] <- "Pr(>Cohen f-sq)"
    class(tab) <- c("anova", class(tab))
    rownames(SS) <- c(pfit$term.labels, "Residuals", "Total")
    out <- list(aov.table = tab, call = match.call(),
                coefficients=pfit$wCoefficients.full[[k]], 
                Y=pfit$Y,  X=pfit$X, 
                QR = pfit$wQRs.full[[k]], fitted=pfit$wFitted.full[[k]],
                residuals = pfit$wResiduals.full[[k]], 
                weights = pfit$weights, Terms = pfit$Terms, term.labels = pfit$term.labels,
                data = pfit$data,
                SS = anova.parts.obs$SS, SS.type = SS.type, df = anova.parts.obs$df, 
                R2 = anova.parts.obs$R2[1:k], F = anova.parts.obs$Fs[1:k], permutations = iter+1,
                random.SS = SS, random.F = Fs, random.cohenf = cohenf, effect.type=effect.type,
                perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"))
  } else {
    SS <- do.call(SS.iter.null, SS.args)
    P.val <- pval(SS)
    Z <- effect.size(SS)
    n <- NROW(Y)
    df <- n - 1
    tab <- data.frame(Df = df,SS = SS[1],
                      MS = SS[1]/df, Rsq = 1,
                      F = NA, P = P.val)
    rownames(tab) <- "Residuals"
    colnames(tab)[NCOL(tab)] <- "Pr(>SS)"
    class(tab) = c("anova", class(tab))
    out <- list(aov.table = tab, call = match.call(),
                coefficients=pfit$wCoefficients.full[[1]], 
                Y=pfit$Y,  X=pfit$X, 
                QR = pfit$wQRs.full[[1]], fitted=pfit$wFitted.full[[1]],
                residuals = pfit$wResiduals.full[[1]], 
                weights = pfit$weights, Terms = pfit$Terms, term.labels = pfit$term.labels,
                data = pfit$data,
                random.SS = SS, 
                GLS = FALSE,
                OLS = TRUE)
  }
  if(!is.null(Pcov)) {
    PY <- as.matrix(crossprod(Pcov, pfit$Y))
    PX <- as.matrix(crossprod(Pcov, pfit$X))
    glsFit <- lm.wfit(PX, PY, pfit$weights)
    out$gls.coefficients <- glsFit$coefficients
    out$gls.fitted <- pfit$X %*% glsFit$coefficients
    out$gls.residuals <- pfit$Y - pfit$X %*% glsFit$coefficients
    GLS.mean = apply(PY, 2, mean)
    out$Cov = Cov
    out$Pcov = Pcov
    out$GLS = TRUE
    out$OLS = FALSE
  }
  class(out) <- "procD.lm"
  out
}