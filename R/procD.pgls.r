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
#' independent variables (discrete or continuous). The response matrix 'Y' can be either in the form of a two-dimensional data 
#' matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that the landmarks have previously 
#' been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].
#' Linear model fits (using the  \code{\link{lm}} function)
#' can also be input in place of a formula.  Arguments for \code{\link{lm}} can also be passed on via this function.
#' The user must also specify a phylogeny describing the evolutionary relationships among species (of class phylo).
#' Note that the specimen labels for both X and Y must match the labels on the tips of the phylogeny.
#'
#'   The function \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'   coordinates; however this step is no longer necessary, as procD.lm can receive 3D arrays as dependent variables.  It is also 
#'   recommended that \code{\link{geomorph.data.frame}} is used to create and input a data frame.  This will reduce problems caused
#'   by conflicts between the global and function environments.  In the absence of a specified data frame, procD.pgls will attempt to 
#'   coerce input data into a data frame, but success is not guaranteed.
#'   
#'   From the phylogeny, a phylogenetic transformation matrix is obtained under a Brownian motion model, and used to 
#'   transform the X and Y variables. Next, the Gower-centered distance matrix is obtained from predicted values from the
#'   model (Y ~ X), from which sums-of-squares, F-ratios, and R-squared are estimated for each factor in the model (see Adams, 2014). 
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
#'   they be present (see Anderson and terBraak 2003). 
#'   
#'   Effect-sizes (Z scores) are computed as standard deviates of either the 
#'   F or Cohen's f-squared sampling distributions generated, which might be more intuitive for P-values than F-values 
#'   (see Collyer et al. 2015).  Values from these distributions are log-transformed prior to effect size estimation,
#'   to assure normally distributed data.  The SS type will influence how Cohen's f-squared values are calculated.  
#'   Cohen's f-squared values are based on partial eta-squared values that can be calculated sequentially or marginally, as with SS.
#'   
#'   In the case  
#'   that multiple factor or factor-covariate interactions are used in the model formula, one can specify whether all main effects should be  
#'    added to the model first, or interactions should precede subsequent main effects 
#'   (i.e., Y ~ a + b + c + a:b + ..., or Y ~ a + b + a:b + c + ..., respectively.)
#'
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{procD.pgls}}.
#'   The generic function, \code{\link{plot}}, produces diagnostic plots for Procrustes residuals of the linear fit.
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
#' @param f1 A formula for the linear model (e.g., y ~ x1 + x2)
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param int.first A logical value to indicate if interactions of first main effects should precede subsequent main effects
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing.
#' @param mod.resid A logical value indicating whether exchangebale units in RRPP are modified residuals (if TRUE)
#' or unmodified residuals (if FALSE).  If RRPP is false, randomizing raw values is the same as randomizing
#' residuals of the overall mean.  Thus, modified residuals are the same as transformed values.
#' @param effect.type One of "F" or "cohen", to choose from which random distribution to estimate effect size.
#' (The default is "F".  The option, "cohen", refers to Cohen's f-squared values. 
#' Values are log-transformed before z-score calculation to assure normally distributed effect sizes.)
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
#' @return procD.lm.pgls returns an object of class "procD.lm".  
#' See \code{\link{procD.lm}} for a description of the list of results generated.  Additionally, procD.pgls provides
#' the phylogenetic correction matrix, Pcor, plus "pgls" adjusted coefficients, fitted values, residuals, and mean.
#' @references Adams, D.C. 2014. A method for assessing phylogenetic least squares models for shape and other high-dimensional 
#' multivariate data. Evolution. 68:2675-2688. 
#' @references Adams, D.C., and M.L. Collyer. 2015. Permutation tests for phylogenetic comparative analyses of high-dimensional 
#' shape data: what you shuffle matters. Evolution. 69:823-829.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
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
#' pleth.pgls$pgls.fitted #the PGLS fitted values 
procD.pgls<-function(f1, phy, iter=999, seed=NULL, int.first = FALSE, 
                     effect.type = c("F", "cohen"),
                     RRPP=TRUE, mod.resid = TRUE, data=NULL, print.progress = TRUE, ...){
  if(int.first==TRUE) ko = TRUE else ko = FALSE
  if(!is.null(data)) data <- droplevels(data)
  pfit <- procD.fit(f1, data=data, keep.order=ko, pca=FALSE, ... )
  n <- dim(pfit$Y)[[1]]
  p <- dim(pfit$Y)[[2]]
  k <- length(pfit$term.labels)
  Y <- as.matrix(pfit$wY)
  if(p > n) pfitr <- procD.fit(f1, data=data, keep.order=ko,  pca=TRUE, ...) else
    pfitr <- pfit
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
  Pcor <- Pcor[rownames(Y),rownames(Y)]
  if(RRPP && mod.resid) method = "rrpp.trans"
  if(RRPP && !mod.resid) method = "rrpp.std"
  if(!RRPP && mod.resid) method = "frpp.trans"
  if(!RRPP && !mod.resid) method = "frpp.std"
  
  if(k > 0) {
    if(print.progress) {
      if(method == "rrpp.std") P <- SS.pgls.iter(pfitr, Pcor=Pcor, iter=iter, seed=seed, Yalt="RRPP")
      if(method == "frpp.std") P <- SS.pgls.iter(pfitr, Pcor=Pcor, iter=iter, seed=seed, Yalt="resample")
      if(method == "rrpp.trans") P <- SS.pgls.iter.trans(pfitr, Pcor=Pcor, iter=iter, seed=seed, Yalt="RRPP")
      if(method == "frpp.trans") P <- SS.pgls.iter.trans(pfitr, Pcor=Pcor, iter=iter, seed=seed, Yalt="resample")
    } else {
      if(method == "rrpp.std") P <- .SS.pgls.iter(pfitr, Pcor, Yalt="RRPP",  iter=iter, seed=seed) 
      if(method == "frpp.std") P <- .SS.pgls.iter(pfitr, Pcor, Yalt="resample",  iter=iter, seed=seed)
      if(method == "rrpp.trans") P <- .SS.pgls.iter.trans(pfitr, Pcor, Yalt="RRPP", iter=iter, seed=seed)
      if(method == "frpp.trans") P <- .SS.pgls.iter.trans(pfitr, Pcor, Yalt="resample", iter=iter, seed=seed)
    }
    anova.parts.obs <- anova.parts(pfitr, P)
    anova.tab <-anova.parts.obs$anova.table 
    df <- anova.parts.obs$df
    SS <- P$SS
    SSE <- P$SSE
    MS <- SS/df[1:k]
    MSE <- SSE/df[k+1]
    SSE.mat <- matrix(SSE, k, length(SSE), byrow = TRUE)
    MSE.mat <- matrix(MSE, k, length(MSE), byrow = TRUE)
    SSY <- P$SSY
    effect.type <- match.arg(effect.type)
    SS.type <- pfit$SS.type
    if(is.matrix(SS)){
      Fs <- (SS/df[1:k])/MSE.mat
      if(SS.type == "III") {
        etas <- SS/(SS+SSE.mat)
        cohenf <- etas/(1-etas)
      } else {
        etas <- SS/SSY
        unexp <- 1 - apply(etas, 2, cumsum)
        cohenf <- etas/unexp
      }
      P.val <- apply(Fs, 1, pval)
      if(effect.type == "F") Z <- apply(log(Fs), 1, effect.size) else
        Z <- apply(log(cohenf), 1, effect.size) 
      rownames(SS) <- rownames(Fs) <- rownames(cohenf) <- pfit$term.labels
      colnames(SS) <- colnames(Fs) <- colnames(cohenf) <- c("obs", paste("iter", 1:iter, sep=":"))
    } else {
      MSE <- SSE/df[2]
      Fs <- (SS/df[1])/MSE
      etas <- SS/SSY
      cohenf <- etas/(1-etas)
      P.val <- pval(Fs)
      if(effect.type != "F" && effect.type != "cohen") {
        effect.type <- "F"
        cat("\nWarning: only F or Cohen's f-squared can be used for effect sizes
            with PGLS.  Effect type has been changed to F.\n")
      }
      if(effect.type == "cohen") Z <- effect.size(log(cohenf)) else
        Z <- effect.size(log(Fs)) 
      names(SS) <- names(Fs) <- names(cohenf) <- c("obs", paste("iter", 1:iter, sep=":"))
      }
    if(effect.type == "SS") effect.type <- "F"
    tab <- data.frame(anova.tab, Z = c(Z, NA, NA), Pr = c(P.val, NA, NA))
    colnames(tab)[1] <- "Df"
    colnames(tab)[ncol(tab)] <- "Pr(>F)"
    if(effect.type == "cohen") colnames(tab)[ncol(tab)] <- "Pr(>Cohen f-sq)"
    class(tab) <- c("anova", class(tab))
    PY <- Pcor%*%pfit$Y; PX <- Pcor%*%pfit$X
    Pfit <- lm.wfit(PX, PY, pfit$weights)
    SS <- rbind(SS, SSE, SSY)
    rownames(SS) <- c(pfit$term.labels, "Residuals", "Total")
    out = list(aov.table = tab, call = match.call(),
               coefficients=pfit$wCoefficients.full[[k]],
               Y=pfit$Y,  X=pfit$X, 
               Pcor=Pcor, 
               QR = pfit$QRs[[k]], fitted=pfit$wFitted.full[[k]], 
               residuals = pfit$wResiduals.full[[k]], 
               weights = pfit$weights, Terms = pfit$Terms, term.labels = pfit$term.labels,
               SS = anova.parts.obs$SS, SS.type = SS.type,
               df = anova.parts.obs$df, R2 = anova.parts.obs$R2[1:k], 
               pgls.coefficients = Pfit$coefficients, 
               pgls.fitted = pfit$X%*%Pfit$coefficients, 
               pgls.residuals = Y - pfit$X%*%Pfit$coefficients,
               phylo.mean = apply(PY, 2, mean),
               F = anova.parts.obs$Fs[1:k], permutations = iter+1, random.SS = SS,
               random.SSE <- SSE,
               random.F = Fs, random.cohenf = cohenf, effect.type=effect.type,
               perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"), PGLS = TRUE)
  } else {
    Y <- pfit$Y
    PY <- crossprod(Pcor, Y) 
    X <- pfit$X
    PX <- crossprod(Pcor, X)
    Pfit <- lm.wfit(PX, PY, w = pfit$weights)
    Pcoef <- Pfit$coefficients
    Pfitted <- X%*%Pfit$coefficients
    Pres <- Y - Pfitted
    SSY <- sum(Pfit$residuals^2)
    n <- NROW(Y)
    df <- n - 1
    tab <- data.frame(Df = df,SS = SSY,
                      MS = SSY/df, Rsq = NA,
                      F = NA, P = NA)
    rownames(tab) <- "Residuals"
    colnames(tab)[NCOL(tab)] <- "Pr(>F)"
    class(tab) = c("anova", class(tab))
    out <- list(aov.table = tab, call = match.call(),
                coefficients=pfit$wCoefficients.full[[1]],
                Y=pfit$Y,  X=pfit$X, 
                Pcor=Pcor, 
                QR = pfit$QRs[[1]], fitted=pfit$wFitted.full[[1]], 
                residuals = pfit$wResiduals.full[[1]], 
                weights = pfit$weights, Terms = pfit$Terms, term.labels = pfit$term.labels,
                pgls.coefficients = Pcoef, 
                pgls.fitted = Pfitted,
                pgls.residuals = Pres,
                phylo.mean = apply(PY, 2, mean)
    )
  }
  class(out) <- "procD.lm"
  out
}
