#' Procrustes ANOVA and pairwise tests for shape data, using complex linear models
#' 
#' The function quantifies the relative amount of shape variation explained by a suite of factors
#' and covariates in a "full" model, after accounting for variation in a "reduced" model. Inputs are
#' formulae for full and reduced models (order is not important, but it is better to list the model
#' with the most terms first or use a geomorph data frame), plus indication if means or slopes
#' are to be compared among groups, with appropriate formulae to define how they should be compared.
#'
#'   This function calculates residual sum of squares either via ordinary least squares (OLS) estimation or
#'   phylogenetic least squares (PGLS) estimation for both full and reduced models.  Residuals from the reduced model are used
#'   in a randomized residual permutation procedure (RRPP) to find the difference in residual sum of squares (trace of the residual
#'   sums of squares and cross-products matrix, SSCP) over many permutations, thus creating a distribution of sum of squares (SS)
#'   for the parameters that differ between models (Collyer et al. 2015).  The SS can be converted to F-values to generate an empirical F-distribution.  
#'   A P-value is estimated as the percentile of the observed value in this distribution.
#'   
#'   The response matrix 'Y' can be in the form of a two-dimensional data
#'   matrix of dimension (n x [p x k]) or a 3D array (p x k x n). It is assumed that the landmarks have previously
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The names specified for the
#'   independent (x) variables in the formula represent one or more
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the
#'   shape matrix matches the order of values in the independent variables. Linear model fits (using the  \code{\link{lm}} function)
#'   can also be input in place of a formula.  Arguments for \code{\link{lm}} can also be passed on via this function.
#'
#'   The SS calculated is the same as the sum of squared Procrustes distances among specimens, as used as a measure of SS in Procrustes ANOVA (see Goodall 1991). 
#'   Procrustes ANOVA, often used in morphometrics applications is equivalent
#'   to distance-based anova designs (Anderson 2001). Unlike \code{\link{procD.lm}}, this function is strictly for comparison
#'   of two nested models. (Use of \code{\link{procD.lm}} will be more suitable in most cases.)
#'   Effect-sizes (Z-scores) are computed as standard deviates of the statistic chosen for ANOVA (see arguments) or for
#'   pairwise statistic sampling distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2015).
#'   For ANOVA Z-scores, a log-transformation is performed first, to assure a normally distributed sampling distribution.
#'   
#'   Pairwise tests have two flavors: 1) tests for differences in group means (based on vector length between
#'   means for pairwise comparisons) and 2) tests for angular differences in slopes between groups.  These tests are 
#'   similar in concept to trajectory analysis (Adams and Collyer 2007; Collyer and Adams 2007; Adams and Collyer 2009;
#'   Collyer and Adams 2013), in that pairwise statistics are either vector lengths or angular differences between vectors.  
#'   These tests are different than trajectory analysis (see\code{\link{trajectory.analysis}}), however, because a factorial model
#'   is not explicitly needed to contrast vectors between point factor levels nested within group factor levels.  For angular differences 
#'   between factor-covariate slopes, either the angle or the vector correlation can be tested.  It should be understood
#'   that a vector correlation of 1 (parallel vectors), not 0, is the null hypothesis, meaning slopes are the same.
#'
#'   Pairwise tests are only performed if formulae are provided to compute such results.
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{advanced.procD.lm}}.
#'   The generic function, \code{\link{plot}}, produces diagnostic plots for residuals of the linear fit.  Note that there is an
#'   argument in print/summary generic functions to print formulas as row names of the ANOVA table.  If
#'   formulas are long, it is recommended to make this argument, \code{formula = FALSE}, in which case 
#'   "reduced" and "full" models will be acknowledged.
#'   
#'  \subsection{Notes for geomorph 3.0.7 and subsequent versions}{
#'  The \code{advanced.procD.lm} function now defers to the R package, \code{RRPP}, specifically the \code{\link[RRPP]{anova.lm.rrpp}} and
#'  \code{\link[RRPP]{pairwise}} functions.  These functions perform all necessary computations needed for \code{advanced.procD.lm}, as well as other
#'  analyses.  Therefore, \code{advanced.procD.lm} is now a wrapper for these other functions. 
#'  The \code{\link[RRPP]{lm.rrpp}} function can be used for multiple models, if one wishes to work directly in \code{RRPP}, prior to using
#'   \code{\link[RRPP]{anova.lm.rrpp}} and \code{\link[RRPP]{pairwise}} functions.  The only difference in results (compared to version 3.0.6 and before) 
#'   should occur when comparing univariate slopes.  Version 3.0.6 and earlier versions appended a vector of 1s to slopes as an ad-hoc strategy to make 
#'   computations work.  This is no longer needed, as the \code{RRPP} functions can better handle univariate data.
#' }
#' 
#'  \subsection{Notes for geomorph 3.0.6 and subsequent versions}{
#'  For pairwise tests, previous versions assumed that pairwise comparisons of least-squares means used models with parallel slopes.
#'  Under most circumstances, this assumption is safe (and preferred), as the estimation of mean differences otherwise would have to 
#'  assume something about the mean values of covariates as appropriate locations for estimating means.  Version 3.0.6 and subseqent versions
#'  find least-squares means that are truer to the model defined.  For example, if a user defines a full model with parallel slopes, e.g.,
#'  shape ~ x + A + B + A:B, where x is a covariate and A and B are factors, results should be no different than before.  However, if a user 
#'  defines a full model which allows unique slopes, e.g., shape ~ x + A + B + x:A + x:B + A:B + x:A:B, least squares means will now be estimated
#'  for mean values of x using the coefficients for x:A, x:B, and x:A:B (previous versions did not).  This change is to made to
#'  be consistent with other least-squares means estimation functions in other packages.
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
#'
#'  An optional argument for including a phylogenetic tree of {class phylo} is included in this function.  ANOVA performed on separate PGLS models is analogous
#'  to a likelihood ratio test between models (Adams and Collyer 2018).  Pairwise tests can also be performed after PGLS estimation of coefficients but users
#'  should be aware that no formal research on the statistical properties (type I error rates and statistical power) of pairwise statistics with PGLS has yet
#'  been performed.  Using PGLS and analysis of pairwise statistics, therefore, assumes some risk.
#' }
#'
#' @param f1 A formula for a linear model, containing the response matrix (e.g., y ~ x1 + x2)
#' @param f2 A formula for another linear model (e.g., ~ x1 + x2 + x3 + a*b). f1 and f2 should be nested.
#' @param groups A formula for grouping factors (e.g., ~ a, or ~ a*b).  This argument should be left NULL unless one wishes to perform pairwise
#' comparisons of different group levels.  Note that this argument is used in conjunction with the argument, slope.  If slope is NULL, a pairwise
#' comparison test is performed on group least squares (LS) means.  If slope is not NULL, this argument will designate the group levels to compare
#' in terms of their slopes.
#' @param slope A formula with one - and only one - covariate (e.g., ~ x3).  This argument must be used in conjunction with the groups argument.  It
#' will not make sense if the groups argument is left NULL.  The groups argument defines the groups; the slope argument defines for which covariate group
#' slopes are compared.  Group slopes can differ in their magnitude and direction of shape change.
#' @param angle.type A value specifying whether directional differences between slopes should be represented by vector
#' correlations (r), radians (rad) or degrees (deg).
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape (optional).
#' @param Cov A covariance matrix to guide GLS computations.  If both a phylogenetic tree and covariance matrix
#' are provided, a Brownian motion covariance matrix will not be calculated, based on the phylogeny (it will be ignored).  
#' Currently this function cannot handle multiple covariance matrices.
#' @param pc.shape An argument for whether analysis should be performed on the principal component scores of shape.  This is a useful
#' option if the data are high-dimensional (many more variables than observations) but will not affect results.
#' @param effect.type An optional argument for which distribution of statistics should be used for calculating effect sizes (
#' and P-values).  The default is "F" for the distribution of random F-statistics, but "SS" and "Rsq" are also
#' possible, for the distributions of random SS between models or R-squared values, respectively.  One should not choose "SS" if a PGLS model is considered.
#' P-values should be similar in most cases, regardless of statistic chosen, as the rank correlations between statistics are either perfect (SS and Rsq for OLS) 
#' or generally large.
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param data A data frame for the function environment; see \code{\link{geomorph.data.frame}}.  If variables
#' are transformed in formulae, they should also be transformed in the geomorph data frame.  (See examples.)
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.
#' This is helpful for long-running analyses.
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function,
#' such as weights or offset).
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{procD.lm}}, \code{\link{procD.pgls}}, \code{\link{trajectory.analysis}},
#' @seealso \code{\link[RRPP]{lm.rrpp}}, \code{\link[RRPP]{anova.lm.rrpp}}, \code{\link[RRPP]{pairwise}} 
#' (used in some internal computations)
#' @return Function returns an ANOVA table of statistical results for model comparison: error df (for each model), SS, MS,
#' F ratio, Z, and Prand.  A list of essentially the same components as \code{\link{procD.lm}} is also returned, and additionally
#' LS means or slopes, pairwise differences comparisons of these, effect sizes, and P-values may also be returned.  If a group formula
#' is provided but slope formula is null, pairwise differences are Procrustes distances between least squares (LS) means for the
#' defined groups.  If a slope formula is provided, two sets of pairwise differences, plus effect sizes and P-values, are provided.
#' The first is for differences in slope vector length (magnitude).  The length of the slope vector corresponds to the amount of shape
#' change per unit of covariate change.  Large differences correspond to differences in the amount of shape change between groups.
#' The second is for slope vector orientation differences.  Differences in the direction of shape change (covariance of shape variables)
#' can be summarized as a vector correlation or angle between vectors.  See \code{\link{summary.advanced.procD.lm}} for summary options.

#' @references Adams, D.C., and M.L. Collyer. 2007. The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D.C., and M.L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C. and M.L. Collyer. 2018. Multivariate phylogenetic comparative methods: evaluations, comparisons, and
#' recommendations. Systematic Biology. 67:14-31.
#' @references Collyer, M.L., and D.C. Adams. 2007. Analysis of two-state multivariate phenotypic change 
#' in ecological studies. Ecology 88:683-692.
#' @references Collyer, M.L., and D.C. Adams. 2013. Phenotypic trajectory analysis: comparison of shape change patterns 
#' in evolution and ecology. Hystrix 24: 75-83.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#'
#' @examples
#'data(plethodon)
#'Y.gpa<-gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment
#'gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, 
#'site = plethodon$site)
#'
#'# Example of a nested model comparison (as with ANOVA with RRPP)
#'ANOVA <-  advanced.procD.lm(f1= coords ~ log(Csize) + species,
#'f2= ~ log(Csize)*species*site, iter=99, data = gdf)
#'summary(ANOVA, formula = FALSE) # formulas too long to print
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons
#'PW.means.test <- advanced.procD.lm(f1= coords ~ site*species, f2= ~ site + species, 
#'groups = ~site*species, iter=99, data = gdf)
#'summary(PW.means.test, formula = TRUE)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons,
#'# accounting for a common allometry
#'PW.ls.means.test <- advanced.procD.lm(f1= coords ~ Csize + site*species,
#'f2= ~ log(Csize) + site + species,
#'groups = ~ site*species, iter = 99, data = gdf)
#'summary(PW.ls.means.test, formula = TRUE)
#'
#'# Example of a test of homogeneity of slopes, plus pairwise slopes comparisons
#'gdf$group <- factor(paste(gdf$species, gdf$site, sep="."))
#'HOS <- advanced.procD.lm(f1= coords ~ log(Csize) + group,
#'f2= ~ log(Csize) * group, groups = ~ group,
#'slope = ~ log(Csize), angle.type = "deg", iter = 99, data = gdf)
#'summary(HOS, formula = FALSE) # formulas too long to print
#'
#'# Example of partial pairwise comparisons, given greater model complexity.
#'# Plus, working with class advanced.procD.lm objects.
#'aov.pleth <- advanced.procD.lm(f1= coords ~ log(Csize)*site*species,
#'f2= ~ log(Csize) + site*species, groups = ~ species, 
#'slope = ~ log(Csize), angle.type = "deg", iter = 99, data = gdf)
#'
#'summary(aov.pleth, formula = FALSE)  # formulas too long to print
#'
#'# Diagnostic plots
#'plot(aov.pleth) 
#'
#'# Extracting objects from results
#'aov.pleth$slopes # extract the slope vectors
#'
#'# GLS Examples (same as procD.gpls example)
#'data(plethspecies)
#'Y.gpa<-gpagen(plethspecies$land)    
#'gdf <- geomorph.data.frame(Y.gpa, tree = plethspecies$phy)
#'procD.pgls(coords ~ Csize, phy = tree, data = gdf, iter = 999)
#'
#'advanced.procD.lm(coords ~ Csize, ~1, phy = gdf$tree, data = gdf, iter = 999)
#'
#' # Could also do this with ape function
#' # phyCov <- vcv.phylo(plethspecies$phy)
#' # advanced.procD.lm(coords ~ Csize, ~1, Cov = phyCov, data = gdf, iter = 999)
#' 
advanced.procD.lm<-function(f1, f2, groups = NULL, slope = NULL,
                            angle.type = c("r", "deg", "rad"),
                            phy = NULL, Cov = NULL,
                            pc.shape = FALSE,
                            effect.type = c("F", "SS", "Rsq"),
                            iter=999,
                            seed = NULL,
                            print.progress = TRUE,
                            data=NULL, ...){
  
  if(!is.null(data)) data <- droplevels(data)
  if(pc.shape) pfit1 <- procD.fit(f1, data=data, pca = TRUE, ...) else
    pfit1 <- procD.fit(f1, data=data, pca=FALSE, ...)
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  Y <- as.matrix(pfit1$Y)
  n <- dim(Y)[1]; p <- dim(Y)[2]
  
  # covariance matrices
  if(is.null(Cov) && is.null(phy)) Cov <- NULL
  if(!is.null(phy)) {
    phy.name <- deparse(substitute(phy))
    if(!inherits(phy, "phylo")) 
      stop(paste("No phylo object called,", phy.name,"found in data or gobal environment."))
    N <- length(phy$tip.label)
    if(is.null(Cov)) {
      if(length(match(rownames(Y), phy$tip.label))!=N)
        stop("Data matrix missing some taxa present on the tree.")
      if(length(match(phy$tip.label,rownames(Y)))!=N)
        stop("Tree missing some taxa in the data matrix.")
      C <- vcv.phylo(phy)
      id <- rownames(Y)
      Cov <- C[id, id]
    }
  }
  if(!is.null(Cov)) {
    id <- rownames(Y)
    Cov <- C[id, id]
 } 

  # initial model evaluations
  if(!is.null(data)) dat.temp <- gdf.to.df(data) else
    dat.temp <- pfit1$data
  dat.temp$Y <- Y
  if(class(f2) == "formula" && length(f2) == 3) f2 <- f2[-2]
  if(class(f2) == "formula") f2 <- update(f2, Y ~.)
  if(pc.shape) pfit2 = procD.fit(f2, pca = TRUE, data = dat.temp, ...) else
    pfit2 = procD.fit(f2, pca = FALSE, data = dat.temp, ...)
  k1 <- pfit1$QRs.full[[length(pfit1$QRs.full)]]$rank
  k2 <- pfit2$QRs.full[[length(pfit2$QRs.full)]]$rank
  if(k1 > k2) pfitf <- pfit1 else pfitf <- pfit2
  if(k1 > k2) pfitr <- pfit2 else pfitr <- pfit1
  if(k1 == k2) stop("Models have same df")
  dat <- pfitf$data
  kr <- length(pfitr$residuals.full)
  kf <- length(pfitf$residuals.full)
  dfr <- pfitr$QRs.full[[length(pfitr$QRs.full)]]$rank
  dff <- pfitf$QRs.full[[length(pfitf$QRs.full)]]$rank
  k.total <- kr+kf
  k.unique <- length(unique(c(pfitf$term.labels, pfitr$term.labels)))
  if(kr >1 & kf > 1 & k.unique == k.total) stop("Models are not nested")
  dfr <- NROW(pfitr$wResiduals.full[[kr]]) - dfr
  dff <- NROW(pfitf$wResiduals.full[[kf]]) - dff
  w1 <- pfitr$weights; w2 <- pfitf$weights
  if(any(is.na(match(w1,w2)))) stop("Weights cannot be different between models.")
  w <- sqrt(w1)
  ind <- perm.index(n, iter, seed=seed)
  if(!is.null(groups) && !is.null(slope)) pairwise.cond <- "slopes"
  if(!is.null(groups) && is.null(slope)) pairwise.cond <-"means"
  if(is.null(groups) && is.null(slope)) pairwise.cond <-"none"
  if(is.null(groups) && !is.null(slope)) {
    cat("\nNo groups for which to compare means or slopes.
        No pairwise tests will be performed\n")
    pairwise.cond <-"none"
  }
  
  # Start RRPP here
  
  if(print.progress) cat("\nObtaining initial full model fit:")
  pfitF <- lm.rrpp(formula(pfitf$Terms), data = pfitf$data, iter=iter, seed = seed, 
                   Cov = Cov, print.progress = print.progress)
  
  if(print.progress) cat("\nObtaining initial reduced model fit:")
  pfitR <- lm.rrpp(formula(pfitr$Terms), data = pfitr$data, iter=iter, seed = seed, 
                   Cov = Cov, print.progress = print.progress)
  

  # ANOVA Table
  effect.type <- match.arg(effect.type)
  
  if(print.progress) cat("\nPerforming ANOVA:")
  anova.res <- anova.lm.rrpp(pfitR, pfitF, 
                    effect.type = effect.type, print.progress = print.progress)
  
  # pairwise means
  if(pairwise.cond == "means") {
    GT <- attr(terms(groups, data = dat.temp), "term.labels")
    nms <- strsplit(GT, ":")
    nml <- sapply(nms, length)
    gp <- GT[which(nml == 1)]
    dat.gp <- dat.temp[names(dat.temp) %in% gp]
    gp <- factor(apply(dat.gp, 1, paste, collapse = "."))
    
    if(print.progress) cat("\nCalculating LS means:")
    PW <-  pairwise(pfitF, fit.null = pfitR, groups = gp, 
                    covariate = NULL, print.progress = print.progress)
  }
  

  # pairwise slopes
  if(pairwise.cond == "slopes") {
    GT <- attr(terms(groups, data = dat.temp), "term.labels")
    nms <- strsplit(GT, ":")
    nml <- sapply(nms, length)
    gp <- GT[which(nml == 1)]
    dat.gp <- dat.temp[names(dat.temp) %in% gp]
    gp <- factor(apply(dat.gp, 1, paste, collapse = "."))
    covariate <- model.frame(slope, data = dat.temp)
    if(NCOL(covariate) > 1) stop("This analysis is limited to one covariate (slope).") else
      covariate <- as.vector(covariate[[1]])
    PW <-  pairwise(pfitF, fit.null = pfitR, groups = gp, 
                    covariate = covariate, print.progress = print.progress)
  }

  # pairwise output
  
  if(pairwise.cond == "means") {
    lsms <- PW$LS.means
    P.dist <- PW$means.dist
    Means.dist <- P.dist[[1]]
    P.dist.s <- simplify2array(P.dist)
    P.Means.dist <- Pval.matrix(P.dist.s)
    Z.Means.dist <- Effect.size.matrix(P.dist.s)
  }

  if(pairwise.cond == "slopes") {
    angle.type = match.arg(angle.type)
    P.cor <- PW$slopes.vec.cor
    slope.lengths <- PW$slopes.length
    P.slopes.dist <- PW$slopes.dist
    P.slopes.dist.s <-simplify2array(P.slopes.dist)
    P.val.slopes.dist <- Pval.matrix(P.slopes.dist.s)
    Z.slopes.dist <- Effect.size.matrix(P.slopes.dist.s)
    P.cor.t <- 1 - simplify2array(P.cor)
    P.val.cor <- Pval.matrix(P.cor.t)
    Z.cor <- Effect.size.matrix(P.cor.t)
    options(warn = -1)
    random.angles <- acos(simplify2array(P.cor))
    angles.obs <-random.angles[,,1]
    diag(angles.obs) <- 0
    cor.obs <- P.cor[[1]]
    if(angle.type == "deg") {
      random.angles <- random.angles*180/pi
      angles.obs <-random.angles[,,1]
      diag(angles.obs) <- 0
    }
    options(warn = 0)
    obs.slope.lengths <- slope.lengths[[1]]
    obs.slope.dist <- as.matrix(dist(obs.slope.lengths))
    dimnames(P.val.slopes.dist) <- dimnames(Z.slopes.dist) <- dimnames(obs.slope.dist)
  }

# output
  out <- list(anova.table = anova.res$table,
              Y=pfitF$LM$Y, X=pfitF$LM$X,
              QR = pfitF$LM$QR,
              coefficients = pfitF$LM$coefficients,
              fitted = pfitF$LM$fitted,
              residuals = pfitF$LM$residuals,
              weights = pfitF$LM$weights, data = pfitF$LM$data, 
              random.SS = anova.res$SS, random.F = anova.res$F,
              effect.type = effect.type,
              OLS = TRUE, GLS = FALSE,
              Terms = pfitF$LM$Terms, term.labels = pfitF$LM$term.labels, 
              permutations = iter+1,
              call= match.call()
  )
  if(!is.null(Cov)) {
    out$PCov = pfitF$LM$Pcov
    out$GLS = TRUE
    out$OLS = FALSE
    gls.fitted = pfitF$LM$gls.fitted
    gls.residuals = pfitF$LM$gls.residuals
    gls.coefficients = pfitF$LM$gls.coefficients
  }
  
  if(pairwise.cond == "means"){
    out$LS.means <- lsms[[1]]
    out$random.LS.means <- lsms
    out$random.means.dist <- P.dist.s
    out$LS.obs.means.dist <- P.dist[[1]]
    out$Z.means.dist <- Z.Means.dist
    out$P.means.dist <- P.Means.dist
  }
  
  if(pairwise.cond == "slopes"){
    if(angle.type == "r"){
      out$slopes <- PW$slopes[[1]]
      out$obs.slope.lengths <- obs.slope.lengths
      out$obs.slopes.dist <- obs.slope.dist
      out$random.slopes <- PW$slopes
      out$random.slope.lengths <- slope.lengths
      out$random.slopes.dist <- P.slopes.dist
      out$P.slopes.dist <- P.val.slopes.dist
      out$Z.slopes.dist <- Z.slopes.dist
      out$obs.slopes.cor <- P.cor[[1]]
      out$random.cor <- P.cor
      out$P.slopes.cor <- P.val.cor
      out$Z.slopes.cor <- Z.cor
    
    } else {
      
      out$slopes <- PW$slopes[[1]]
      out$obs.slope.lengths <- obs.slope.lengths
      out$obs.slopes.dist <- obs.slope.dist
      out$random.slopes <- PW$slopes
      out$random.slope.lengths <- slope.lengths
      out$random.slopes.dist <- P.slopes.dist
      out$P.slopes.dist <- P.val.slopes.dist
      out$Z.slopes.dist <- Z.slopes.dist
      out$obs.slopes.angles <- angles.obs 
      out$random.angles <- random.angles
      out$P.angles <- P.val.cor
      out$Z.angles <- Z.cor
    }
  }
  rownames(out$anova.table) <- c(formula(pfitr$Terms)[-2],
                                 formula(pfitf$Terms)[-2], "Total")
  rownames(out$anova.table)[1] <- paste(rownames(out$anova.table)[1], "(Null)")
  Rsq.fix <- out$anova.table$Rsq 
  out$anova.table$Rsq <- c(NA, Rsq.fix[2] - Rsq.fix[1], NA)
  out$anova.table <- out$anova.table[,-2]
  
  class(out) <- "advanced.procD.lm"
  out
}
