#' Procrustes ANOVA and pairwise tests for shape data, using complex linear models
#'
#' The function quantifies the relative amount of shape variation explained by a suite of factors
#' and covariates in a "full" model, after accounting for variation in a "reduced" model. Inputs are 
#' formulae for full and reduced models (order is not important, but it is better to list the model 
#' with the most terms first or use a geomorph data frame), plus indication if means or slopes 
#' are to be compared among groups, with appropriate formulae to define how they should be compared.
#' 
#'   The response matrix 'y' can be in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]) or a 3D array (p x k x n). It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The names specified for the 
#'   independent (x) variables in the formula represent one or more 
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the 
#'   shape matrix matches the order of values in the independent variables. Linear model fits (using the  \code{\link{lm}} function)
#'   can also be input in place of a formula.  Arguments for \code{\link{lm}} can also be passed on via this function.
#'
#'   The function performs statistical assessment of the terms in the model using Procrustes distances among 
#'   specimens, rather than explained covariance matrices among variables. With this approach, the sum-of-squared 
#'   Procrustes distances are used as a measure of SS (see Goodall 1991). The SS between models is evaluated through 
#'   permutation. In morphometrics this approach is known as a Procrustes ANOVA (Goodall 1991), which is equivalent
#'   to distance-based anova designs (Anderson 2001). Unlike \code{\link{procD.lm}}, this function is strictly for comparison
#'   of two nested models. (Use of \code{\link{procD.lm}} will be more suitable in most cases.)  
#'   A residual randomization permutation procedure (RRPP) is utilized 
#'   for reduced model residuals to evaluate the SS between models (Collyer et al. 2015).  Effect-sizes (Z-scores) are 
#'   computed as standard deviates of the SS or pairwise statistic sampling 
#'   distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2015).  If a phylogeny is 
#'   used, the ANOVA Z-score is calculated from the sampling distributions of the F value, as the total SS will vary among permutations.
#'   For ANOVA Z-scores, a log-transformation is performed first, to assure a norammly distributed sampling distribution.
#'   
#'   Pairwise tests are only performed if formulae are provided to compute such results.
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{advanced.procD.lm}}.
#'   The generic function, \code{\link{plot}}, produces diagnostic plots for Procrustes residuals of the linear fit.
#'   
#'  \subsection{Notes for geomorph 3.0.4 and subsequent versions}{ 
#'  Compared to previous versions of geomorph, users might notice differences in effect sizes.  Previous versions used z-scores calculated with 
#'  expected values of statistics from null hypotheses (sensu Collyer et al. 2015); however Adams and Collyer (2016) showed that expected values 
#'  for some statistics can vary with sample size and variable number, and recommended finding the expected value, empirically, as the mean from the set 
#'  of random outcomes.  Geomorph 3.0.4 and subsequent versions now center z-scores on their empirically estimated expected values and where appropriate, 
#'  log-transform values to assure statistics are normally distributed.  This can result in negative effect sizes, when statistics are smaller than 
#'  expected compared to the avergae random outcome.  For ANOVA-based functions, the option to choose among different statistics to measure effect size 
#'  is now a function argument.
#'  
#'  An optional argument for including a phylogenetic tree of {class phylo} is included in this function.  ANOVA performed on separate PGLS models is analagous
#'  to a likelihood ratio test between models (Adams and Collyer 2017).  Pairwise tests can also be performed after PGLS estimation of coefficients but users
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
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape (optional)
#' @param pc.shape An argument for whether analysis should be performed on the principal component scores of shape.  This is a useful
#' option if the data are high-dimensional (many more variables than observations) but will not affect results
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
#' @seealso \code{\link{procD.lm}}
#' @return Function returns an ANOVA table of statistical results for model comparison: error df (for each model), SS, MS,
#' F ratio, Z, and Prand.  A list of essentially the same components as \code{\link{procD.lm}} is also returned, and additionally
#' LS means or slopes, pairwise differences comparisons of these, effect sizes, and P-values may also be returned.  If a group formula
#' is provided but slope formula is null, pairwise differences are Procrustes distances between least squares (LS) means for the 
#' defined groups.  If a slope formula is provided, two sets of pairwise differences, plus effect sizes and P-values, are provided.  
#' The first is for differences in slope vector length (magnitude).  The length of the slope vector corresponds to the amount of shape
#' change per unit of covariate change.  Large differences correspond to differences in the amount of shape change between groups.  
#' The second is for slope vector orientation differences.  Differences in the direction of shape change (covariance of shape variables)
#' can be summarized as a vector correlation or angle between vectors.  See \code{\link{summary.advanced.procD.lm}} for summary options.

#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C. and M.L. Collyer. 2017. Multivariate comparative methods: evaluations, comparisons, and
#' recommendations. Systematic Biology. In press.
#' 
#' @examples
#'data(plethodon)
#'Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)
#'
#'# Example of a nested model comparison (as with ANOVA with RRPP)
#'advanced.procD.lm(coords ~ log(Csize) + species, 
#'~ log(Csize)*species*site, iter=499, data = gdf)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons 
#'advanced.procD.lm(coords ~ site*species, ~ site + species, groups = ~site*species, 
#'    iter=499, data = gdf)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons, 
#'# accounting for a common allometry  
#'advanced.procD.lm(coords ~ Csize + site*species, 
#'~ log(Csize) + site + species, 
#'groups = ~ site*species, slope = ~log(Csize), iter = 499, data = gdf)
#'
#'# Example of a test of homogeneity of slopes, plus pairwise slopes comparisons
#'gdf$group <- factor(paste(gdf$species, gdf$site, sep="."))
#'advanced.procD.lm(coords ~ log(Csize) + group, 
#'~ log(Csize) * group, 
#'groups = ~ group, 
#'slope = ~ log(Csize), angle.type = "deg", iter = 499, data = gdf)
#'
#'# Example of partial pairwise comparisons, given greater model complexity.
#'# Plus, working with class advanced.procD.lm objects.
#'aov.pleth <- advanced.procD.lm(coords ~ log(Csize)*site*species, 
#'~ log(Csize) + site*species, 
#'groups = ~ species, slope = ~ log(Csize), angle.type = "deg", iter = 499, data = gdf)
#'
#'summary(aov.pleth) # ANOVA plus pairwise tests
#'plot(aov.pleth) # diagnostic plots
#'aov.pleth$slopes # extract the slope vectors

advanced.procD.lm<-function(f1, f2, groups = NULL, slope = NULL, 
                            angle.type = c("r", "deg", "rad"), 
                            phy = NULL,
                            pc.shape = FALSE,
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
  if(!is.null(data)) dat.temp <- gdf.to.df(data) else 
    dat.temp <- pfit1$data
  dat.temp$Y <- Y
  if(class(f2) == "formula" && length(f2) == 3) f2 <- f2[-2]
  if(class(f2) == "formula") f2 <- update(f2, Y ~.) 
  if(pc.shape) pfit2 = procD.fit(f2, pca = TRUE, data = dat.temp, ...) else 
    pfit2 = procD.fit(f2, pca = FALSE, data = dat.temp, ...)
  if(!is.null(phy)){
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
  } else Pcor <- NULL
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
  Xf <- as.matrix(pfitf$Xfs[[kf]])
  if(is.null(pfitr$Xrs)) Xr <- matrix(1, n) else Xr <- as.matrix(pfitr$Xrs[[kr]])
  Er <- pfitr$residuals.full[[kr]]
  Ef <- pfitf$residuals.full[[kf]]
  wEr <- pfitr$wResiduals.full[[kr]]
  wEf <- pfitf$wResiduals.full[[kf]]
  Yhr <- pfitr$fitted.full[[kr]]
  Yhf <- pfitf$fitted.full[[kf]]
  wYhr <- pfitr$wFitted.full[[kr]]
  wYhf <- pfitf$wFitted.full[[kf]]
  w1 <- pfitr$weights; w2 <- pfitf$weights
  if(any(is.na(match(w1,w2)))) stop("Weights cannot be different between models.")
  w <- w1
  Qr <-qr.Q(qr(Xr*sqrt(w)))
  Qf <-qr.Q(qr(Xf*sqrt(w)))
  if(!is.null(phy)){
    PXf <- Pcor%*%Xf
    PXr <- Pcor%*%Xr
    PQr <-qr.Q(qr(PXr*sqrt(w)))
    PQf <-qr.Q(qr(PXf*sqrt(w)))
    PY <- Pcor%*%Y
    SSEr <- sum(fastLM(PQr, PY*sqrt(w))$residuals^2)
    SSEf <- sum(fastLM(PQf, PY*sqrt(w))$residuals^2)
    SSY <- sum(.lm.fit(Pcor%*%matrix(sqrt(w)), PY*sqrt(w))$residuals^2)
  } else {
    SSEr <- sum(wEr^2)
    SSEf <- sum(wEf^2)
    SSY <- sum(.lm.fit(matrix(sqrt(w)),  pfitf$wY)$residuals^2)
  }
  SSm <- SSEr - SSEf
  Fs <- (SSm/(dfr-dff))/(SSEf/dff)
  ind <- perm.index(n, iter, seed=seed)
  if(!is.null(groups) && !is.null(slope)) pairwise.cond <- "slopes"
  if(!is.null(groups) && is.null(slope)) pairwise.cond <-"means"
  if(is.null(groups) && is.null(slope)) pairwise.cond <-"none"
  if(is.null(groups) && !is.null(slope)) {
    cat("\nNo groups for which to compare means or slopes.  
        No pairwise tests will be performed\n")
    pairwise.cond <-"none" 
  }
  data.types <- lapply(data, class)
  keep = sapply(data.types, function(x) x != "array" & x != "phylo"  & x != "dist")
  if(!is.null(data)) dat2 <- as.data.frame(data[keep]) else dat2 <- pfitf$data
  
  if(!is.null(groups)){
    g.match <- match(names(dat2), attr(terms(groups), "term.labels"))
    if(!all(is.na(g.match))) gps <- dat2[,which(!is.na(g.match))] else
      gps <- model.frame(groups, data = dat2)
  } else gps <- NULL
  if(!is.null(slope)){
    s.match <- match(names(dat2), attr(terms(slope), "term.labels"))
    if(!all(is.na(s.match))) slp <- dat2[,which(!is.na(s.match))] else
      slp <- model.frame(slope, data = dat2)
  } else slp <- NULL
  
  if(pairwise.cond == "none"){
    if(print.progress)
      pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
    P <- NULL
    jj <- iter+1
    step <- 1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    while(jj > 0){
      ind.j <- ind[j]
      P <- c(P, lapply(1:length(j), function(i){
        if(!is.null(phy)) {
          y <- Pcor%*%(Er[ind[[i]],] + Yhr)*sqrt(w)
          ssr <- sum(fastLM(PQr,y)$residuals^2)
          ssf <- sum(fastLM(PQf,y)$residuals^2)
          ((ssr-ssf)/(dfr-dff))/(ssf/dff)
        } else {
          y <- (Er[ind[[i]],] + Yhr)*sqrt(w)
          sum((fastFit(Qf, y, n, p)- fastFit(Qr, y, n, p))^2)
        }
      }))
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
      if(print.progress) setTxtProgressBar(pb,step)
      step <- step+1
    }
    P <- simplify2array(P)
    P.val <- pval(P) 
    Z.score <- effect.size(log(P))
  }
  
  if(pairwise.cond == "means") {
    P <- lsms <- NULL
    if(print.progress)
      pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3) 
    jj <- iter+1
    step <- 1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    while(jj > 0){
      ind.j <- ind[j]
      Yr <- lapply(1:length(j), function(i) (Er[ind.j[[i]],] + Yhr)*sqrt(w))
      if(!is.null(phy)) {
        Yr <- lapply(1:length(j), function(i) Pcor%*%Yr[[i]])
        P <- c(P, lapply(1:length(j), function(i) {
          y <- Yr[[i]]
          ssr <- sum(fastLM(PQr,y)$residuals^2)
          ssf <- sum(fastLM(PQf,y)$residuals^2)
          ((ssr-ssf)/(dfr-dff))/(ssf/dff)
        }))
      } else 
        P <- c(P, lapply(1:length(j), function(i) 
          sum((fastFit(Qf, Yr[[i]], n, p)- fastFit(Qr, Yr[[i]], n, p))^2)))
      lsms <- c(lsms, apply.ls.means(pfitf, Yr, g = gps, data = dat2, Pcor = Pcor)) 
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
      if(print.progress) setTxtProgressBar(pb,step)
      step <- step+1
    }
    if(print.progress) close(pb)
    P <- simplify2array(P)
    P.dist <- lapply(1:length(lsms), function(j){as.matrix(dist(lsms[[j]]))})
    P.val <- pval(P) 
    Z.score <- effect.size(log(P))
    Means.dist <- P.dist[[1]]
    P.dist.s <- simplify2array(P.dist)
    P.Means.dist <- Pval.matrix(P.dist.s)
    Z.Means.dist <- Effect.size.matrix(P.dist.s)
  }
  
  if(pairwise.cond == "slopes") {
    P <- g.slopes <- NULL
    if(print.progress)
      pb <- txtProgressBar(min = 0, max = ceiling(iter/100), initial = 0, style=3)
    jj <- iter+1
    step <- 1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    while(jj > 0){
      ind.j <- ind[j]
      Yr <- lapply(1:length(j), function(i) (Er[ind.j[[i]],] + Yhr)*sqrt(w))
      if(!is.null(phy)) {
        Yr <- lapply(1:length(j), function(i) Pcor%*%Yr[[i]])
        P <- c(P, lapply(1:length(j), function(i) {
          y <- Yr[[i]]
          ssr <- sum(fastLM(PQr,y)$residuals^2)
          ssf <- sum(fastLM(PQf,y)$residuals^2)
          ((ssr-ssf)/(dfr-dff))/(ssf/dff)
        }))
      } else P <- c(P, lapply(1:length(j), function(i) 
        sum((fastFit(Qf, Yr[[i]],n,p)- fastFit(Qr, Yr[[i]],n,p))^2)))
      g.slopes <- c(g.slopes, apply.slopes(pfitf, g=gps, slope=slp, Yr, data=dat2, Pcor = if(is.null(Pcor)) NULL else Pcor)) 
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
      if(print.progress) setTxtProgressBar(pb,step)
      step <- step+1
    }
    if(print.progress) close(pb)
    P <- simplify2array(P)
    slope.lengths <- Map(function(y) sqrt(diag(tcrossprod(y))), g.slopes) 
    P.slopes.dist <- Map(function(y) as.matrix(dist(matrix(y))), slope.lengths) 
    P.cor <- Map(function(y) vec.cor.matrix(y), g.slopes) 
    P.val <- pval(P) 
    Z.score <- effect.size(log(P))
    P.slopes.dist.s <-simplify2array(P.slopes.dist)
    P.val.slopes.dist <- Pval.matrix(P.slopes.dist.s)
    Z.slopes.dist <- Effect.size.matrix(P.slopes.dist.s)
    P.cor.t <- 1 - simplify2array(P.cor)
    P.val.cor <- Pval.matrix(P.cor.t)
    Z.cor <- Effect.size.matrix(P.cor.t)
  }
  
  anova.table <- data.frame(df = c(dfr,dff), SSE = c(SSEr, SSEf), SS = c(NA, SSm),
                            R2 = c(NA, SSm/SSY), F = c(NA, Fs), Z = c(NA, Z.score), P = c(NA,P.val))
  rownames(anova.table) <- c(formula(pfitr$Terms), formula(pfitf$Terms))
  colnames(anova.table)[1] <- "Df"
  colnames(anova.table)[ncol(anova.table)] <- "Pr(>F)"
  class(anova.table) <- c("anova", class(anova.table))
  
  if(pairwise.cond == "slopes") {
    angle.type = match.arg(angle.type)
    random.angles <- acos(simplify2array(P.cor))
    angles.obs <-random.angles[,,1]
    diag(angles.obs) <- 0
    cor.obs <- P.cor[[1]]
    if(angle.type == "deg") {
      random.angles <- random.angles*180/pi
      angles.obs <-random.angles[,,1]
      diag(angles.obs) <- 0
    } 
    obs.slope.lengths <- slope.lengths[[1]]
    obs.slope.dist <- as.matrix(dist(obs.slope.lengths))
    dimnames(P.val.slopes.dist) <- dimnames(Z.slopes.dist) <- dimnames(obs.slope.dist)
  }
  
  if(pairwise.cond == "none"){
    if(is.null(phy)) {
      random.SS = P 
      random.F = NULL
    } else {
      random.SS = NULL
      random.F = P
    }
    out <- list(anova.table = anova.table, 
                coefficients=pfitf$wCoefficients.full,
                Y=pfitf$Y, X=pfitf$X,
                QR = pfitf$wQRs.full[[kf]], 
                fitted = wYhf,
                residuals = wEf,
                weights = w, data = dat2, random.SS = random.SS, random.F = random.F,
                Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
                call= match.call()
    )
  }
  if(pairwise.cond == "means"){
    if(is.null(phy)) {
      random.SS = P 
      random.F = NULL
    } else {
      random.SS = NULL
      random.F = P
    }
    out <- list(anova.table = anova.table, LS.means = lsms[[1]], 
                LS.means.dist = Means.dist, Z.means.dist = Z.Means.dist, P.means.dist = P.Means.dist, 
                coefficients=pfitf$wCoefficients.full, 
                Y=pfitf$Y, X=pfitf$X, 
                QR = pfitf$wQR.full[[kf]], fitted = wYhf, 
                residuals = wEf, 
                weights = w, data = dat2, random.SS = random.SS, random.F = random.F, random.means.dist = P.dist,
                Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
                call= match.call()
    )
  }
  if(pairwise.cond == "slopes"){
    if(is.null(phy)) {
      random.SS = P 
      random.F = NULL
    } else {
      random.SS = NULL
      random.F = P
    }
    if(angle.type == "r"){
      out <- list(anova.table = anova.table, slopes = g.slopes[[1]], slope.lengths = obs.slope.lengths,
                  slopes.dist = obs.slope.dist, P.slopes.dist = P.val.slopes.dist,
                  Z.slopes.dist = Z.slopes.dist,
                  slopes.cor = P.cor[[1]], P.slopes.cor = P.val.cor, Z.slopes.cor = Z.cor,
                  random.slopes = g.slopes, random.slopes.dist = P.slopes.dist, random.slopes.cor = P.cor,
                  coefficients=pfitf$wCoefficients.full, 
                  Y=pfitf$Y, X=pfitf$X, 
                  QR = pfitf$wQR.full[[kf]], fitted = wYhf, 
                  residuals = wEf, 
                  weights = w, data = dat2, random.SS = random.SS, random.F = random.F, 
                  Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
                  call= match.call()
      )
    } else {
      out <- list(anova.table = anova.table, slopes = g.slopes[[1]], slope.lengths = obs.slope.lengths,
                  slopes.dist = obs.slope.dist, P.slopes.dist = P.val.slopes.dist,
                  Z.slopes.dist = Z.slopes.dist,
                  slopes.angles = angles.obs, P.angles = P.val.cor, Z.angles = Z.cor,
                  random.slopes = g.slopes, random.slopes.dist = P.slopes.dist, random.angles = random.angles,
                  coefficients=pfitf$wCoefficients.full, 
                  Y=pfitf$Y, X=pfitf$X, 
                  QR = pfitf$wQR.full[[kf]], fitted = wYhf, 
                  residuals = wEf, 
                  weights = w, data = dat2, random.SS = random.SS, random.F = random.F, 
                  Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
                  call= match.call()
      )
    }
  }
  class(out) <- "advanced.procD.lm"
  out
}
