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
#'   computed as standard deviates of the SS sampling 
#'   distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2015).  
#'   
#'   Pairwise tests are only performed if formulae are provided to compute such results.
#'   The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{advanced.procD.lm}}.
#'   The generic function, \code{\link{plot}}, produces diagnostic plots for Procrustes residuals of the linear fit.
#'
#' @param f1 A formula for a linear model, containing the response matrix (e.g., y ~ x1 + x2)
#' @param f2 A formula for another linear model (e.g., ~ x1 + x2 + x3 + a*b). f1 and f2 should be nested.
#' @param groups A formula for grouping factors (e.g., ~ a, or ~ a*b)
#' @param slope A formula with one covariate (e.g., ~ x3)
#' @param angle.type A value specifying whether differences between slopes should be represented by vector
#' correlations (r), radians (rad) or degrees (deg)
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param data A data frame for the function environment; see \code{\link{geomorph.data.frame}}.  If variables
#' are transformed in formulae, they should also be transformed in the geomorph data frame.  (See examples.)
#' @param ... Arguments passed on to procD.fit (typically associated with the lm function)
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{procD.lm}}
#' @return Function returns an ANOVA table of statistical results for model comparison: error df (for each model), SS, MS,
#' F ratio, Z, and Prand.  A list of essentially the same components as \code{\link{procD.lm}} is also returned, and additionally
#' LS means or slopes, pairwise differences comparisons of these, effect sizes, and P-values may also be returned.

#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @examples
#'data(plethodon)
#'Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'gdf <- geomorph.data.frame(Y.gpa, logcs = log(Y.gpa$Csize), 
#'species = plethodon$species, site = plethodon$site)
#'
#'# Example of a nested model comparison (as with ANOVA with RRPP)
#'advanced.procD.lm(coords ~ logcs + species, 
#'~ logcs*species*site, iter=999, data = gdf)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons 
#'advanced.procD.lm(coords ~ site*species, ~site + species, groups = ~site*species, 
#'    iter=199, data = gdf)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons, 
#'# accounting for a common allometry  
#'advanced.procD.lm(coords ~ Csize + site*species, 
#'~Csize + site + species, 
#'groups = ~site*species, slope = ~Csize, iter = 199, data = gdf)
#'
#'# Example of a test of homogeneity of slopes, plus pairwise slopes comparisons
#'advanced.procD.lm(coords ~ logcs, 
#'~logcs + site*species, 
#'groups = ~site*species, 
#'slope = ~logcs, angle.type = "deg", iter = 199, data = gdf)
#'
#'# Example of partial pairwise comparisons, given greater model complexity.
#'# Plus, working with class advanced.procD.lm objects.
#'aov.pleth <- advanced.procD.lm(coords ~ logcs*site*species, 
#'~logcs + site*species, 
#'groups = ~species, slope = ~logcs, angle.type = "deg", iter = 199, data = gdf)
#'
#'summary(aov.pleth) # ANOVA plus pairwise tests
#'plot(aov.pleth) # diagnostic plots
#'aov.pleth$slopes # extract the slope vectors

advanced.procD.lm<-function(f1, f2, groups = NULL, slope = NULL, 
                            angle.type = c("r", "deg", "rad"), iter=999, 
                            seed = NULL, data=NULL, ...){
  pfit1 <- procD.fit(f1, data=data)
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  Y <- as.matrix(pfit1$Y)
  if(!is.null(pfit1$weights)) w <- pfit1$weights else w <- rep(1,n)
  if(any(w < 0)) stop("Weights cannot be negative")
  n <- nrow(Y)
  if(!is.null(data)) data2 <- geomorph.data.frame(Y=Y, data) else
    data2 <- geomorph.data.frame(Y=Y, pfit1$data[,-(1:ncol(Y))])
  if(any(class(f2)=="lm")) pfit2 = procD.fit(f2) else {
      f2 <- update(f2, Y ~.) 
      if(!is.null(data)) pfit2= procD.fit(f2, data=data2) else pfit2= procD.fit(f2)
  }
  k1 <- pfit1$QRs[[length(pfit1$QRs)]]$rank
  k2 <- pfit2$QRs[[length(pfit2$QRs)]]$rank
  if(k1 > k2) pfitf <- pfit1 else pfitf <- pfit2
  if(k1 > k2) pfitr <- pfit2 else pfitr <- pfit1
  if(k1 == k2) stop("Models have same df")
  dat <- pfitf$data
  kr <- length(pfitr$residuals)
  kf <- length(pfitf$residuals)
  dfr <- pfitr$QRs[[length(pfitr$QRs)]]$rank
  dff <- pfitf$QRs[[length(pfitf$QRs)]]$rank
  k.total <- kr+kf-2
  k.unique <- length(unique(c(pfitf$term.labels, pfitr$term.labels)))
  if(kr >1 & kf > 1 & k.unique == k.total) stop("Models are not nested")
  dfr <- nrow(pfitr$wResiduals[[kr]]) - dfr
  dff <- nrow(pfitf$wResiduals[[kf]]) - dff
  SSEr <- sum(pfitr$wResiduals[[kr]]^2)
  SSEf <- sum(pfitf$wResiduals[[kf]]^2)
  SSY <- sum(qr.resid(pfitf$wQRs[[1]], pfitf$wY)^2)
  SSm <- SSEr - SSEf
  Fs <- (SSm/(dfr-dff))/(SSEf/dff)
  Qrf <- pfitf$wQRs[[kf]]
  Qrr <- pfitr$wQRs[[kr]]
  ind <- perm.index(n, iter, seed=seed)
  if(!is.null(groups) && !is.null(slope)) pairwise.cond <- "slopes"
  if(!is.null(groups) && is.null(slope)) pairwise.cond <-"means"
  if(is.null(groups) && is.null(slope)) pairwise.cond <-"none"
  if(is.null(groups) && !is.null(slope)) {
    print("No groups for which to compare means or slopes.  No pairwise tests will be performed")
    pairwise.cond <-"none"
  }
  data.types <- lapply(data, class)
  keep = sapply(data.types, function(x) x != "array" & x != "phylo")
  if(!is.null(data)) dat2 <- as.data.frame(data[keep]) else dat2 <- pfitf$data
    
  if(!is.null(groups)){
    g.match <- match(names(dat2), attr(terms(groups), "term.labels"))
    if(!all(is.na(g.match))) gps <- dat2[,which(!is.na(g.match))] else
      gps <- model.frame(groups)
    } else gps <- NULL
  if(!is.null(slope)){
    s.match <- match(names(dat2), attr(terms(slope), "term.labels"))
    if(!all(is.na(s.match))) slp <- dat2[,which(!is.na(s.match))] else
      slp <- model.frame(slope)
    } else slp <- NULL

  if(pairwise.cond == "none"){
    P <- sapply(1:(iter+1), function(j){
      y <- (pfitr$residuals[[kr]][ind[[j]],] + pfitr$fitted[[kr]])*w
      sum(qr.resid(Qrr,y*w)^2)-sum(qr.resid(Qrf,y*w)^2)
    })
    P.val <- pval(P) 
    Z.score <- effect.size(P)
  }
  
  if(pairwise.cond == "means") {
    P <- lsms <- NULL
    jj <- iter+1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    while(jj > 0){
      ind.j <- ind[j]
      Yr <- lapply(1:length(j), function(j) (pfitr$residuals[[kr]][ind.j[[j]],] + pfitr$fitted[[kr]])*w)
      P <- c(P, lapply(1:length(j), function(j) sum(qr.resid(Qrr,Yr[[j]])^2)-sum(qr.resid(Qrf,Yr[[j]])^2)))
      lsms <- c(lsms, apply.ls.means(pfitf, Yr, g = gps, data = dat2)) 
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
    }
      P <- simplify2array(P)
      P.dist <- lapply(1:length(lsms), function(j){as.matrix(dist(lsms[[j]]))})
      P.val <- pval(P) 
      Z.score <- effect.size(P)
      Means.dist <- P.dist[[1]]
      P.dist.s <- simplify2array(P.dist)
      P.Means.dist <- Pval.matrix(P.dist.s)
      Z.Means.dist <- Effect.size.matrix(P.dist.s)
  }
  
  if(pairwise.cond == "slopes") {
    P <- g.slopes <- NULL
    jj <- iter+1
    if(jj > 100) j <- 1:100 else j <- 1:jj
    while(jj > 0){
      ind.j <- ind[j]
      Yr <- lapply(1:length(j), function(j) (pfitr$residuals[[kr]][ind.j[[j]],] + pfitr$fitted[[kr]])*w)
      P <- c(P, lapply(1:length(j), function(j) sum(qr.resid(Qrr,Yr[[j]])^2)-sum(qr.resid(Qrf,Yr[[j]])^2)))
      g.slopes <- c(g.slopes, apply.slopes(pfitf, g=gps, slope=slp, Yr, data=dat2)) 
      jj <- jj-length(j)
      if(jj > 100) kk <- 1:100 else kk <- 1:jj
      j <- j[length(j)] +kk
    }
    P <- simplify2array(P)
    P.slopes.dist <- Map(function(y) as.matrix(dist(y)), g.slopes) 
    P.cor <- Map(function(y) vec.cor.matrix(y), g.slopes) 
    P.val <- pval(P) 
    Z.score <- effect.size(P)
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
  }
  
  if(pairwise.cond == "none"){
    out <- list(anova.table = anova.table, 
                coefficients=pfitf$coefficients,
                Y=pfitf$Y, X=pfitf$X,
                QR = pfitf$QRs[[kf]], fitted=pfitf$fitted[[kf]],
                residuals = pfitf$residuals[[kf]],
                weights = w, data = dat2, random.SS = P,
                Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
                call= match.call()
                )
  }
  if(pairwise.cond == "means"){
    out <- list(anova.table = anova.table, LS.means = lsms[[1]], 
               LS.means.dist = Means.dist, Z.means.dist = Z.Means.dist, P.means.dist = P.Means.dist, 
               coefficients=pfitf$coefficients, 
               Y=pfitf$Y, X=pfitf$X, 
               QR = pfitf$QR[[kf]], fitted=pfitf$fitted[[kf]], 
               residuals = pfitf$residuals[[kf]], 
               weights = w, data = dat2, random.SS = P, random.means.dist = P.dist,
               Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
               call= match.call()
              )
  }
  if(pairwise.cond == "slopes"){
    if(angle.type == "r"){
      out <- list(anova.table = anova.table, slopes = g.slopes[[1]],
      slopes.dist = P.slopes.dist[[1]], P.slopes.dist = P.val.slopes.dist,
      Z.slopes.dist = Z.slopes.dist,
      slopes.cor = P.cor[[1]], P.slopes.cor = P.val.cor, Z.slopes.cor = Z.cor,
      random.slopes = g.slopes, random.slopes.dist = P.slopes.dist, random.slopes.cor = P.cor,
      coefficients=pfitf$coefficients, 
      Y=pfitf$Y, X=pfitf$X, 
      QR = pfitf$QR[[kf]], fitted=pfitf$fitted[[kf]], 
      residuals = pfitf$residuals[[kf]], 
      weights = w, data = dat2, random.SS = P, 
      Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
      call= match.call()
      )
    } else {
      out <- list(anova.table = anova.table, slopes = g.slopes[[1]],
      slopes.dist = P.slopes.dist[[1]], P.slopes.dist = P.val.slopes.dist,
      Z.slopes.dist = Z.slopes.dist,
      slopes.angles = angles.obs, P.angles = P.val.cor, Z.angles = Z.cor,
      random.slopes = g.slopes, random.slopes.dist = P.slopes.dist, random.angles = random.angles,
      coefficients=pfitf$coefficients, 
      Y=pfitf$Y, X=pfitf$X, 
      QR = pfitf$QR[[kf]], fitted=pfitf$fitted[[kf]], 
      residuals = pfitf$residuals[[kf]], 
      weights = w, data = dat2, random.SS = P, 
      Terms = pfitf$Terms, term.labels = pfitf$term.labels, permutations = iter+1,
      call= match.call()
      )
    }
  }
class(out) <- "advanced.procD.lm"
out
}
  