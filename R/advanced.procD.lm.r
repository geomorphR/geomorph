#' Procrustes ANOVA and pairwise tests for shape data, using complex linear models
#'
#' The function quantifies the relative amount of shape variation explained by  a suite of factors
#' and covariates in a "full" model, after accounting for variation in a "reduced" model. Inputs are 
#' formulae for full and reduced models (order is not important), plus indication if means or slopes 
#' are to be comapred among groups, with appropriate formulae to define how they should be compared.
#' 
#'   The response matrix 'y' can be in the form of a two-dimensional data 
#'   matrix of dimension (n x [p x k]) or a 3D array (p x k x n). It is assumed that the landmarks have previously 
#'   been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The names specified for the 
#'   independent (x) variables in the formula represent one or more 
#'   vectors containing continuous data or factors. It is assumed that the order of the specimens in the 
#'   shape matrix matches the order of values in the independent variables.
#'
#'   The function performs statistical assessment of the terms in the model using Procrustes distances among 
#'   specimens, rather than explained covariance matrices among variables. With this approach, the sum-of-squared 
#'   Procrustes distances are used as a measure of SS (see Goodall 1991). The SS betwen models is evaluated through 
#'   permutation. In morphometrics this approach is known as a Procrustes ANOVA (Goodall 1991), which is equivalent
#'   to distance-based anova designs (Anderson 2001). Unlike \code{\link{procD.lm}}, this function is strictly for comparison
#'   of two nested models. The function will readily accept non-nested models, but the results will not be meaningful.
#'   (Use of \code{\link{procD.lm}} will be more suitable in most cases.)  
#'   A residual randomization permutation procedure (RRPP) is utilized 
#'   for reduced model residuals to evalute the SS between models (Collyer et al. 2015).  Effect-sizes (Z-scores) are 
#'   computed as standard deviates of the SS sampling 
#'   distributions generated, which might be more intuitive for P-values than F-values (see Collyer et al. 2015).  
#'   
#'   Pairwise tests are only performed if formulae are provided to compute such results.
#' @param f1 A formula for a linear model, containing the response matrix (e.g., y ~ x1 + x2)
#' @param f2 A formula for another linear model (e.g., ~ x1 + x2 + x3 + a*b) (f1 and f2 should be nested)
#' @param groups A formula for grouping factors (e.g., ~ a, or ~ a*b)
#' @param slope A formula with one covariate (e.g., ~ x3)
#' @param angle.type A value specifying whether differences between slopes should be represented by vector
#' correlations (r), radians (rad) or degrees (deg)
#' @param iter Number of iterations for significance testing
#' @param verbose A logical value specifying whether additional output should be displayed (see Value below)
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{procD.lm}}
#' @return Function returns an ANOVA table of statistical results for model comparison: error df (for each model), SS, MS,
#' F ratio, Z, and Prand.  The following may also be returned.
#'   \item{Means.dist}{Pairwise distance between means, if applicable}
#'   \item{LS.Means.dist}{Pairwise distance between LS means, if applicable}
#'   \item{Prob.Means.dist}{P-values for pairwise distances between means}
#'   \item{Slopes.dist}{Pairwise distance between slope vectors (difference in amount of shape change), if applicable}
#'   \item{Prob.Slopes.dist}{P-values for pairwise distances between slopes}
#'   \item{Slopes.correlation}{Pairwise vector correlations between slope vectors, if applicable}
#'   \item{Prob.Slopes.cor}{P-values for pairwise correlations between slope vectors (high correlation less significant)}
#'   \item{Slopes.angle}{Angles between between slope vectors, if applicable}
#'   \item{Prob.Slopes.angle}{P-values for pairwise angles between slope vectors}
#'   \item{SS.rand}{Random SS from RRPP permutations (when {verbose=TRUE})}
#'   \item{random.mean.dist}{random pairwise distances between means from RRPP permutations (when {verbose=TRUE})}
#'   \item{random.slope.dist}{random pairwise distances between slopes from RRPP permutations (when {verbose=TRUE})}
#'   \item{random.slope.comp}{random pairwise slope direction comparisons (r or angle) from RRPP permutations (when {verbose=TRUE})}
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 113: doi:10.1038/hdy.2014.75.
#' @examples
#'data(plethodon)
#'Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#'Y<-two.d.array(Y.gpa$coords)
#'CS <- Y.gpa$Csize
#'sp<- plethodon$species
#'st<- plethodon$site
#'
#'# Example of a nested model comparison (as with ANOVA with RRPP)
#'advanced.procD.lm(Y ~ log(CS) + sp, ~ log(CS)*sp*st, iter=19)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons (replaces pairwiseD.test)
#'advanced.procD.lm(Y ~ st*sp, ~st + sp, groups = ~st*sp, iter=19)
#'
#'# Example of a test of a factor interaction, plus pairwise comparisons, 
#'# accounting for a common allomtry  (replaces pairwiseD.test)
#'advanced.procD.lm(Y ~ log(CS) + st*sp, 
#'~log(CS) + st + sp, 
#'groups = ~st*sp, slope = ~log(CS), iter=19)
#'
#'# Example of a test of homogeneity of slopes, plus pairwise slopes comparisons
#'# (replaces pairwise.slope.test)
#'advanced.procD.lm(Y ~ log(CS)*st*sp, 
#'~log(CS) + st*sp, 
#'groups = ~st*sp, slope = ~log(CS), angle.type = "deg", iter=19)
#'
#'# Example of partial pairwise comparisons, given greater model complexity
#'advanced.procD.lm(Y ~ log(CS)*st*sp, 
#'~log(CS) + st*sp, 
#'groups = ~sp, slope = ~log(CS), angle.type = "deg", iter=19)

advanced.procD.lm<-function(f1, f2, groups = NULL, slope = NULL, angle.type = c("r", "deg", "rad"), iter=999, verbose = FALSE){
  data=NULL
  f1 <- as.formula(f1)
  Y <- eval(f1[[2]], parent.frame())
  if(length(dim(Y)) == 3)  Y <- two.d.array(Y) else Y <- as.matrix(Y)
  f1 <- as.formula(paste(c("Y",f1[[3]]),collapse="~"))
  f2 <- as.formula(f2)
  k1 <- length(attr(terms(f1), "term.labels"))
  k2 <- length(attr(terms(f2), "term.labels"))
  if (any(is.na(Y)) == T) stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  if(k1 > k2) ff <- f1 else ff <- f2
  if(k1 > k2) fr <- f2 else fr <- f1
  if(k1 == k2) stop("Models have same df")
  full.terms <- terms(ff)
  red.terms <- terms(fr)
  dfr <- nrow(Y) - qr(model.matrix(red.terms))$rank
  dff <- nrow(Y) - qr(model.matrix(full.terms))$rank 
  SSEr <- SSE(lm(Y ~ model.matrix(red.terms) - 1))
  SSEf <- SSE(lm(Y ~ model.matrix(full.terms) - 1))  
  SSm <- SSEr - SSEf
  Fs <- (SSm/(dfr-dff))/(SSEf/dff)
  
  R <- as.matrix(resid(lm(Y ~ model.matrix(red.terms) - 1)))
  P <- array(,iter+1)
  P[1] <- SSm
  m <-Bslopes <-pairwise.cond <- NULL
  if(is.null(groups)==FALSE && is.null(slope)==FALSE) pairwise.cond <- "slopes"
  if(is.null(groups) == FALSE && is.null(slope)==TRUE) pairwise.cond <-"means"
  if(is.null(groups) && is.null(slope)) pairwise.cond <-"none"
  if(is.null(groups) == TRUE && is.null(slope)==FALSE) {
    print("No groups for which to compare means or slopes.  No pairwise tests will be performed")
    pairwise.cond <-"none"
  }
  
  if(pairwise.cond == "none"){
    for(i in 1:iter){
      Rr <- R[sample(nrow(R)),]
      pseudoY = predict(lm(Y ~ model.matrix(red.terms) - 1)) + Rr
      P[i+1] <- SSE(lm(pseudoY ~ model.matrix(red.terms) - 1)) - SSE(lm(pseudoY ~ model.matrix(full.terms) - 1)) 
    }
    P.val <- pval(P)
    Z.score <- effect.size(P)
    anova.tab <- data.frame(df = c(dfr,dff), SSE = c(SSEr, SSEf), SS = c(NA, SSm),
                            F = c(NA, Fs), Z = c(NA, Z.score), P = c(NA,P.val))
    rownames(anova.tab) <- c(paste(attr(red.terms, "term.labels"), collapse="+"),
                             paste(attr(full.terms, "term.labels"), collapse="+"))
    attr(anova.tab, "heading") <- "\nANOVA with RRPP\n"
    class(anova.tab) <- c("anova", class(anova.tab))
    if(verbose == TRUE) out = list(anova.table = anova.tab, SS.rand = P) else out = anova.tab
  }
  
  if(pairwise.cond == "means") {
    gr <- as.factor(single.factor(groups, keep.order=FALSE))
    m <- ls.means(gr, cov.mf = NULL, Y)
    P.dist <- array(, c(nrow(m), nrow(m),iter+1))
    P.dist[,,1] <- as.matrix(dist(m))
    for(i in 1:iter){
      Rr <- R[sample(nrow(R)),]
      pseudoY = predict(lm(Y ~ model.matrix(red.terms) - 1)) + Rr
      P[i+1] <- SSE(lm(pseudoY ~ model.matrix(red.terms) - 1)) - SSE(lm(pseudoY ~ model.matrix(full.terms) - 1)) 
      mr <- ls.means(gr, cov.mf = NULL, pseudoY)
      P.dist[,,i+1] <- as.matrix(dist(mr))  
    }
    P.val <- pval(P)
    Z.score <- effect.size(P)

    anova.tab <- data.frame(df = c(dfr,dff), SSE = c(SSEr, SSEf), SS = c(NA, SSm),
                            F = c(NA, Fs), Z = c(NA, Z.score), P = c(NA,P.val))
    rownames(anova.tab) <- c(paste(attr(red.terms, "term.labels"), collapse="+"),
                             paste(attr(full.terms, "term.labels"), collapse="+"))
    attr(anova.tab, "heading") <- "\nANOVA with RRPP\n"
    class(anova.tab) <- c("anova", class(anova.tab))
    Means.dist <- as.matrix(dist(m))
    P.Means.dist <- Pval.matrix(P.dist)
    dimnames(P.Means.dist) = dimnames(Means.dist)
    if(verbose == TRUE) {
      out = list(anova.table = anova.tab, Means.dist = Means.dist, Prob.Means.dist = P.Means.dist, SS.rand = P)
    } else out = list(anova.table = anova.tab, Means.dist = Means.dist, Prob.Means.dist = P.Means.dist)
  }
  
  if(pairwise.cond == "slopes") {
    gr <- as.factor(single.factor(groups, keep.order=FALSE))
    cov <- model.frame(as.formula(slope))
    m <- ls.means(gr, cov, Y)
    Bslopes <- slopes(gr, cov, Y)
    P.mean.dist <- P.slope.dist <- P.cor <- array(, c(nrow(m), nrow(m),iter+1))
    P.mean.dist[,,1] <- as.matrix(dist(m))
    P.slope.dist[,,1] <- as.matrix(dist(Bslopes))
    P.cor[,,1] <- 1-vec.cor.matrix(Bslopes)
    for(i in 1: iter){
      Rr <- R[sample(nrow(R)),]
      pseudoY = predict(lm(Y ~ model.matrix(red.terms) - 1)) + Rr
      P[i+1] <- SSE(lm(pseudoY ~ model.matrix(red.terms) - 1)) - SSE(lm(pseudoY ~ model.matrix(full.terms) - 1)) 
      mr <- ls.means(gr, cov, pseudoY)  
      Bslopes.r <- slopes(gr, cov, pseudoY)
      P.mean.dist[,,i+1] <- as.matrix(dist(mr))     
      P.slope.dist[,,i+1] <- as.matrix(dist(Bslopes.r))
      P.cor[,,i+1] <- 1 - vec.cor.matrix(Bslopes.r)     
    }
    P.val <- pval(P)
    Z.score <- effect.size(P)
    anova.tab <- data.frame(df = c(dfr,dff), SSE = c(SSEr, SSEf), SS = c(NA, SSm),
                            F = c(NA, Fs), Z = c(NA, Z.score), P = c(NA,P.val))
    rownames(anova.tab) <- c(paste(attr(red.terms, "term.labels"), collapse="+"),
                             paste(attr(full.terms, "term.labels"), collapse="+"))
    attr(anova.tab, "heading") <- "\nANOVA with RRPP\n"
    class(anova.tab) <- c("anova", class(anova.tab))
    Means.dist <- as.matrix(dist(m))
    Slopes.dist <- as.matrix(dist(Bslopes))
    angle.type = match.arg(angle.type)
    Angles.dist = vec.ang.matrix(Bslopes, angle.type)
    P.Means.dist <- Pval.matrix(P.mean.dist)
    P.Slopes.dist <- Pval.matrix(P.slope.dist)
    P.Cor <- Pval.matrix(P.cor)
    dimnames(Angles.dist) <- dimnames(P.Means.dist) <- dimnames(P.Slopes.dist) <- dimnames(P.Cor) <- dimnames(Means.dist)
    if(verbose==TRUE) out = list(anova.table = anova.tab, LS.means = m, Group.slopes = Bslopes) else out = list(anova.table = anova.tab)
    if(angle.type == "r") {out = c(out, list(LS.Means.dist = Means.dist, Prob.Means.dist = P.Means.dist, Slopes.dist = Slopes.dist,
                  Prob.Slopes.dist = P.Slopes.dist, Slopes.correlation = Angles.dist, Prob.Slopes.cor = P.Cor))
    } else {out = c(out, list(LS.Means.dist = Means.dist, Prob.Means.dist = P.Means.dist, Slopes.dist = Slopes.dist,
                 Prob.Slopes.dist = P.Slopes.dist, Slopes.angle = Angles.dist, Prob.Slopes.angle = P.Cor))
    }
    if(verbose == TRUE) out = c(out, list(SS.rand = P, random.mean.dist=P.mean.dist, 
                                          random.slopes.dist=P.slope.dist, random.slope.comp=P.cor))
  }
  out
}
