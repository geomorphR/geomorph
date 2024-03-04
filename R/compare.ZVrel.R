#' Comparisons of Effect Sizes from Overall Integration Analyses
#'
#' Function performs an analysis to compare the effect sizes of two or more Vrel effects
#'
#' The function statistically compares the effect sizes of two or more Vrel analyses.  Typically, this
#' function might be used to compare the strength of integration in one dataset as compared with another 
#' (see Conaway and Adams 2022). 
#' 
#' The analysis performs two-sample z-tests based on effect sizes (Z-scores) of Vrel. The method 
#' follows that of Conaway Adams (2022) used to compare the strength of integration across datasets.
#'  
#' To use this function, simply perform \code{\link{integration.Vrel}} on as many samples or as desired.  
#' Any number of objects of class rel.eig can be input. One may perform the comparison as either a 
#' one-tailed or a two-tailed (default) test.
#' 
#' @param ... saved analyses of class rel.eig
#' @param two.tailed A logical value to indicate whether a two-tailed test (typical and default) should be performed.
#' @keywords analysis
#' @export
#' @author Dean Adams
#' @return An object of class compare.rel.eig, returns a list of the following
#' 
#' \item{sample.Re.obs}{A vector of observed Vrel for each sample.}
#' \item{sample.Z.obs}{A vector of effect sizes for each sample.}
#' \item{sample.Z.var}{A vector of variances for each effect size.}
#' \item{pairwise.z}{A matrix of pairwise, two-sample z scores between all pairs of effect sizes.}
#' \item{pairwise.p}{A matrix of corresponding P-values.}
#' @references Conaway, M.A., and D.C. Adams. 2022. An effect size for comparing the strength of 
#'   morphological integration across studies. Evolution. 76: 2244-2259.
#' @examples
#' \dontrun{
#'  data("plethodon")
#'  Y.gpa <- gpagen(plethodon$land)
#'  
#'  coords.gp <- coords.subset(Y.gpa$coords, plethodon$species)
#'  Vrel.gp <- Map(function(x) integration.Vrel(x), coords.gp) 
#'  
#'  out <- compare.ZVrel(Vrel.gp$Jord, Vrel.gp$Teyah)
#'  
#'  summary(out)
#'  }
compare.ZVrel <- function(...,two.tailed = TRUE){
  dots <- list(...)
  tails <- if(two.tailed) 2 else 1
  if(length(dots) < 2) stop("At least two objects of class rel.eig are needed")
  is.rel.eig <- function(x) class(x) == "rel.eig"
  list.check <- sapply(1:length(dots), function(j) any(is.rel.eig(dots[[j]])))
  if(any(list.check == FALSE)) stop("Not all objects are class rel.eig")
  k <- length(list.check)
  list.names <- as.list(substitute(list(...)))[-1L]
  k.combn <- combn(k,2)
  list.re.obs <- sapply(1:k, function(j) dots[[j]]$Re.obs)
  list.ZRs <- sapply(1:k, function(j) dots[[j]]$Z.obs)
  list.vars <- sapply(1:k, function(j) dots[[j]]$ZR.var) 
  z12 <- sapply(1:ncol(k.combn), function(j){
    a <- k.combn[1,j]; b <- k.combn[2,j]
    r1 <- list.ZRs[a]; r2 <- list.ZRs[b]; var1 <- list.vars[a]; var2 <- list.vars[b]
    abs(r1-r2)/sqrt( var1+var2)
  })
  z12.p <- sapply(1:length(z12), function(j) pnorm(abs(z12[[j]]), lower.tail = FALSE) * tails)
  d <- rep(0,k); names(d) <- list.names
  D <-dist(d)
  z12.pw <- p12.pw <- D
  for(i in 1:length(z12)) z12.pw[i] <-z12[i]
  for(i in 1:length(z12)) p12.pw[i] <-z12.p[i]
  names(list.ZRs) <-list.names
  pairwise.z <- as.matrix(z12.pw)
  pairwise.P <- as.matrix(p12.pw)
  diag(pairwise.P) <- 1
  out <- list(sample.Re.obs = list.re.obs, sample.Z.obs = list.ZRs, 
              sample.Z.var = list.vars,
              pairwise.z = pairwise.z,
              pairwise.P = pairwise.P)
  class(out) <- "compare.ZVrel"
  out
}