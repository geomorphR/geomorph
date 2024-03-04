#' Quantify integration in a set of traits
#'
#' Function quantifies the morphological integration in a set of traits
#'
#' The function quantifies the strength of morphological integration in a set of variables. 
#' Here the set of traits are treated as a single unit, and the overall degree of covariation in them
#' is quantified using the relative eigenvalue index: Vrel (Pavlicev et al. 2009). Following 
#' Conaway and Adams (2022), only the non-trivial dimensions of variation are used in the calculation of Vrel. 
#' The measure is then converted to an effect size (Z-score), based on the procedures in Conaway and Adams (2022). 
#' These may be used in subsequent comparisons of the strength of integration across datasets. 
#' Input for the analysis may be a 3D array of Procrustes coordinates, of a matrix of variables. If the observations
#'  are species related by a phylogeny, the phylogeny may also be included. 
#'  
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param phy A phylogenetic tree of class = "phylo" - see \code{\link[ape]{read.tree}} in library ape
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Objects of class "rel.eig" from integration.Vrel return a list of the following:
#'  \item{Re.obs}{The observed relative eigenvalue index (Vrel).}
#'  \item{Z.obs}{The associated Z-score, which represents the effect size of Vrel.}
#'  \item{ZR}{The effect size translated to a positive scale (so that no integration is ZR = 0).}
#'  \item{ZR.var}{The variance of the effect size.}
#' @references  Pavlicev, M., J. M. Cheverud, and G. P. Wagner. 2009. Measuring morphological 
#' integration using eigenvalue variance. Evolutionary Biology 36:157-170.
#' @references Conaway, M.A., and D.C. Adams. 2022. An effect size for comparing the strength of 
#'   morphological integration across studies. Evolution. 76: 2244-2259.

#' @seealso \code{\link{compare.ZVrel}}
#' @examples
#' \dontrun{
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)    #GPA-alignment    
#' integration.Vrel(Y.gpa$coords)
#' }

integration.Vrel <- function(A,phy = NULL){ 
  if(length(dim(A))==3){ 
    A <- two.d.array(A)
  }
  n <- nrow(A)
  if (!is.null(phy)) {
    phy.parts <- phylo.mat(A, phy)
    Ptrans <- phy.parts$D.mat %*% (diag(n) - matrix(1, n) %*% 
                                     crossprod(matrix(1, n), phy.parts$invC)/sum(phy.parts$invC))
    A <- Ptrans %*% A
  }
  eig.obs <- eigen(cov(A))$values
  eig.obs <- eig.obs[which(zapsmall(eig.obs)>0)] 
  p <- length(eig.obs)
  VNull <- (p+2)/((p*(n-1))+2)
  Vtrans <- 2*VNull-1
  ZN <- 0.5*log((1+Vtrans) / (1-Vtrans))
  Re.obs <- var(eig.obs) / (mean(eig.obs)^2*p)
  Re.obs.trans <- ((2*Re.obs)-1)
  if(Re.obs.trans==1){Re.obs.trans=0.999}
  if(Re.obs.trans==-1){Re.obs.trans=-0.999}
  Z.obs <- 0.5*log((1+Re.obs.trans) / (1-Re.obs.trans))
  ZR <- Z.obs + abs(ZN)
  if (ZR < 0) warning("ZR less than zero (likely because p>N). Use Z.obs for interpretation.")
  ZR.var <- 1/(n-3)
  out <- list(Re.obs = Re.obs, Z.obs = Z.obs, ZR = ZR, ZR.var = ZR.var)
  class(out) <- "rel.eig"
  out
}