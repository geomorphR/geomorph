#' Assessing dimensions of phylogenetic signal in Procrustes shape variables
#'
#' Function calculates the degree of phylogenetic signal across dimensions of shape space from 
#' Procrustes shape variables
#'
#' The function estimates the degree of phylogenetic signal across dimensions of shape space, based 
#' on Procrustes shape variables for a given phylogeny. The method is intended to interrogate patterns 
#' of phylogenetic signal that are not uniformly distributed across shape space, but rather are 
#' concentrated in one or a few dimensions. The degree of phylogenetic signal in data is estimated 
#' using a matrix generalization of Blomberg's Kappa (K), and the phylogenetic signal represented by this
#' matrix is decomposed into linear combinations that describe directions of phylogenetic signal within the 
#' dataspace containing decreasing levels of phylogenetic signal (Mitteroecker et al. 2024). Two 
#' summary measures based on the eigenvalues of K are calculated: the determinant of K (detK) and the 
#' trace of K (traceK). Both describe the degree of phylogenetic signal in multivariate data. In addition, 
#' the  multivariate version of the K-statistic (Kmult: Adams 2014) is also provided. All three 
#' statistics are evaluated via permutation, where the shape data is permuted among the tips of the 
#' phylogeny. Effect sizes are also provided.
#' 
#' It should be noted that by default the function uses OLS centering when calculating the phylogenetic 
#' signal matrix, unlike most authors who have used GLS centering. Both OLS and GLS centering were 
#' proposed in Blomberg et al.s (2003) original exploration of phylogenetic signal. However, for both 
#' mathematical and computational reasons OLS is used here (see justification in Mitteroecker et al. 
#' 2024). Both measures are highly rank correlated with one another (Mitteroecker et al. 2024). Additionally,
#' using using BLOMBERG = TRUE will result in GLS centering, in which case the statistic Kmult 
#' obtained with this function will be identical to that obtained when using \code{\link{physignal}}.
#' 
#' Importantly, because detK and traceK are based on the covariance matrix K, singularity can become
#' an issue. In particular, geometric morphometric shape data will not be of full rank, as several 
#' dimensions are standardized during the Generalized Procrustes Analysis (one may also have fewer 
#' observations than variables, which will also generate redundancies). For this reason, a 
#' principal components analysis of the data is performed, and the redundant dimensions are 
#' removed so that detK and traceK may be computed (see Mitteroecker et al. 2024). Additionally, if 
#' n< (p X k), the last nontrivial PC dimension is also removed, as in this case, using 100% of the 
#' variation results in invariant K-statistics across permutations. 
#'   
#' The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{physignal.eigen}}.
#'   
#' @param Y A matrix (n x [p x k]) or 3D array (p x k x n) containing Procrustes shape variables for a 
#' set of specimens (it is assumed that the data have been subjected to a Generalized Procrustes Analysis)
#' @param phy A phylogenetic tree of class = "phylo" - see \code{\link[ape]{read.tree}} in library ape
#' @param Cov An object covariance matrix describing the phylogenetic relationships among species, which may be 
#' used in place of `tree
#' @param Blomberg A logical value to indicate whether GLS-centering or OLS-centering (the default) should be used
#' @param unit.tree A logical value to indicate whether the tree should be scaled to unit height 
#' @param lambda An optional value for scaling the tree by its phylogenetic signal (as represented by the 
#' parameter lambda)
#' @param test A logical value to indicate whether significance testing should be performed 
#' @param iter Number of iterations for significance testing
#' @param tol A value indicating the magnitude below which 
#' components should be omitted, following projection.  See \code{\link[RRPP]{ordinate}} 
#' for details.
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @keywords analysis
#' @author Dean Adams and Michael Collyer
#' @seealso \code{\link{gm.prcomp}}, \code{\link{physignal}}
#' @export
#' @return Function returns a list with the following components: 
#'   \item{eib.obs}{The observed eigenvalues of the phylogenetic signal matrix, K.}
#'   \item{rand.eigen.values}{The set of eigenvalues from the permuted datasets.}
#'   \item{traceK.obs}{The observed traceK statistic.}
#'   \item{traceK}{The set of traceK statistics from the permuted datasets.}
#'   \item{p.traceK}{The p-value of traceK from permutation.}
#'   \item{Z.traceK}{The effect size of traceK.}
#'   \item{detK.obs}{The observed detK statistic.}
#'   \item{detK}{The set of detK statistics from the permuted datasets.}
#'   \item{p.detK}{The p-value of detK from permutation.}
#'   \item{Z.detK}{The effect size of detK.}
#'   \item{Kmult.obs}{The observed Kmult statistic.}
#'   \item{Kmult}{The set of Kmult statistics from the permuted datasets.}
#'   \item{p.Kmult}{The p-value of multK from permutation.}
#'   \item{Z.Kmult}{The effect size of Kmult.}   
#'   \item{permutations}{The number of random permutations used in the resampling procedure.}
#'   \item{call}{The matched call}
#'   
#' @references Mitteroecker, P., M.L. Collyer, and D.C. Adams. 2024. Exploring phylogenetic signal in 
#' multivariate phenotypes by maximizing Blomberg's K. Systematic Biology. (In Press).
#' @references Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative 
#' data: behavioral traits are more labile. Evolution, 57:717-745.
#' @references Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and 
#' other high-dimensional multivariate data. Systematic Biology.  63:685-697.
#' @examples
#' \dontrun{
#' data(plethspecies) 
#' Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment    
#'
#' #Test for phylogenetic signal in shape
#' PSe.shape <- physignal.eigen(Y = Y.gpa$coords, phy = plethspecies$phy)
#' summary(PSe.shape)
#' plot(PSe.shape)
#' plot(PSe.shape, type = "vectors")
#' 
#' }
physignal.eigen <- function(Y, phy = NULL, Cov = NULL,
                   Blomberg = FALSE,
                   unit.tree = TRUE,
                   lambda = 1, test = TRUE,
                   iter = 999, seed = NULL,
                   tol = 0.001){

  if(tol < 0.001) {tol = 0.001}
  
  if (length(dim(Y)) == 3){ 
    if(is.null(dimnames(Y)[[3]]))
      stop("Data array does not include taxa names as dimnames for 3rd dimension.\n", call. = FALSE)  
    Y <- two.d.array(Y)
    
  }
  if (length(dim(Y))==2){ 
    if(is.null(rownames(Y))) stop("Data matrix does not include taxa names as dimnames for rows.\n", 
                                  call. = FALSE)  
  }
  Y <- center(as.matrix(Y))
  n <- NROW(Y)
  p <- ncol(Y)
  PCA <- ordinate(Y, tol = tol)
  Y <- PCA$x
  if(n<p && Blomberg == FALSE){Y <- Y[,-ncol(Y)]}

  if(is.null(phy) && is.null(Cov))
    stop("Either a tree or covariance matrix is needed.\n",
         call. = FALSE)
  if(is.null(Cov))
    Cov <- fast.phy.vcv(phy)
  Cov.nms <- rownames(Cov)
  if(unit.tree) {
    if(length(unique(diag(Cov))) == 1)
      Cov <- Cov <- Cov / Cov[1,1] else {
        D <- diag(1 / sqrt(diag(Cov)))
        Cov <- D %*% Cov %*% D
        D <- NULL
      }
    dimnames(Cov) <- list(Cov.nms, Cov.nms)
  }
  
  if(lambda != 1) {
    if(lambda < 0 || lambda > 1)
      stop("lambda must be between 0 and 1\n",
           call. = FALSE)
    D <- diag(diag(Cov))
    Cov  <- lambda * (Cov - D) + D
  }
  
  Pcov <- Cov.proj(Cov, id = rownames(Y))
  TY <- as.matrix(Pcov %*% Y)
  X <- matrix(1, nrow(Y))
  TX <- as.matrix(Pcov %*% X)
  fit <- lm.fit(TX, TY)
  R <- if(Blomberg) Y - X %*% fit$coefficients else Y
  MSE0 <- crossprod(R)
  MSE <- crossprod(fit$residuals)
  K <- fast.solve(MSE) %*% MSE0
  eig.ob <- eigen(K)
  
  if(test){
    ind <- perm.index(n = n, iter = iter, block = NULL, seed = seed)
    eigs <- lapply(1:length(ind), function(j){
      Yj <- as.matrix(as.matrix(Y)[ind[[j]], ])
      rownames(Yj) <- rownames(Pcov)
      TYj <- as.matrix(Pcov %*% Yj)
      fitj <- lm.fit(TX, TYj)
      if(Blomberg){
        Rj <- Yj - X %*% fitj$coefficients
        MSE0j <- crossprod(Rj)
      } else MSE0j <- MSE0
      MSEj <- crossprod(fitj$residuals)
      Kj <- fast.solve(MSEj) %*% MSE0j
      eigj<- eigen(Kj, only.values = TRUE)$values
      eigj
    })
    
    minp <- min(sapply(eigs, function(x) length(x)))
    eigs <- sapply(eigs, function(x) x[1:minp])
    
  } else eigs <- NULL
  
    if(!Blomberg){
      if(!is.null(eigs)) {
        Kmult <- sapply(1:length(ind), function(j){
          Yj <- as.matrix(as.matrix(Y)[ind[[j]], ])
          rownames(Yj) <- rownames(Y)
          TYj <- as.matrix(Pcov %*% Yj)
          fitj <- lm.fit(TX, TYj)
          MSE0j <- MSE0
          MSEj <- crossprod(fitj$residuals)
          sum(MSE0j^2) / sum(MSEj^2)
        })
      } else Kmult <- sum(MSE0^2) / sum(MSE^2)
      
    } else{
        if(is.null(phy))
          stop("An input phylogeny is required for Blomberg = TRUE.\n",
               call. = FALSE)
      if(test){iter = iter} else iter = 0
        Kmult <- physignal(A = Y, phy = phy, iter = iter)$random.K
      } 
 
  if(!is.null(eigs)){
    traceK <- apply(eigs, 2, sum)
    detK <- apply(eigs, 2, prod)
  } else {
    
    traceK <- sum(eig.ob$values)
    detK <- prod(eig.ob$values)
    
    }
  
    ptrace <- pval(traceK)
    pdet <- pval(detK)
    pKm <- pval(Kmult)
    ztrace <- effect.size(traceK)
    zdet <- effect.size(detK)
    zKm <- effect.size(Kmult)

  out <- list(eig.obs = eig.ob, rand.eigen.values = eigs,
              traceK.obs = traceK[1], traceK = traceK, 
              detK.obs = detK[1], detK = detK,
              Kmult.obs = Kmult[1], Kmult = Kmult,
              p.traceK = ptrace, p.detK = pdet, p.Kmult = pKm,
              Z.traceK = ztrace, Z.detK = zdet, Z.Kmult = zKm,
              permutations = iter+1)
  class(out) <- "physignal.eigen"
  out

}
