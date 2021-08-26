#' Eigen-analysis for different covariance matrices defined by a modular hypothesis
#'
#' Function performs eigen analysis on covariance matrices for modularity only, integration only, and
#' both integration and modularity (total), as a descriptive way to understand whether modularity or integration are more
#' prominent in a data set, for a particular modular hypothesis.
#'
#' The function performs eigen-analysis on a covariance matrix (of shape data) four different ways: the observed covariance
#' matrix, a covariance matrix that has only within-module covariances (modularity) , a covariance matrix that has only 
#' among-module covariances (integration), and a diagonal matrix, matching the variable variances found in each other 
#' covariance matrix (independence).  The trace of each of these matrices is the same, meaning the sum of eigenvalues is also
#' the same among matrices.  However, they will differ in the distribution of eigenvalues.  How these distributions differ can provide
#' descriptive support for whether modularity, integration, or both are prominent components of variation in the data.
#' 
#' An eigen-analysis of the diagonal (independence) matrix is akin to a broken-stick model, providing eigen structure
#' for the case that variables are all independent.  One can determine the dimensions for which eigenvalues
#' exceed independent eigenvalues as "relevant dimensions".  Fewer relevant dimensions suggest stronger signal in the data,
#' whether due to modularity, integration, or both.  A similar number of relevant dimensions for the total covariance matrix
#' and either the modularity or integration covariance matrix would suggest that shape diversity is consistent
#' with one of these signals.  Comparison of eignevalues, or their proportion of the trace, by component, is also a valuable way to
#' discern whether variation is more consistent with modularity or integration.  
#' 
#' It is also possible to consider the consistency in eigenvector (principal component) orientation 
#' among the different types of covariance matrices.  Results are provided both for choice of component vector
#' correlations and Krzanowski (1979) correlations (mean of squared vector cross-products), using \code{\link{summary.module.eigen}}.
#' 
#' If a phylogeny is provided, a generalized least-squares (GLS) approach is used to estimate the mean vector (tree root values).
#' One has the option whether to transform residuals to obtain a GLS covariance matrix or simply center data on the GLS mean
#' to estimate a covariance matrix.  The latter (transform = FALSE) centers data on the GLS mean rather than the ordinary least-squares 
#' (OLS) mean but does not condition residuals on the phylogenetic covariances in the estimation of variances and covariances.
#' 
#'   
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param A2 An optional 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables) for a second partition
#' @param partition.gp A list of which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C). Required when only 1 dataset provided.
#' @param phy Optional argument to include a class \code{phylo} phylogenetic tree.  Tip labels must match data names.
#' This argument instructs the function to estimate a phylogenetic covariance matrix based on a Brownian motion model
#' of evolutionary divergence.  The Cov argument allows a user to define a hypothetical covariance matrix, if different 
#' than a BM model.
#' @param Cov Optional argument to include a hypothetical covariance matrix used for non-independence of observations.  
#' Row and column names must match data names.  If both a phy and Cov are provided, Cov will override phy.
#' @param transform A logical argument for whether to use transformed residuals, if a phylogeny is provided.  If TRUE,
#' a GLS covariance matrix will be estimated; if FALSE, data will be centered on GLS mean but an OLS covariance matrix 
#' will be estimated.  The former is more representative of covariances independent of phylogeny;  the latter 
#' is more representative of dispersion in the tangent space.
#' @param only.values A logical argument for whether to only return eigenvalues.  If TRUE, vector correlations
#' using \code{\link{summary.module.eigen}} will not be possible.
#' @param tol A value between 0 and 1 to indicate a tolerance for relevant dimensions (via a broken stick model).
#' This value is the fraction of the standard deviation of the first eigenvector.  It helps to cut off dimensions that are 
#' no different in eigenvalue than independent vectors.  If tol = 0, a strict criterion is used to retain all vectors
#' that have larger eigenvalues than expected from independent variables. 
#' 
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @return Objects of class "model.eigen" return a list of the following:
#'  \item{evals}{The eigenvalues of covariance matrices.}
#'  \item{evecs}{The eigenvectors of all covariance matrices. This is NULL if only values are returned.}
#'  \item{rel.dim}{The relevant dimensions based on a broken stick model.}
#'  \item{n.modules}{The number of modules considered.}
#'  \item{prop.mod.cells}{The proportion of elements (cells) in the covariance matrix corresponding to 
#'  modularity covariances.}
#'  \item{prop.int.cells}{The proportion of elements (cells) in the covariance matrix corresponding to 
#'  integration covariances.}
#' @references  Krzanowski, W. J. 1979. Between-groups comparison of principal components. 
#' Journal of the American Statistical Association, 74, 703–707.
#' @references Collyer et al. In review.
#' @seealso \code{\link{summary.module.eigen}}, \code{\link{plot.module.eigen}}
#' \code{\link{two.b.pls}}, \code{\link{modularity.test}}, 
#' \code{\link{phylo.integration}}, \code{\link{phylo.modularity}},
#' and \code{\link{compare.pls}}
#' @examples
#' 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
#' # OLS approach
#' 
#' ME.ols <- module.eigen(Y.gpa$coords, partition.gp = land.gps,
#' only.values = FALSE)
#' summary(ME.ols)
#' plot(ME.ols)
#' 
#' # GLS-centered approach.  This approach find the GLS mean but does not
#' # transform residuals.  Thus, it is produces OLS-like covariance matrcies
#' # that are centered on the GLS mean rather than the OLS mean
#' 
#' ME.glsc <- module.eigen(Y.gpa$coords, partition.gp = land.gps,
#' only.values = FALSE, phy = plethspecies$phy, transform = FALSE)
#' summary(ME.glsc)
#' plot(ME.glsc)
#' 
#' # GLS-transformed approach.  This approach find the GLS mean and
#' # transforms residuals.  Thus, it produces GLS covariance matrices,
#' # from residuals rendered independent of phylogenetic covariances.
#' 
#' ME.glst <- module.eigen(Y.gpa$coords, partition.gp = land.gps,
#' only.values = TRUE, phy = plethspecies$phy, transform. = TRUE)
#' summary(ME.glst)
#' plot(ME.glst)

module.eigen <- function(A, A2 = NULL, partition.gp = NULL, 
                         phy = NULL, Cov = NULL,
                         transform. = TRUE, only.values = TRUE,
                         tol = 0.001){
  
  if(any(is.na(A)))
    stop("\nData matrix contains missing values. Estimate these first (see 'estimate.missing').",
         call. = FALSE)  
  
  x <- try(two.d.array(A), silent = TRUE)
  if(inherits(x, "try-error")) x <- try(as.matrix(A), silent = TRUE)
  if(inherits(x, "try-error"))
    stop("\nA is not a suitable data array for analysis. ", call. = FALSE)
  
  namesX <- rownames(x)
  if (is.null(namesX)) namesX <- 1:NROW(x)
  
  if(!is.null(A2)) {
    if(any(is.na(A2)))
      stop("\nData matrix 2 contains missing values. Estimate these first (see 'estimate.missing').",
           call. = FALSE) 
    
    y <- try(two.d.array(A2), silent = TRUE)
    if(inherits(y, "try-error")) y <- try(as.matrix(A2), silent = TRUE)
    if(inherits(y, "try-error"))
      stop("\nA2 is not a suitable data array for analysis. ", call. = FALSE)
    
    namesY <- rownames(y)
    if (is.null(namesY)) {
      namesY <- namesX
      cat("\nNo names for A2.  Data are assumed to be ordered as in A.\n")
    }
    
    if(!setequal(namesX, namesY))
      stop("\nSpecimen names do not match between A and A2. ", call. = FALSE)
    
    gps <- factor(c(rep(1, ncol(x)), rep(2, ncol(y))))
    x <- cbind(x, y)
  }
  
  if(!is.null(partition.gp) && is.null(A2)){
    partition.gp <- as.factor(partition.gp)
    if (length(dim(A)) == 3){ 
      dims <- dim(A)
      p <- dims[1]
      k <- dims[2]
      
      if(length(partition.gp) != p) 
        stop("\nNot all landmarks are assigned to a partition.", call. = FALSE)
      
      gps <- as.factor(rep(partition.gp, k, each = k, length = p * k))  
      A1.new <- A[which(partition.gp == levels(partition.gp)[1]),,]
      A2.new <- A[which(partition.gp != levels(partition.gp)[1]),,]
    }
    
    if (length(dim(A)) == 2){ 
      
      if(length(partition.gp) != ncol(x))
        stop("\nNot all variables are assigned to a partition.", call. = FALSE)
      
      gps <- as.factor(partition.gp) 
      A1.new <- A[, which(partition.gp == levels(partition.gp)[1])]
      A2.new <- A[, which(partition.gp != levels(partition.gp)[1])]
      
    }
  }
    
    ngps <- nlevels(gps)
    ind.levels <- levels(gps)
    x <- as.matrix(x)
    n <- nrow(x)
    
    if(!is.null(phy) || !is.null(Cov)) {
      if(!is.null(phy) && is.null(Cov)) {
        Cov <- fast.phy.vcv(phy)
        Pcov <- Cov.proj(Cov, rownames(x))}
      
      if(is.null(phy) && !is.null(Cov)) {
        Pcov <- try(Cov.proj(Cov, rownames(x)), silent = TRUE)
        if(inherits(Pcov, "try-error"))
          stop("The names of Covariance matrix do not seem to match data names.\n",
               call. = FALSE)
      }
      
      if(!is.null(phy) && !is.null(Cov)) {
        cat("Both phy and Cov were provided; only Cov will be used\n")
        Pcov <- try(Cov.proj(Cov, rownames(x)), silent = TRUE)
        if(inherits(Pcov, "try-error"))
          stop("The names of Covariance matrix do not seem to match data names.\n",
               call. = FALSE)
      }
      
      ones <- matrix(1, n)
      B <- lm.fit(Pcov %*% ones, Pcov %*% x)$coefficients
      R <- x - ones %*% B
      V <- if(transform.) crossprod(Pcov %*% R) / (n-1) else
        crossprod(R) / (n-1)
      
    } else V <- var(x)
    
    Ind <- diag(diag(V))
    M <- Ind
    for(i in 1:ngps) {
      keep <- which(gps == ind.levels[i])
      M[keep, keep] <- V[keep, keep]
    }
    Int <- V - M + Ind
    
    mcells <- length(which(M != 0)) - nrow(V)
    icells <- length(which(Int != 0)) - nrow(V)
    mprop <- mcells / (length(V) - nrow(V))
    iprop <- icells / (length(V) - nrow(V))
    
    S <- svd(V)
    Sint <- svd(Int)
    Sm <- svd(M)
    Sind <- svd(Ind)
    
    tol.d <- sqrt(S$d[1])*tol
    
    bsV <- which(sqrt(S$d) > (sqrt(Sind$d) + tol.d))
    bsV <- which(cumsum(bsV) <= cumsum(1:length(bsV)))
    bsM <- which(sqrt(Sm$d) > (sqrt(Sind$d) + tol.d))
    bsM <- which(cumsum(bsM) <= cumsum(1:length(bsM)))
    bsI <- which(sqrt(Sint$d) > (sqrt(Sind$d) + tol.d))
    bsI <- which(cumsum(bsI) <= cumsum(1:length(bsI)))
    
    if(length(bsV) == 0) bsV <- 1
    if(length(bsM) == 0) bsM <- 1
    if(length(bsI) == 0) bsI <- 1
    
    
    out <- list(
      evals = list(total = S$d, mod = Sm$d, int = Sint$d,
                    ind = Sind$d),
      evecs = NULL,
      rel.dims = list(total = bsV, mod = bsM, int = bsI),
      n.modules = ngps,
      prop.mod.cells = mprop,
      prop.int.cells = iprop
    )
  
    if(!only.values)  
      out$evecs <- list(total = S$u, mod = Sm$u, int = Sint$u,
           ind = Sind$u)
    
  class(out) <- "module.eigen"
  out
  
  }
  



  