#' K-modules simulation
#'
#' Function finds many possible permutations of K modules and ranks them in terms of 
#' greatest maximum eigenvalue,
#' based on covariance matrices for modules (see \code{\link{module.eigen}}).
#'
#' The function simulates a desired number of modular hypotheses for a set of morphometric data, and 
#' ranks them in terms of their first eigenvalue, for an eigenanalysis performed on a covariance matrix 
#' containing only modular covariances (see \code{\link{module.eigen}}).  An a priori modular hypothesis can 
#' be included (same input as in \code{\link{modularity.test}}, \code{\link{phylo.modularity}}, and 
#' \code{\link{module.eigen}}), 
#' and the minimum number of landmarks in a module can be adjusted.  Input data should be
#' a p x k x n array of coordinates for p points in k dimensions, for n specimens.  A matrix of data can also
#' be used but then variables rather than landmarks will be assigned to modules (not advised for 
#' coordinate data).  A
#' phylogeny can also be used for phylogenetic generalized least-squares calculation of covariances.  If a class
#' phylo object is used, a covariance matrix based on a Brownian motion (BM) model of evolutionary divergence will
#' be estimated to account for the non-independence among observations.  Alternatively, one can
#' assert the covariance structure based on an alternative model, overriding BM covariance estimation.  All arguments
#' used in \code{\link{module.eigen}} can be utilized in this function. 
#' 
#' The function simulates many modular hypotheses and partitions covariance matrices by modules 
#' (see \code{\link{module.eigen}}),
#' performs eigen-analysis, and retains the first eigenvalue along with the hypothesis.  All simulated 
#' outcomes are then 
#' rank-ordered to ascertain which modular hypotheses produce the strongest modularity.  Note that the 
#' analysis is "blind"
#' in terms of anatomical feasibility; thus, a good solution might identify modules comprising integrative 
#' rather than modular 
#' covariances.  Plotting tools found in \code{\link{plot.K.modules}} should help elucidate if this is the case.  
#' This function 
#' is probably only appealing if used with \code{\link{plot.K.modules}} and \code{\link{summary.K.modules}}.  
#' Examples below
#' should demonstrate how these tools can be used with real data.
#' 
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param K The number of modules into which the data are partitioned.
#' @param hyp An optional vector for an a priori, modular  hypothesis, explaining
#' which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C). This is the same as partition.gp in, e.g., \code{\link{module.eigen}}.  If provided, it
#' will be one of simulated hypotheses considered.
#' @param nsims The number of simulations to use.
#' @param phy Optional argument to include a class \code{phylo} phylogenetic tree.  Tip labels must match data names.
#' This argument instructs the function to estimate a phylogenetic covariance matrix based on a Brownian motion model
#' of evolutionary divergence.  The Cov argument allows a user to define a hypothetical covariance matrix, if different 
#' than a BM model.
#' @param Cov Optional argument to include a hypothetical covariance matrix used for non-independence of observations.  
#' Row and column names must match data names.  If both a phy and Cov are provided, Cov will override phy.
#' @param min.lmk A numeric value to indicate the minimum number of landmarks (or variables) for a module.  
#' If this value is larger than p/K, for p landmarks in K partitions, it will be adjusted to be the largest integer 
#' equal to or less than p/k.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen. 
#' This is helpful for long-running analyses.
#' @param seed An optional value for seed control of simulations.  The value can be NULL, "random",
#' or an integer.  See  \code{\link{set.seed}}
#' @param ... Arguments passed onto \code{\link{module.eigen}}.
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @return Objects of class "K.modules" return a list of the following:
#'  \item{eigs}{A rank-ordered vector of first eigen-values.}
#'  \item{modules}{A rank-ordered list of modular hypotheses, consistent with eigenvalues.  Modules are 
#'  categories presented numerically, irrespective of factor levels used for an a priori modular hypothesis.}
#'  \item{hypothesis}{Whether an a priori hypothesis was used.}
#'  \item{hypothesis.rank}{If an a priori hypothesis was used, where it ranks among all simulated.}
#'  \item{mean}{If coordinate data are input, the mean configuration for use in plotting tools.}
#'  \item{VCV}{The (variance-) covariance matrix produced, for use in plotting tools.}
#' @references Collyer et al. In review.
#' @seealso \code{\link{summary.K.modules}}, \code{\link{plot.K.modules}}
#' \code{\link{module.eigen}}, \code{\link{modularity.test}}, 
#' \code{\link{phylo.integration}}, \code{\link{phylo.modularity}},
#' @examples
#' 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
K.modules <- function(A, K = 2, hyp = NULL, nsims = 1000, 
                      phy = NULL, Cov = NULL, min.lmk = 2, seed = NULL,
                      print.progress = FALSE,
                      ...){
  
  if(print.progress) {
    cat("Acquiring covariance matrix and", nsims, "modular hypotheses:\n")
    pb <- txtProgressBar(min = 0, max = nsims, initial = 0, style=3)
  }
    
  p <- if(is.matrix(A)) dim(A)[2] else dim(A)[1]
  
  L <- list(...)
  vcv.args <- list(A = A, phy = NULL, Cov = NULL, 
                   transform. = TRUE)
  pass <- match(names(L), c("A", "phy", "Cov", "tranform."))
  vcv.args[pass] <- L[pass]
  VCV <- do.call(get.VCV, vcv.args)
  eig1.ref <- La.svd(VCV, 0, 0)$d[1]

  if(!is.numeric(seed)) {
    if(is.null(seed)) seed <-nsims else
      seed <- sample(1:nsims, 1) 
  }
  set.seed(seed)
  
  if(min.lmk > p/K) min.lmk <- floor(p/K)
  
  sims <- lapply(1:nsims, function(j){
    tol <- 0
    while(tol < min.lmk){
      res <- sample(1:K, size = p, replace = TRUE)
      res <- factor(res, levels = 1:K)
      gps <- by(res, res, length)
      gps[is.na(gps)] <- 0
      tol <- min(gps)
    }
    if(print.progress) setTxtProgressBar(pb, j)
    res
  })
  
  names(sims) <- paste("sim", 1:nsims, sep = ".")
  
  rm(.Random.seed, envir=globalenv())
  attr(sims, "seed") <- seed
  if(print.progress) close(pb)
  
  if(!is.null(hyp)) {
    hyp <- as.factor(hyp)
    if(length(hyp) != p)
      stop("The modular hypothesis is not consistent with the number of landmarks.\n",
           call. = FALSE)
    if(nlevels(hyp) != K)
      stop("The number of levels in the modular hypothesis does not equal K.\n",
           call. = FALSE)
    sims[[1]] <- as.numeric(hyp)
    names(sims)[[1]] <- "hypothesis"
  }
  
  modVCV <- function(VCV, sim.i) {
    dims <- dim(VCV)
    M <- matrix(0, dims[1], dims[2])
    for(i in 1:K) {
      keep <- which(sim.i == i)
      M[keep, keep] <- VCV[keep, keep]
    }
    M
  }
  
  get.eig <- function(M) La.svd(M, 0, 0)$d[1]
  
  if(print.progress) {
    cat("\nEigen-analysis for", nsims, "covariance matrices:\n")
    pb <- txtProgressBar(min = 0, max = nsims, initial = 0, style=3)
  }
   
  result <- sapply(1:nsims, function(j){
    sim.j <- sims[[j]]
    M <- modVCV(VCV, sim.j)
    if(print.progress) setTxtProgressBar(pb, j)
    get.eig(M)
  })
  
  if(print.progress) close(pb)
  ranks <- order(result, decreasing = TRUE)
  sims.sorted <- sims[ranks]
  result <- result[ranks]
  names(result) <- names(sims.sorted)
  
  out <- list(eigs = result, modules = sims.sorted)
  out$hypothesis <- !is.null(hyp) 
  
  out$hypothesis.rank <- which(names(sims.sorted) == "hypothesis")
  if(length(out$hypothesis.rank) == 0) out$hypothesis.rank <- NULL
  
  if(!is.matrix(A)){
    mn <- mshape(A)
    class(mn) <- "matrix"
  } else mn <- NULL
  out$mean <- mn
  
  out$VCV <- VCV
  out$A <- A
  out$eig1.ref <- eig1.ref
  
  class(out) <- "K.modules"
  out
}