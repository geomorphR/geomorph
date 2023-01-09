#' K-modules simulation
#'
#' Function finds many possible permutations of K modules and ranks them in terms of 
#' greatest maximum eigenvalue or summed eigenvalues for relevant vectors,
#' based on covariance matrices for modules (see \code{\link{module.eigen}}).
#'
#' The function simulates a desired number of modular hypotheses for a set of morphometric data, and 
#' ranks them in terms of either their first eigenvalue or sum of a number of eigenvalues, which can
#' be directed or assumed from a module eigen-analysis performed on a covariance matrix 
#' containing only modular covariances (see \code{\link{module.eigen}}).  An a priori modular hypothesis can 
#' be included (same input as in \code{\link{modularity.test}}, \code{\link{phylo.modularity}}, and 
#' \code{\link{module.eigen}}), 
#' and the minimum number of landmarks in a module can be adjusted.  A module.eigen object
#'  \code{\link{module.eigen}} is required for input. 
#' 
#' The function simulates many modular hypotheses and partitions covariance matrices by modules 
#' (see \code{\link{module.eigen}}),
#' performs eigen-analysis, and retains the eigenvalues along with the hypothesis.  All simulated 
#' outcomes (summed eigenvalues for relevant dimensions) are then 
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
#' @param ME a module.eigen object, see \code{\link{module.eigen}}.
#' @param K The number of modules into which the data are partitioned.
#' @param hyp An optional vector for an a priori, modular  hypothesis, explaining
#' which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C). This is the same as partition.gp in, e.g., \code{\link{module.eigen}}.  If provided, it
#' will be one of simulated hypotheses considered.
#' @param eig.no The number of eigenvalues to sum for ranking covariance strength.
#' @param rel.dims The number of dimensions (eigenvectors) to sum for ranking outcomes.  If NULL, the relevant
#' dimensions from the module.eigen results will be used.  If more dimensions than ar epossible are chosen, the
#' number will be truncated.
#' @param min.lmk A numeric value to indicate the minimum number of landmarks (or variables) for a module.  
#' If this value is larger than p/K, for p landmarks in K partitions, it will be adjusted to be the largest integer 
#' equal to or less than p/k.
#' @param nsims The number of simulations to use.
#' @param seed An optional value for seed control of simulations.  The value can be NULL, "random",
#' or an integer.  See  \code{\link{set.seed}}
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen. 
#' This is helpful for long-running analyses.
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @return Objects of class "K.modules" return a list of the following:
#'  \item{eig.sums}{A rank-ordered vector of the sum of relevant eigenvalues.}
#'  \item{modules}{A rank-ordered list of modular hypotheses, consistent with eigenvalues.  Modules are 
#'  categories presented numerically, irrespective of factor levels used for an a priori modular hypothesis.}
#'  \item{hypothesis}{Whether an a priori hypothesis was used.}
#'  \item{hypothesis.rank}{If an a priori hypothesis was used, where it ranks among all simulated.}
#'  \item{mean}{If coordinate data are input, the mean configuration for use in plotting tools.}
#'  \item{VCV}{The (variance-) covariance matrix produced, for use in plotting tools.}
#'  \item{A}{The array of shape coordinates, if available, passed on for plotting purposes.}
#'  \item{eig1.ref}{The first eigenvalue or sum of eignvalues for relevant dimensions, from eigenanalysis of the total coavriance matrix.
#'  This value is used to relative the eigenvalue sums in random simulations.}
#'  \item{rel.dims}{The relevant dimensions used in the analysis.}
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
K.modules <- function(ME, 
                      rel.dims = NULL, 
                      min.lmk = 2, 
                      nsims = 1000,
                      seed = NULL,
                      print.progress = FALSE){
  
  if(!inherits(ME, "module.eigen"))
    stop("\nNot a class module.eigen object.  Please use module.eigen first.\n",
         call. = FALSE)
  
  if(print.progress) {
    cat("Acquiring covariance matrix and", nsims, "modular hypotheses:\n")
    pb <- txtProgressBar(min = 0, max = nsims, initial = 0, style=3)
  }
  
  hyp <- ME$lm.hypothesis
  var.hyp <- ME$var.hypothesis
  k <- if(length(var.hyp) > length(hyp)) length(var.hyp)/ length(hyp) else 1
  K <- nlevels(hyp <- as.factor(hyp))
  if(is.null(rel.dims)) rel.dims <- max(ME$rel.dims$total)
  V <- ME$VCV
  eig1.ref <- La.svd(V, 0, 0)$d
  eig1.ref <- eig1.ref[1:(min(length(eig1.ref), rel.dims))]
  eig.no <- length(eig1.ref)
  eig1.ref <- sum(eig1.ref)
  p <- length(hyp)
  if(min.lmk > p/K) min.lmk <- floor(p/K)
  A <- ME$A
  
  if(!is.numeric(seed)) {
    if(is.null(seed)) seed <-nsims else
      seed <- sample(1:nsims, 1) 
  }
  set.seed(seed)
  
  sims <- lapply(1:nsims, function(j){
    tol <- 0
    while(tol < min.lmk){
      res <- sample(1:K, size = p, replace = TRUE)
      res <- factor(res, levels = 1:K)
      gps <- by(res, res, length)
      gps[is.na(gps)] <- 0
      tol <- min(gps)
    }
    res <- rep(res, each = k)
    if(print.progress) setTxtProgressBar(pb, j)
    res
  })
  
  names(sims) <- paste("sim", 1:nsims, sep = ".")
  sims[[1]] <- as.numeric(var.hyp)
  names(sims)[[1]] <- "hypothesis"
  
  rm(.Random.seed, envir=globalenv())
  attr(sims, "seed") <- seed
  if(print.progress) close(pb)
  
  dims <- dim(V)
  ind.levels <- levels(var.hyp)
  modVCV <- function(VCV, sim.i) {
    M <- matrix(0, dims[1], dims[2])
    for(i in 1:K) {
      keep <- which(sim.i == ind.levels[i])
      M[keep, keep] <- VCV[keep, keep]
    }
    M
  }
  
  get.eig <- function(M) sum(La.svd(M, 0, 0)$d[1:eig.no])
  
  if(print.progress) {
    cat("\nEigen-analysis for", nsims, "covariance matrices:\n")
    pb <- txtProgressBar(min = 0, max = nsims, initial = 0, style=3)
  }
  
  result <- sapply(1:nsims, function(j){
    sim.j <- sims[[j]]
    M <- modVCV(V, sim.j)
    if(print.progress) setTxtProgressBar(pb, j)
    get.eig(M)
  })
  
  if(print.progress) close(pb)
  ranks <- order(result, decreasing = TRUE)
  sims.sorted <- sims[ranks]
  result <- result[ranks]
  names(result) <- names(sims.sorted)
  
  out <- list(eig.sums = result, modules = sims.sorted)
  out$hypothesis <- !is.null(hyp) 
  
  out$hypothesis.rank <- which(names(sims.sorted) == "hypothesis")
  if(length(out$hypothesis.rank) == 0) out$hypothesis.rank <- NULL
  
  if(!is.matrix(A)){
    mn <- mshape(A)
    class(mn) <- "matrix"
  } else mn <- NULL
  out$mean <- mn
  
  out$VCV <- V
  out$A <- A
  out$eig1.ref <- eig1.ref
  out$rel.dims = eig.no
  
  class(out) <- "K.modules"
  out
}
