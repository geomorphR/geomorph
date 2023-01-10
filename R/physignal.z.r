#' Assessing phylogenetic signal effect size in Procrustes shape variables
#'
#' Function calculates the degree of phylogenetic signal from Procrustes shape variables, using
#' effect size estimated from likelihood.
#'
#' The function estimates the degree of phylogenetic signal present in Procrustes shape variables for a given
#' phylogeny, based on the standardized location of the log-likelihood in a distribution generated from
#' randomization of residuals in a permutation procedure (RRPP). Unlike the analysis performed in
#' \code{\link{physignal}}, this function finds the optimal branch-length transformation, lambda, which maximizes 
#' likelihood.  This step results in better statistical properties (Collyer et al. 2022), and the effect
#' size (Z) -- the standardized location of the observed log-likelihood in its sampling distribution -- has a linear 
#' relationship with lambda.  Although, Pagel's lambda (Pagel et al. 1999) and Blomberg's K (Blomberg et al. 2003;
#' Adams 2014) are frequently used as measures of the amount of phylogenetic signal (and both are also reported), 
#' the Z-score is an effect size that can be compared among different traits 
#' (see \code{\link{compare.physignal.z}}), and it more precisely describes the fit of data to a tree 
#' (Collyer et al. 2022).
#' 
#' It is assumed that the landmarks have previously been aligned 
#' using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].  Multivariate data are subjected
#' to a phylogenetically aligned component analysis to assure that data dimensions do not exceed the number of 
#' observations, making estimation of multivariate log-likelihoods possible.  The tolerance argument can be used to assure
#' an appropriate number of data dimensions are used.  The number of phylogenetically aligned components (PAC.no) 
#' can also be used to estimate log-likelihood in fewer dimensions.  This might be necessary for some data sets,
#' as residual covariance matrices can be ill-conditioned, especially if the number of variables exceeds the number
#' of observations.

#' It is assumed that lambda is constant for all shape variables.  It would be possible to consider, for example,
#' whether phylogenetic signal is consistent among different modules via separate analyses, followed by
#' analysis with \code{\link{compare.physignal.z}}.  This could be facilitated by using 
#' \code{\link{coords.subset}} to partition shape data. 
#' This function can also be used with univariate data (i.e. centroid size) if imported as matrix with rownames
#' giving the taxa names. 
#' 
#' Optimization of the tree branch scaling parameter, lambda, is not a process with a single, resolved solution.  For
#' univariate data, optimization is fairly simple: lambda is optimized at the value that yields maximum likelihood.  For
#' multivariate data, generalization is not straightforward.  Although this function currently assumes all variables 
#' in multivariate data have the same lambda, maximizing likelihood over all variables can ignore strong phylogenetic signal
#' restricted to just some of the variables.  Therefore, optimization can can consider the latent strength of phylogenetic 
#' signal in different ways.  Below is a summary of optimization methods, including advantages or disadvantages that might be
#' incurred.  Future versions of this function might update optimization methods, as more research affirms better methods.
#' 
#' \itemize{
#' \item{\bold{burn}}{  A burn-in optimization samples a spectrum of lambda from 0 to 1, finding the effect
#' size, Z, from several distributions of log-likelihoods, in order to find the lambda value that maximizes the
#' effect size rather than the log-likelihood.  Once this value of lambda is obtained, the requested number
#' of permutations are run and the final effect size is estimated from the distribution of random log-likelihoods
#' generated.  This is perhaps the best optimization method to assure that the result effect size is 
#' maximized, but it requires considerably more time than the other methods.}
#' \item{\bold{mean}}{  Lambda is optimized for every variable and the mean of the lambdas is used as the 
#' optimized lambda for estimating multivariate log-likelihoods.  This method will guard against a tendency
#' for multivariate optimization to be biased toward 0 or 1, but will also tend to find intermediate values
#' of lambda, even if there is strong phylogenetic signal for some variables.}
#' \item{\bold{front}}{  This method performs phylogenetically aligned component analysis (PACA) and
#' determines the relevant number of components for which covariation between data and phylogeny is 
#' possible (as many or fewer than the number of variables).  PACA front-loads phylogenetic signal in first 
#' components. Thus, lambda is optimized with multivariate likelihood
#' but only for a portion of the full data space where phylogenetic signal might be found.  This method 
#' will tend to bias lambda toward 1, making analysis more similar to using multivariate K, as in
#' \code{\link{physignal}}.}
#' \item{\bold{all}}{  This method simply optimizes lambda using a multivariate likelihood for optimization.  
#' This is the simplest generalization of univariate lambda optimization, but it will tend to optimize lambda
#' at 0 or 1.}
#' \item{\bold{numeric override}}{  Users can override optimization, by choosing a fixed value.  Choosing
#' lambda = 1 will be fundamentally the same as performing the analysis on multivariate K, as in
#' \code{\link{physignal}}.}
#' 
#' }
#' 
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{physignal}}.
#'  The generic function, \code{\link{plot}}, produces a histogram of random log-likelihoods, 
#'  associated with the resampling procedure.
#' 
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing Procrustes shape variables for a set of specimens
#' @param lambda An indication for how lambda should be optimized.  This can be a numeric value between 0 and 1 to override 
#' optimization.  Alternatively, it can be one of four methods: burn-in ("burn"); mean lambda across all data dimensions
#' ("mean"); use a subset of all data dimensions where phylogenetic signal is possible, based on a front-loading
#' of phylogenetic signal in the lowest phylogenetically aligned 
#' components ("front"); or use of all data dimensions to calculate the log-likelihood in the optimization process ("all").  See Details
#' for more explanation.
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param tol A value indicating the magnitude below which phylogenetically aligned components should 
#' be omitted. (Components are omitted if their standard deviations are less than or equal to tol times 
#' the standard deviation of the first component.)  See \code{\link[RRPP]{ordinate}} for more details.
#' @param PAC.no The number of phylogenetically aligned components (PAC) on which to project the data.  Users can
#' choose fewer PACs than would be found, naturally.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @param verbose Whether to include verbose output, including phylogenetically-aligned
#' component analysis (PACA), plus K, lambda,
#' and the log-likelihood profiles across PAC dimensions.
#' @keywords analysis
#' @author Michael Collyer
#' @seealso \code{\link{gm.prcomp}}, \code{\link{physignal}}
#' @export
#' @return Function returns a list with the following components: 
#'   \item{Z}{Phylogenetic signal effect size.}
#'   \item{pvalue}{The significance level of the observed signal.}
#'   \item{random.detR}{The determinants of the rates matrix for all random permutations.}
#'   \item{random.logL}{The log-likelihoods for all random permutations.}
#'   \item{lambda}{The optimized lambda value.}
#'   \item{K}{The multivariate K value.}  
#'   \item{PACA}{A phylogenetically aligned component analysis, based on OLS residuals.}
#'   \item{K.by.p}{The multivariate K value in 1, 1:2, 1:3, ..., 1:p dimensions, for the 
#'   p components from PACA.}   
#'   \item{lambda.by.p}{The optimized lambda in 1, 1:2, 1:3, ..., 1:p dimensions, for the 
#'   p components from PACA.}   
#'   \item{logL.by.p}{The log-likelihood in 1, 1:2, 1:3, ..., 1:p dimensions, for the 
#'   p components from PACA.}
#'   \item{permutations}{The number of random permutations used in the resampling procedure.}
#'   \item{opt.method}{Method of optimization used for lambda.}
#'   \item{opt.dim}{Number of data dimensions for which optimization was performed.}
#'   \item{call}{The matched call}
#'   
#' @references Collyer,  M.L., E.K. Baken, & D.C. Adams.  2022. A standardized effect size for evaluating
#' and comparing the strength of phylogenetic signal. Methods in Ecology and Evolution. 13:367-382.
#' @references Blomberg SP, Garland T, Ives AR. 2003. Testing for phylogenetic signal in comparative 
#' data: behavioral traits are more labile. Evolution, 57:717-745.
#' @references Adams, D.C. 2014. A generalized K statistic for estimating phylogenetic signal from shape and 
#' other high-dimensional multivariate data. Systematic Biology.  63:685-697.
#' @references Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
#' alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
#' Evolution. 73:2352-2367.
#' @references Pagel, M. D. (1999). Inferring the historical patterns of biological evolution. Nature. 401:877-884.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' # Test for phylogenetic signal in shape
#' PS.shape <- physignal.z(A = Y.gpa$coords, phy = plethspecies$phy, 
#' lambda = "front", iter=999)
#' summary(PS.shape)
#' 
#' # Problem with ill-conditioned residual covariance matrix; try shaving one dimension
#' 
#' PS.shape <- physignal.z(A = Y.gpa$coords, phy = plethspecies$phy, 
#' lambda = "front", PAC.no = 7, iter=999)
#' summary(PS.shape)
#' plot(PS.shape)
#' plot(PS.shape$PACA, phylo = TRUE)
#' 
#' # Test for phylogenetic signal in size
#' PS.size <- physignal.z(A = Y.gpa$Csize, phy = plethspecies$phy, 
#' lambda = "front", iter=999)
#' summary(PS.size)
#' plot(PS.size)
physignal.z <- function(A, phy, lambda = c("burn", "mean", "front", "all"), iter = 999, seed = NULL, 
                        tol = 1e-4, PAC.no = NULL, print.progress = FALSE,
                        verbose = FALSE){
  if(any(is.na(A)))
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').\n",
         call. = FALSE)  
  if (length(dim(A)) == 3){ 
    if(is.null(dimnames(A)[[3]]))
      stop("Data array does not include taxa names as dimnames for 3rd dimension.\n", call. = FALSE)  
    Y <- two.d.array(A)
  }
  
  
  if (length(dim(A))==2){ 
    if(is.null(rownames(A))) stop("Data matrix does not include taxa names as dimnames for rows.\n", 
                                  call. = FALSE)  
    Y <- A 
  }
  
  if (is.vector(A)){ 
    if(is.null(names(A))) stop("Data vector does not include taxa names as names.\n", call. = FALSE)  
    Y <- as.matrix(A) 
  }
  
  if (!inherits(phy, "phylo"))
    stop("tree must be of class 'phylo.'\n", call. = FALSE)
  
  N <- length(phy$tip.label)
  if(N != dim(Y)[1]) 
    stop("Number of taxa in data matrix and tree are not not equal.\n", call. = FALSE)  
  
  if(length(match(rownames(Y), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.\n", call. = FALSE)
  
  if(length(match(phy$tip.label,rownames(Y))) != N) 
    stop("Tree missing some taxa in the data matrix.\n", call. = FALSE)
  
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(Y)))))) 
    stop("Names do not match between tree and data matrix.\n", call. = FALSE) 
  
  if(is.null(dim(Y))) Y <- matrix(Y, dimnames = list(names(Y))) 
  
  Y <- as.matrix(Y)
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  
  Cov <- fast.phy.vcv(phy)
  Pcov <- Cov.proj(Cov, rownames(Y))
  
  PaCA <- ordinate(Y, A = Cov, tol = tol,
                   newdata = anc.BM(phy, Y))
  Y <- PaCA$x
  p <- NCOL(Y)
  if(!is.null(PAC.no) && PAC.no < p) {
    Y <- as.matrix(Y)[, 1:PAC.no]
  }
  
  class(PaCA) <- c("gm.prcomp", class(PaCA))
  PaCA$phy <- phy
  names(PaCA)[[which(names(PaCA) == "xn")]] <- "anc.x"

  Y <- as.matrix(Y)
  dims <- dim(Y)
  n <- dims[1]
  p <- dims[2]
  x <- matrix(1, n)
  rownames(x) <- rownames(Y)
  U <- qr.Q(qr(Pcov %*% x))
  
  opt <- lambda.opt(Y, phy)
  
  opt.method <- if(is.numeric(lambda)) "override" else match.arg(lambda)
  opt.dim <- p
  
  if(opt.method == "override") {
    opt <- lambda
    if(opt < 0 || opt > 1)
      stop("A value of lambda was chosen that is not between 0 and 1.\n",
           call. = FALSE)
  }
  
  if(opt.method == "burn"){
    if(print.progress) {
      cat("Performing burn-in iterations:\n")
      pb <- txtProgressBar(min = 0, max = 10, initial = 0, style=3) 
    }
      
    lambdas <- seq(0.1, 1, 0.1)
    b.iter <- max(round((iter + 1) / 10), 100)
    ind.b <- perm.index(n, b.iter, seed = seed)
    res <- lapply(1:length(lambdas), function(j){
      lamb <- lambdas[[j]]
      cov.j <- updateCov(Cov, lamb)
      pcov.j <- Cov.proj(cov.j, rownames(Y))
      U.j <- qr.Q(qr(pcov.j %*% x))
      detR <- sapply(ind.b, function(jj){
        PY <- pcov.j %*% Y[jj, ]
        PR <- as.matrix(PY - fastFit(U.j, PY, n, p))
        Sig <- crossprod(PR) / n
        determinant(Sig, logarithm = TRUE)$modulus[1]
      })
      if(print.progress) setTxtProgressBar(pb,j)
      
      list(detR = detR, 
        detC = determinant(cov.j, logarithm = TRUE)$modulus[1])
      
    })
    
    if(print.progress) close(pb)
    
    const <- const <- n * p * log(2 * pi)
    
    Zs <- lapply(1:length(res), function(i){
      Res <- res[[i]]
      detC <- Res$detC
      detR <- Res$detR
      ll <- 0.5 * (const + p * detC +
                     n * detR + n * p) 
      effect.size(ll)
    })
    
    spln <- spline(lambdas, Zs)
    
    opt <- spln$x[which.max(spln$y)]
    
  }
  
  if(opt.method == "mean") {
    
    lambdas <- sapply(1:p, function(j) lambda.opt(Y[,j], phy))
    opt <- mean(lambdas)
  }
  
  if(opt.method == "PACA") {
    RV <- PaCA$RV
    prRV <- cumsum(RV)/sum(RV)
    pp <- max(c(1, length(which(prRV < 0.995))))
    pp <- min(pp, p)
    opt <- lambda.opt(Y[, 1:pp], phy)
    opt.dim <- pp
  }
  
  Cov.o <- updateCov(Cov, opt)
  Pcov.o <- Cov.proj(Cov.o, rownames(Y))
  U.o <- qr.Q(qr(Pcov.o %*% x))
  
  ind <- perm.index(n, iter, seed)
  perms <- length(ind)
  
  if(print.progress) {
    cat("Calculating log-likelihoods for", perms, "permutations:\n")
    pb <- txtProgressBar(min = 0, max = perms, initial = 0, style=3) 
  }
    
  detSig <- unlist(lapply(1:perms, function(j) {
    PY <- Pcov.o %*% Y[ind[[j]], ]
    PR <- as.matrix(PY - fastFit(U.o, PY, n, p))
    Sig <- crossprod(PR) / n
    if(print.progress) setTxtProgressBar(pb,j)
    determinant(Sig, logarithm = TRUE)$modulus[1]
  }))
  
  if(print.progress) close(pb)
  
  detCov <- determinant(Cov.o, logarithm = TRUE)$modulus[1]
  const <- const <- n * p * log(2 * pi)
  
  logLs <- -0.5 * (const + p * detCov +
                     n * detSig + n * p) 
  
  invC <- fast.solve(Cov)
  a.adj <- x %*% crossprod(x, invC) / sum(invC)
  
  Kmult <- function(Y){
    Y.c <- Y - a.adj %*% Y
    MSEobs.d <- sum(Y.c^2)  
    Y.a <- Pcov %*% Y.c
    MSE.d <- sum(Y.a^2)  
    K.denom <- (sum(diag(Cov)) - N/sum(invC)) / (N - 1)
    (MSEobs.d/MSE.d) / K.denom
  }
  
  K <- Kmult(Y)
  
  if(verbose && p > 1) {

    K.by.p <- sapply(1:p, function(j) {
      Kmult(as.matrix(Y[, 1:j]))
    })
    
    lambda.by.p <- sapply(1:p, function(j) {
      lambda.opt(Y[,1:j], phy)$lambda
    })
    
    logL.by.p <-sapply(1:p, function(j) {
      Cov.j <- updateCov(Cov, lambda.by.p[j])
      Pcov.j <- Cov.proj(Cov.j)
      detCov.j <- determinant(Cov.j, logarithm = TRUE)$modulus[1]
      U.j <- Pcov.j %*% x
      PY <- Pcov.j %*% Y[, 1:j]
      PR <- as.matrix(PY - fastFit(U.j, PY, n, j))
      Sig.j <- crossprod(PR) / n
      detSig.j <- determinant(Sig.j, logarithm = TRUE)$modulus[1]
      -0.5 * (const + j * detCov.j +
                n * detSig.j + n * j) 
    })
    
  } else {
    K.by.p <- lambda.by.p <- logL.by.p <- NULL
  }
  
  out <- list(Z = effect.size(logLs), pvalue = pval(logLs),
              rand.detR = detSig, rand.logL = logLs, 
              lambda = opt,
              K = K, PACA = PaCA, K.by.p = K.by.p,
              lambda.by.p = lambda.by.p, logL.by.p = logL.by.p,
              permutations = perms, 
              opt.method = opt.method, opt.dim = opt.dim,
              call = match.call())
  
  class(out) <- "physignal.z"
  out
}


