#' Extended Phylogenetic ANOVA/regression for Procrustes shape variables
#'
#' Function performs extended PGLS that enables more than one individual per species to be included
#' in a phylogenetic least squares analysis. 
#' 
#' The function performs linear models in a phylogenetic context in a manner that can 
#' accommodate multiple individuals per species and can accommodate high-dimensional datasets.
#' The approach utilizes a hierarchical linear model, an expanded phylogenetic covariance matrix, 
#' and permutation procedures to obtain empirical sampling distributions and effect sizes for 
#' model effects that can evaluate differences in intraspecific trends across species for 
#' both univariate and multivariate data, while conditioning them on the phylogeny (Adams and Collyer 2024).
#' 
#' Data input is specified by several components, including minimally: 1) a formula (e.g., y ~ X), where 
#' 'y' specifies the response variables (shape data), and 'X' contains one or more independent variables 
#' (discrete or continuous).  The response matrix 'Y' can be either in the form of a two-dimensional data 
#' matrix of dimension (n x [p x k]), or a 3D array (p x n x k).  It is assumed that the landmarks 
#' have previously been aligned using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].
#' 
#' 2) An argument for species is also required, which specifies the species to which each 
#' observation belongs. 
#' 
#' 3) Either phy or COV must be included, which specifies the expected phylogenetic 
#' non-independence among species. Here, 'phy' represents a phylogenetic tree of class 'phylo'. In this case,
#' the expected covariances among species will be estimated under a Brownian motion model of evolution. Alternatively, 
#' 'COV' is a phylogenetic covariance matrix obtained previously by the user, and may represent the expected
#' covariances among species under any model of evolutionary change (BM, OU, etc.).   
#' 
#' From the phylogenenetic covariance matrix, a phylogenetic transformation matrix is obtained and used to 
#' transform the Y variables and the model design matrices constructed from X, and parameter estimates
#' are obtained. Next, variance components for all model terms are calculated, using a combination of 
#' type II and type III sums of squares, so that within-species effects may be isolated from interspecies
#' effects (see Adams and Collyer 2024). Likewise, residuals from the model are permuted in such a manner
#' as to isolate intra- and inter-species effects, to account for fact that multiple observations are included
#' for each species. Additional statistical and philosophical details are found in Adams and Collyer 2024. 
#' 
#' @param f1 A formula for the linear model (e.g., y~x1+x2).  
#' @param species A variable that can be found in the data frame indicating the species 
#' to which each individual belongs. This variable must be in the data frame.  It is 
#' imperative that these names
#' match the species names in the phylogeny or phylogenetic covariance matrix. The data do not need 
#' to have row names but the species variable has to be provided.
#' @param phy A phylogenetic tree of class = "phylo" - see \code{\link[ape]{read.tree}} in library ape
#' @param Cov An argument for including a phylogenetic covariance matrix if a phylogeny is not 
#' provided. This may be obtained under any evolutionary model, and if included, 
#' any weights are ignored.  This matrix must match in dimensions the number of species.
#' @param delta A within-species scaling parameter for covariances, ranging from 
#' 0 to 1.  If delta = 0, a sight value (0.001) is added to assure variances of the 
#' covariance matrix are 0.1 percent larger than covariances.
#' @param gamma A sample-size scaling parameter that is adjusted to be 1 ("equal")
#' scaling or the square-root of the sample size for species observations ("sample").
#' @param data A data frame for the function environment, see 
#' \code{\link[RRPP]{rrpp.data.frame}}.  A data frame is required for this analysis.
#' @param print.progress A logical value to indicate whether a progress 
#' bar should be printed to the screen.
#' This is helpful for long-running analyses.
#' @param iter Number of iterations for significance testing
#' @param turbo A logical value that if TRUE, suppresses coefficient estimation 
#' in every random permutation.  This will affect subsequent analyses that 
#' require random coefficients (see \code{\link[RRPP]{coef.lm.rrpp}})
#' but might be useful for large data sets for which only ANOVA is needed.
#' @param seed An optional argument for setting the seed for random 
#' permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found 
#' for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  
#' One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param int.first A logical value to indicate if interactions of first 
#' main effects should precede subsequent main effects
#' @param verbose A logical value to indicate if all possible output from an analysis 
#' should be retained. Generally this should be FALSE, unless one wishes to extract, 
#' e.g., all possible terms, model matrices, QR decomposition, or random permutation 
#' schemes.
#' @param Parallel Either a logical value to indicate whether parallel processing 
#' should be used, a numeric value to indicate the number of cores to use, or a predefined
#' socket cluster.  
#' @param ... Arguments typically used in \code{\link{lm}}, such as 
#' weights or offset, passed on to
#' \code{\link[RRPP]{lm.rrpp}} for estimation of coefficients.  If both weights and 
#' a covariance matrix are included,
#' weights are ignored (since inverses of weights are the diagonal elements 
#' of weight matrix, used in lieu
#' of a covariance matrix.)
#' @keywords analysis
#' @export
#' @author Dean Adams   
#' 
#' @return An object of class \code{\link[RRPP]{lm.rrpp.ws}} is a list containing the 
#' following
#' \item{call}{The matched call.}
#' \item{LM}{Linear Model objects, including data (Y), coefficients, 
#' design matrix (X), sample size
#' (n), number of dependent variables (p), dimension of data space (p.prime),
#' QR decomposition of the design matrix, fitted values, residuals,
#' weights, offset, model terms, data (model) frame, random coefficients 
#' (through permutations),
#' random vector distances for coefficients (through permutations), 
#' whether OLS or GLS was performed, 
#' and the mean for OLS and/or GLS methods. Note that the data returned 
#' resemble a model frame rather than 
#' a data frame; i.e., it contains the values used in analysis, which 
#' might have been transformed according to 
#' the formula.  The response variables are always labeled Y.1, Y.2, ..., 
#' in this frame.}
#' \item{ANOVA}{Analysis of variance objects, including the SS type, 
#' random SS outcomes, random MS outcomes,
#' random R-squared outcomes, random F outcomes, random Cohen's f-squared 
#' outcomes, P-values based on random F
#' outcomes, effect sizes for random outcomes, sample size (n), number of 
#' variables (p), and degrees of freedom for
#' model terms (df).  These objects are used to construct ANOVA tables.}
#' \item{PermInfo}{Permutation procedure information, including the number 
#' of permutations (perms), The method
#' of residual randomization (perm.method), and each permutation's sampling 
#' frame (perm.schedule), which
#' is a list of reordered sequences of 1:n, for how residuals were 
#' randomized.}
#' @seealso \code{\link[RRPP]{lm.rrpp.ws}}
#' @references Adams, D.C and M.L Collyer. 2024. Extending phylogenetic regression models for 
#' comparing within-species patterns across the tree of life. 
#' Methods in Ecology and Evolution. 15:2234-2246.
#' @examples 
#' \dontrun{
#' data(pupfish.ws)
#' 
#' # With phylogeny
#' fit <- extended.pgls(f1 = coords~Species * Sex + Population, 
#'   data = pupfish.ws, species = "Species",
#'   phy = pupfish.ws$phy) 
#'   
#' anova(fit) 
#'  
#' # multivariate stats 
#' fit.mult <- manova.update(fit, PC.no = 40)
#' summary(fit.mult, test = "Wilks") 
#' 
#' # With phylogenetic covariance matrix
#' 
#' fit2 <- extended.pgls(f1 = coords ~ Species * Sex + Population, 
#'   data = pupfish.ws, species = "Species",
#'   Cov = pupfish.ws$Cov)
#'   
#' anova(fit2)
#' 
#' #sultivariate stats
#' fit2.mult <- manova.update(fit2, PC.no = 40)
#' summary(fit2.mult, test = "Wilks") 
#' }

extended.pgls<-function(f1, phy = NULL, Cov = NULL, species = NULL,
                 delta = 0.001, gamma = c("sample", "equal"), iter=999, 
                 seed=NULL, int.first = FALSE, 
                 turbo = FALSE, Parallel = FALSE,
                 verbose = FALSE,
                 data = NULL, print.progress = TRUE, ...){
  
  if(!inherits(data, "geomorph.data.frame")){
    if(!inherits(data, "data.frame") && 
       !inherits(data, "list"))
      stop(paste("\nThe data argument is neither an geomorph.data.frame",
                 "object, a data.frame object, nor a list.\n",
                 "Please see function details.\n", sep = " "), 
           call. = FALSE)
  }
  class(data) <- "rrpp.data.frame"
  
  if (is.null(phy) && is.null(Cov)) {
    stop("Must provide a phylogeny or a phylogenetic covariance matrix.\n", call. = FALSE)  }

  if(is.null(Cov)) { Cov <- fast.phy.vcv(phy)}

  if(inherits(f1, "formula")){
    Y <- try(eval(f1[[2]], envir = data , enclos = parent.frame()), silent = TRUE)
    if(inherits(Y, "try-error"))
      Y <- try(eval(f1[[2]], envir = parent.frame), silent = TRUE)
    if(inherits(Y, "try-error")) stop("Cannot find data in data frame or global environment.\n",
                                      call. = FALSE)
    nms <- get.names(Y)
    dims.Y <- dim(Y)
    f <- update(f1, Y ~ .)
    if(length(dims.Y) == 3) {
      GM <- TRUE
      Y <- two.d.array(Y) 
      rownames(Y) <- nms
      p <- dims.Y[[1]]
      k <- dims.Y[[2]]
      n <- dims.Y[[3]]
    } else {
      GM <- FALSE
      Y <- as.matrix(Y)
      rownames(Y) <- nms
      if(isSymmetric(Y)) colnames(Y) <- nms
    }
    data$Y <- Y
    
  } else {
    f <- f1
    GM <- FALSE
  }
  subjects <- species
  
  epgls.args <- list(
    f = f,
    subjects = subjects, 
    Cov = Cov, 
    delta = delta, gamma = gamma, 
    data = data,
    iter = iter, 
    seed = seed,
    Parallel = Parallel,
    turbo = turbo,
    print.progress = print.progress,
    verbose = verbose
  )
  
  out<- suppressWarnings(do.call(lm.rrpp.ws, epgls.args))
  out$call <- match.call()
  out

}
