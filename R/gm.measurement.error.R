#' This function is a wrapper for the function measurement.error in RRPP
#' 
#' Evaluation of measurement error for two or more multivariate measurements, 
#' for common research subjects.
#'
#' Function performs analyses concerned with the repeatability (reliability) of multivariate data 
#' (measurements) collected from the same research subjects.  Although there is no
#' requirement for repeated measurements on all research subjects, the analysis assumes
#' that multiple observations are made. 
#' 
#' This function performs analyses as described in Collyer and Adams (in press)
#'  to assess systematic and random components of 
#' measurement error (ME).  It basically performs ANOVA with RRPP,
#' but with different restricted randomization strategies.  The reliability of research subject variation 
#' can be considered by restricting randomization within replicates; the consistency of replicate measures
#' can be considered by restricting randomization within subjects.
#' Inter-subject variation remains constant across all random permutations within subjects and 
#' inter-replicate variation remains constant across all random permutations within replicates.  Type II
#' sums of squares and cross-products (SSCP) are calculated to assure conditional estimation.
#' 
#' The results include univariate-like statistics based on dispersion of values and
#' eigenanalysis performed on a signal to noise matrix product of SSCP matrices 
#' (sensu Bookstein and Mitteroecker, 2014) 
#' including the inverse of the random component of ME and the systematic
#' component of ME.  The multivariate test is a form of multivariate ANOVA (MANOVA), using
#' RRPP to generate sampling distributions of the major eigenvalue (Roy's maximum root).
#' 
#' @param coords A 3D array (p x k x n) containing Procrustes shape variables for all specimens
#' @param subjects A vector or factor of research subjects (each subject should occur twice or more).  
#' The length of the vector must equal the number of observations and will be coerced into a factor.
#' @param replicates A vector or factor for replicate measurements for research subjects.  
#' The length of the vector must equal the number of observations and will be coerced into a factor.
#' @param groups An optional vector, coercible to factor, to be included in the linear model
#' (as an interaction with replicates)..
#' This would be of interest if one were concerned with systematic ME occurring perhaps differently among 
#' certain strata within the data.  For example, systematic ME because of an observer bias might
#' only be observed with females or males.  
#' @param data An data frame of class \code{\link{geomorph.data.frame}}.
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random 
#' permutations of the resampling procedure.
#' If left NULL (the default), the exact same P-values will be found 
#' for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  
#' One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param multivariate Logical value for whether to include multivariate analyses.  Intraclass correlation 
#' matrices and relative eigenanalysis are based on products of sums of squares and cross-products (SSCP)
#' matrices, some of which must be inverted and potentially require
#'  significant computation time.  If FALSE, only statistics based on dispersion of values are calculated.
#' @param use.PCs A logical argument for whether to use the principal components of the data.  
#' This might be helpful for relative eigenanalysis, and if p > n, 
#' in which case inverting singular covariance matrices would not be possible.
#' @param tol A value indicating the magnitude below which 
#' components should be omitted., if use.PCs is TRUE. (Components are omitted if their 
#' standard deviations are less than or equal to tol times the 
#' standard deviation of the first component.)  See \code{\link[RRPP]{ordinate}} for more details.
#' @param Parallel The same argument as in \code{\link[RRPP]{lm.rrpp}} to govern parallel processing (
#' either a logical vale -- TRUE or FALSE -- or the number of threaded cores to use).  See \code{\link[RRPP]{lm.rrpp}} 
#' for additional details.
#' @param turbo Logical value for whether to suppress coefficient estimation in RRPP iteration,
#' thus turbo-charging RRPP.
#' @param print.progress A logical value to indicate whether a progress 
#' bar should be printed to the screen.
#' @param verbose A logical value to indicate if all the output from an
#' \code{\link[RRPP]{lm.rrpp}} analysis should be retained.  If FALSE, only the needed
#' output for summaries and plotting is retained.
#' @param ... Arguments passed on to \code{\link[RRPP]{lm.rrpp.ws}}.
#' @export
#' @keywords analysis
#' @author Michael Collyer and Dean Adams
#' @return Objects of class "measurement.error" return the same objects
#' as a \code{\link[RRPP]{lm.rrpp}} fit, plus a list of the following:
#'  \item{AOV}{Analysis of variance to test for systematic error, based on dispersion of values.}
#'  \item{mAOV}{Multivariate AOV based on product of the inverse of the random component (SSCP) of ME
#'  times the systematic component of ME.}
#'  \item{SSCP}{The sums of squares and cross-products matrices for model effects.}
#'  \item{SSCP.ME.product}{The products of the inverse of the random ME SSCP and the SSCP matrices
#'  for systematic ME,.  These are the same matrix products used for eigenanalysis.  
#'  This is the observed matrix.}
#'  \item{SSCP.ME.product.std}{A list of the symmetric forms of standardized SSCP.ME.products 
#'  that yield orthogonal eigenvectors.}

#' @references Collyer, M.L. and D.C. Adams.  2024. Interrogating Random and Systematic Measurement Error 
#' in Morphometric Data. Evolutionary Biology.
#' @references Bookstein, F.L., & Mitteroecker, P. (2014). Comparing covariance matrices by relative eigenanalysis, 
#' with applications to organismal biology. Evolutionary biology, 41(2), 336-350.

#' @seealso \code{\link[RRPP]{lm.rrpp.ws}}, \code{\link[RRPP]{manova.update}}

#' @examples
#' \dontrun{
#' # Measurement error analysis on simulated data of fish shapes
#' 
#' data(fishy)
#' fishy$coordsarray <- arrayspecs(fishy$coords, p = 11, k = 2)  #make 3D array
#' 
#' # Example two digitization replicates of the same research subjects
#' rep1 <- matrix(fishy$coords[1,], 11, 2, byrow = TRUE)
#' rep2 <- matrix(fishy$coords[61,], 11, 2, byrow = TRUE)
#' plot(rep1, pch = 16, col = gray(0.5, alpha = 0.5), cex = 2, asp = 1)
#' points(rep2, pch = 16, col = gray(0.2, alpha = 0.5), cex = 2, asp = 1)
#' 
#' # Analysis unconcerned with groups 
#' 
#' ME1 <- gm.measurement.error(
#'   coords = "coordsarray",
#'   subjects = "subj",
#'   replicates = "reps",
#'   data = fishy)
#' 
#' anova(ME1)
#' ICCstats(ME1, subjects = "Subjects", with_in = "Systematic ME")
#' plot(ME1)
#' 
#' # Analysis concerned with groups 
#' 
#' ME2 <- gm.measurement.error(
#'   coords = "coordsarray",
#'   subjects = "subj",
#'   replicates = "reps",
#'   groups = "groups",
#'   data = fishy)
#'   
#' anova(ME2)
#' ICCstats(ME2, subjects = "Subjects", 
#'   with_in = "Systematic ME", groups = "groups")
#' P <- plot(ME2)
#' focusMEonSubjects(P, subjects = 18:20, shadow = TRUE)
#' 
#' #heat map of inter-subject variability
#' int.var <- interSubVar(ME2, type = "var")
#' plot(int.var)
#' }
#' 
gm.measurement.error <- function(coords, 
                                 subjects, 
                                 replicates, 
                                 groups = NULL,
                                 data,
                                 iter = 999, 
                                 seed = NULL,
                                 multivariate = FALSE,
                                 use.PCs = TRUE, 
                                 tol = 0.001, 
                                 Parallel = FALSE,
                                 turbo = TRUE,
                                 print.progress = FALSE,
                                 verbose = FALSE, ...) {
  
  if(!inherits(data, "geomorph.data.frame")){
    if(!inherits(data, "data.frame") && 
       !inherits(data, "list"))
      stop(paste("\nThe data argument is neither an geomorph.data.frame",
                 "object, a data.frame object, nor a list.\n",
                 "Please see function details.\n", sep = " "), 
           call. = FALSE)
  }
  class(data) <- "rrpp.data.frame"
  
  Y <- as.character(coords)
  Yslot <- which(names(data) %in% Y)
  if(length(Yslot) == 0)
    stop(paste("\nThe Y argument must be a character",
               "value, like 'myData',",
               "which can be found in the RRPP data frame.\n", 
               sep = " "), call. = FALSE)
  
  Y <- data[Yslot][[1]]
  dims <- dim(Y)    
  if(length(dims) == 3) {
    data$Y.2D <- two.d.array(Y)      
  } else{
    stop("Data not a 3D array.\n")
  } 
  me.args <- c(as.list(environment()), list(...))
  me.args <- list(
    Y = "Y.2D",
    subjects = subjects, 
    replicates = replicates, 
    groups = groups,
    data = data,
    iter = iter, 
    seed = seed,
    multivariate = multivariate,
    use.PCs = use.PCs, 
    tol = tol, 
    Parallel = Parallel,
    turbo = turbo,
    print.progress = print.progress,
    verbose = verbose
  )

  out<- suppressWarnings(do.call(measurement.error, me.args))
  out
}