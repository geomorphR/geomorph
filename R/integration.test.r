#' Quantify morphological integration between modules
#'
#' Function quantifies the degree of morphological integration between modules of Procrustes shape variables
#'
#' The function quantifies the degree of morphological integration between modular partitions of shape data as 
#'   defined by Procrustes shape variables. It is assumed that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The function may be used to assess
#'   the degree of morphological integration between two or more sets of variables. 
#'   
#'   The function estimates the degree of morphological integration using a two-block partial least squares 
#'   analysis (PLS). When used with landmark data, this analysis is referred to as singular warps analysis 
#'   (Bookstein et al. 2003). If more than two partitions are defined, the average pairwise PLS correlation is utilized as
#'   the test statistic. The observed test value is then compared to a distribution of values obtained by randomly permuting 
#'   the individuals (rows) in one partition relative to those in the other. A significant result is found when the 
#'   observed PLS correlation is large relative to this distribution, and implies that the structures are integrated with one
#'   another (see Bookstein et al. 2003).  If only two partitions are specified, a plot of PLS scores along the first 
#'   set of PLS axes is optionally displayed, and thin-plate spline deformation grids along these axes are also shown if data were 
#'   input as a 3D array.  
#'   
#'  Input for the analysis can take one of two forms. First, one can input a single dataset (as a matrix or 3D array, along with 
#'  a vector describing which variables correspond to which partitions (for the case of a 3D array, which landmarks belong to which 
#'  partitions is specified). Alternatively, when evaluating the integration between two structures or partitions, two datasets may be provided.
#'
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{integration.test}}.
#'  The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which has two additional
#'  arguments (with defaults): label = NULL, warpgrids = TRUE.  These arguments allow one to include a vector to label points and a logical statement to
#'  include warpgrids, respectively.  Warpgrids can only be included for 3D arrays of Procrustes residuals. The plot is a plot of PLS scores from 
#'  Block1 versus Block2 performed for the first set of PLS axes. 
#'  
#' \subsection{Similarity to \code{\link{two.b.pls}} and \code{\link{compare.pls}} }{ 
#' Note that \code{integration.test} performed on two matrices or arrays returns the same results as \code{\link{two.b.pls}}.  
#' However,  \code{\link{two.b.pls}} is limited to only two modules.  It might be of interest with 3+ modules to perform integration tests
#' between all pairwise comparisons of modules.  This can be done, test by test, and the levels of integration can be compared with
#' \code{\link{compare.pls}}.  Such results are different than using the average amount of integration, as performed by \code{integration.test}
#' when more than two modules are input.
#' }
#' 
#'  \subsection{Using phylogenies and PGLS}{ 
#' If one wishes to incorporate a phylogeny, \code{\link{phylo.integration}} is the function to use.  This function is exactly the same as \code{integration.test}
#' but allows PGLS estimation of PLS vectors.  Because \code{\link{integration.test}} can be used on two blocks, \code{\link{phylo.integration}} likewise allows one to
#' perform a phylogenetic two-block PLS analysis.
#' }
#'  
#'  \subsection{Notes for geomorph 3.0.4 and subsequent versions}{ 
#'  Compared to previous versions of geomorph, users might notice differences in effect sizes.  Previous versions used z-scores calculated with 
#'  expected values of statistics from null hypotheses (sensu Collyer et al. 2015); however Adams and Collyer (2016) showed that expected values 
#'  for some statistics can vary with sample size and variable number, and recommended finding the expected value, empirically, as the mean from the set 
#'  of random outcomes.  Geomorph 3.0.4 and subsequent versions now center z-scores on their empirically estimated expected values and where appropriate, 
#'  log-transform values to assure statistics are normally distributed.  This can result in negative effect sizes, when statistics are smaller than 
#'  expected compared to the average random outcome.  For ANOVA-based functions, the option to choose among different statistics to measure effect size 
#'  is now a function argument.
#' }
#' 
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param A2 An optional 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables) for a second partition
#' @param partition.gp A list of which landmarks (or variables) belong in which partition (e.g. A,A,A,B,B,B,C,C,C) (required when only 1 dataset provided)
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Objects of class "pls" from integration.test return a list of the following:
#'  \item{r.pls}{The estimate of morphological integration: PLS.corr. The mean of pairwise
#'  PLS correlations between partitions is used when there are more than two partitions.}
#'    \item{r.pls.mat}{The pairwise r.pls, if the number of partitions is greater than 2.}
#'    \item{P.value}{The empirically calculated P-value from the resampling procedure.}
#'    \item{left.pls.vectors}{The singular vectors of the left (x) block (for 2 modules only).}
#'    \item{right.pls.vectors}{The singular vectors of the right (y) block (for 2 modules only).}
#'    \item{random.r}{The correlation coefficients found in each random permutation of the 
#'   resampling procedure.}
#'    \item{XScores}{Values of left (x) block projected onto singular vectors 
#'   (for 2 modules only).}
#'    \item{YScores}{Values of right (y) block projected onto singular vectors
#'   (for 2 modules only).}
#'    \item{svd}{The singular value decomposition of the cross-covariances (for 2 modules only).}
#'    \item{A1}{Input values for the left block (for 2 modules only).}
#'    \item{A2}{Input values for the right block (for 2 modules only).}
#'    \item{A1.matrix}{Left block (matrix) found from A1 (for 2 modules only).}
#'    \item{A2.matrix}{Right block (matrix) found from A2 (for 2 modules only).}
#'    \item{permutations}{The number of random permutations used in the resampling procedure.}
#'    \item{call}{The match call.}
#' @references  Bookstein, F. L., P. Gunz, P. Mitteroecker, H. Prossinger, K. Schaefer, and H. Seidler. 
#'   2003. Cranial integration in Homo: singular warps analysis of the midsagittal plane in ontogeny and 
#'   evolution. J. Hum. Evol. 44:167-187.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
#' @seealso \code{\link{two.b.pls}}, \code{\link{modularity.test}}, 
#' \code{\link{phylo.integration}}, and \code{\link{compare.pls}}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' #landmarks on the skull and mandible assigned to partitions
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
#' IT <- integration.test(Y.gpa$coords, partition.gp=land.gps, iter=999)
#' summary(IT) # Test summary
#' plot(IT) # PLS plot
#' IT$left.pls.vectors # extracting just the left (first block) singular vectors

integration.test<-function(A, A2=NULL,partition.gp=NULL,iter=999, seed=NULL, print.progress=TRUE){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(!is.null(partition.gp)){
    partition.gp<-as.factor(partition.gp)
    if (length(dim(A))==3){ x<-two.d.array(A)
    p<-dim(A)[1]; k<-dim(A)[2];n<-dim(A)[3]
    if(length(partition.gp)!=p){stop("Not all landmarks are assigned to a partition.")}
    gps<-as.factor(rep(partition.gp,k,each = k, length=p*k))  
    A.new<-A[which(partition.gp==levels(partition.gp)[1]),,]
    A2.new<-A[which(partition.gp!=levels(partition.gp)[1]),,]
    }
    if (length(dim(A))==2){ x<-A; A.new<-A
    if(length(partition.gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
    gps<-as.factor(partition.gp) ;n<-dim(x)[2] 
    A.new<-x[,which(gps==levels(gps)[1])]; A2.new<-x[,which(gps==levels(gps)[2])]
    }
    ngps<-nlevels(gps)
    if(ngps==2){
      y<-x[,which(gps==levels(gps)[2])]
      x<-x[,which(gps==levels(gps)[1])]
    }
  }
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  if(!is.null(A2)){
    A.new<-A; A2.new<-A2
    if(any(is.na(A2))==T){
      stop("Data matrix 2 contains missing values. Estimate these first (see 'estimate.missing').")  }
    if (length(dim(A))==2){ x<-A }
    if (length(dim(A))==3){ x<-two.d.array(A)}
    if (length(dim(A2))==2){ y<-A2}
    if (length(dim(A2))==3){ y<-two.d.array(A2)}
    ngps=2; n<-dim(x)[2]
  }
  if(ngps==2){
    pls.obs <- pls(x, y, verbose=TRUE)
    if(NCOL(x) > NROW(x)){
      pcax <- prcomp(x)
      d <- which(zapsmall(pcax$sdev) > 0)
      x <- pcax$x[,d]
    }
    if(NCOL(y) > NROW(y)){
      pcay <- prcomp(y)
      d <- which(zapsmall(pcay$sdev) > 0)
      y <- pcay$x[,d]
    }
    if(print.progress) pls.rand <- apply.pls(center(x), center(y), iter=iter, seed=seed) else
      pls.rand <- .apply.pls(center(x), center(y), iter=iter, seed=seed)
    p.val <- pval(abs(pls.rand))
    XScores <- pls.obs$XScores
    YScores <- pls.obs$YScores
  }
  if(ngps>2){
    pls.obs <- plsmulti(x, gps)  
    if(print.progress) pls.rand <- apply.plsmulti(center(x), gps, iter=iter, seed=seed) else
      pls.rand <- .apply.plsmulti(center(x), gps, iter=iter, seed=seed)
    p.val <- pval(abs(pls.rand))
  }  
  ####OUTPUT
  if(ngps > 2) r.pls.mat <- pls.obs$r.pls.mat else r.pls.mat <- NULL
  if(ngps==2){
    out <- list(r.pls = pls.obs$r.pls, r.pls.mat = r.pls.mat, P.value = p.val,
                left.pls.vectors = pls.obs$left.vectors,
                right.pls.vectors = pls.obs$right.vectors,
                random.r = pls.rand, 
                XScores = pls.obs$XScores,
                YScores = pls.obs$YScores,
                svd = pls.obs$pls.svd,
                A1 = A.new, A2 = A2.new,
                A1.matrix = x, A2.matrix =y,
                permutations = iter+1, call=match.call(),
                method = "PLS")
  }
  if(ngps>2){
    out <- list(r.pls = pls.obs$r.pls, r.pls.mat = r.pls.mat, P.value = p.val,
                random.r = pls.rand, 
                permutations = iter+1, call=match.call(),
                method = "PLS")
  }
  
  class(out) <- "pls"
  out  
}