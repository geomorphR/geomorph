#' Quantify morphological integration between modules
#'
#' Function quantifies the degree of morphological integration between modules of Procrustes-aligned 
#'   coordinates
#'
#' The function quantifies the degree of morphological integration between modular partitions of shape data as 
#'   defined by landmark coordinates. It is assumed that the landmarks have previously been aligned using 
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
#'   Input for the analysis can take one of two forms. First, one can input a single dataset (as a matrix or 3D array, along with 
#'  a vector describing which variables correspond to which partitions (for the case of a 3D array, which landmarks belong to which 
#'  partitions is specified). Alternatively, when evaluating the integration between two structures or partitions, two datasets may be provided.
#'
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{modularity.test}}.
#'  The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which has two additional
#'  arguments (with defaults): label = NULL, warpgrids = TRUE.  These arguments allow one to include a vector to label points and a logical statement to
#'  include warpgrids, respectively.  Warpgrids can only be included for 3D arrays of Procrustes residuals. The plot is a plot of PLS scores from 
#'  Block1 versus Block2 performed for the first set of PLS axes. 
#'  
#' @param A A 3D array (p x k x n) containing GPA-aligned coordinates for all specimens, or a matrix (n x variables)
#' @param A2 An optional 3D array (p x k x n) containing GPA-aligned coordinates for all specimens, or a matrix (n x variables) for a second partition
#' @param partition.gp A list of which landmarks (or variables) belong in which partition (e.g. A,A,A,B,B,B,C,C,C) (required when only 1 dataset provided)
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
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
#'    \item{A1}{Input values for the left block (for 2 modules only).}
#'    \item{A2}{Input values for the right block (for 2 modules only).}
#'    \item{A1.matrix}{Left block (matrix) found from A1 (for 2 modules only).}
#'    \item{A2.matrix}{Right block (matrix) found from A2 (for 2 modules only).}
#'    \item{permutations}{The number of random permutations used in the resampling procedure.}
#'    \item{call}{The match call.}
#' @references  Bookstein, F. L., P. Gunz, P. Mitteroecker, H. Prossinger, K. Schaefer, and H. Seidler. 
#'   2003. Cranial integration in Homo: singular warps analysis of the midsagittal plane in ontogeny and 
#'   evolution. J. Hum. Evol. 44:167-187.
#' @seealso \code{\link{two.b.pls}}, \code{\link{modularity.test}}, \code{\link{phylo.pls}}, and 
#' \code{\link{phylo.integration}}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#' #landmarks on the skull and mandible assigned to partitions
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
#' IT <- integration.test(Y.gpa$coords, partition.gp=land.gps, iter=999)
#' summary(IT) # Test summary
#' plot(IT) # PLS plot
#' IT$left.pls.vectors # extracting just the left (first block) singular vectors

integration.test<-function(A, A2=NULL,partition.gp=NULL,iter=999, seed=NULL){
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
    pls.rand <- apply.pls(x, y, iter=iter, seed=seed)
    pls.obs <- pls(x, y, verbose=TRUE)
    p.val <- pval(pls.rand)
    XScores <- pls.obs$XScores
    YScores <- pls.obs$YScores
  }
  if(ngps>2){
    pls.obs <- plsmulti(x, gps)  
    pls.rand <- apply.plsmulti(x, gps, iter=iter, seed=seed)
    p.val <- pval(pls.rand)
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
                A1 = A.new, A2 = A2.new,
                A1.matrix = x, A2.matrix =y,
                permutations = iter+1, call=match.call(),
                method = "PLS")
  }
  if(ngps>2){
    out <- list(r.pls = pls.obs$r.pls, r.pls.mat = r.pls.mat, P.value = p.val,
                permutations = iter+1, call=match.call(),
                method = "PLS")
  }
  
  class(out) <- "pls"
  out  
}