#' Quantify phylogenetic morphological integration between two or more sets of variables under Brownian motion
#' 
#' Function quantifies the degree of phylogenetic morphological covariation between two or more sets of
#' Procrustes-aligned coordinates using partial least squares. 
#' 
#' The function quantifies the degree of phylogenetic morphological integration between two or more sets of shape data as 
#'   defined by landmark coordinates. The approach is based on a Brownian motion model of evolution. It is 
#'   assumed that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].   
#' 
#' The function estimates the degree of  morphological covariation between two or sets of variables 
#' while accounting for phylogeny using partial least squares (Adams and Felice 2014), and under a Brownian
#' motion model of evolution. If more than two partitions are defined, the average pairwise PLS correlation is 
#' utilized as the test statistic. The observed value is statistically assessed using permutation, where data for 
#' one partition are permuted relative to the other partitions. Note that this permutation is performed on phylogenetically-
#' transformed data, so that the probability of phylogenetic association of A vs. B is similar to that of B vs. A: 
#' i.e., prob(A,B|phy)~prob(B,A|phy).  
#' 
#'   Input for the analysis can take one of two forms. First, one can input a single dataset (as a matrix or 3D array, along with 
#'  a vector describing which variables correspond to which partitions (for the case of a 3D array, which landmarks belong to which 
#'  partitions is specified). Alternatively, when evaluating the integration between two structures or partitions, two datasets may be provided.
#'
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{phylo.integration}}.
#'  The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which has two additional
#'  arguments (with defaults): label = NULL, warpgrids = TRUE.  These arguments allow one to include a vector to label points and a logical statement to
#'  include warpgrids, respectively.  Warpgrids can only be included for 3D arrays of Procrustes residuals. The plot is a plot of PLS scores from 
#'  Block1 versus Block2 performed for the first set of PLS axes. 
#'  
#' @param A A 2D array (n x [p1 x k1]) or 3D array (p1 x k1 x n) containing landmark coordinates for the first block
#' @param A2 An optional 2D array (n x [p2 x k2]) or 3D array (p2 x k2 x n) containing landmark coordinates for the second block 
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
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
#' @references  Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological 
#' integration and trait covariation in morphometric data using evolutionary covariance 
#' matrices. PLOS ONE. 9(4):e94335.
#' @seealso \code{\link{integration.test}}, \code{\link{modularity.test}}, \code{\link{phylo.pls}}, and 
#' \code{\link{two.b.pls}}
#' @examples
#' 
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
#' IT<- phylo.integration(Y.gpa$coords,partition.gp=land.gps,phy=plethspecies$phy,iter=999)
#' summary(IT) # Test summary
#' plot(IT) # PLS plot
#' 
phylo.integration <-function(A, A2=NULL, phy, partition.gp=NULL,iter=999, seed=NULL){ 
  if(any(is.na(A))==T){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  } 
  if (!is.phylo(phy))
    stop("phy must be of class 'phylo.'") 
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
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
    Nspec<-num.taxa.X<-num.taxa.Y<-nrow(x)
    namesX<-namesY<-rownames(x)
  }
  if(!is.null(A2)){
    A.new<-A; A2.new<-A2
    if(any(is.na(A2))==T){
      stop("Data matrix 2 contains missing values. Estimate these first (see 'estimate.missing').")  }
    if (length(dim(A))==2){ x<-A }
    if (length(dim(A))==3){ x<-two.d.array(A)}
    if (length(dim(A2))==2){ y<-A2}
    if (length(dim(A2))==3){ y<-two.d.array(A2)}
    ngps=2; n<-dim(x)[2]
    Nspec<-num.taxa.X<-nrow(x)
    namesX<-rownames(x)
    num.taxa.Y<-nrow(y)
    namesY<-rownames(y)
    y<-y[rownames(x),] 
  }
  if (is.null(namesX)){
    stop("No specimen names in data matrix. Please assign specimen names.")  } 
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && length(phy$tip.label) < num.taxa.X)
    stop("Tree is missing some taxa present in the data matrix") 
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && num.taxa.X < length(phy$tip.label)) 
    stop("Tree contains some taxa not present in present in the data matrix")  
  if (is.null(namesY)){
    stop("No specimen names in data matrix 2. Please assign specimen names")  } 
  if (is.null(namesX) == FALSE && is.null(namesY) == FALSE) {
    mtch.A <- namesX[is.na(match(namesX, namesY))]
    if (length(mtch.A) > 0) {
      stop("Specimen names in data sets are not the same.")   } 
  }
  mtch.B <- namesX[is.na(match(namesX, phy$tip.label))]
  if (length(mtch.B) > 0) {
    stop("Taxa labels on tree and taxa matrix are not the same.")} 
#PhyloPrep  
  phy.parts<-phylo.mat(x,phy)
  invC<-phy.parts$invC; D.mat<-phy.parts$D.mat
#Analysis  
  if(ngps==2){
    pls.obs <- pls.phylo(x, y, invC,D.mat,verbose=TRUE)
    pls.rand <- apply.pls.phylo(x, y,invC,D.mat, iter=iter, seed=seed)
    p.val <- pval(pls.rand)
    XScores <- pls.obs$XScores
    YScores <- pls.obs$YScores
  }
  if(ngps>2){
    pls.obs <- plsmulti.phylo(x, gps, invC,D.mat)  
    pls.rand <- apply.plsmulti.phylo(x, gps, invC,D.mat, iter=iter, seed=seed)
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

