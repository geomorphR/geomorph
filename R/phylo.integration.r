#' Quantify phylogenetic morphological integration between two or more sets of variables under Brownian motion
#' 
#' Function quantifies the degree of phylogenetic morphological covariation between two or more sets of
#' Procrustes shape variables using partial least squares. 
#' 
#' The function quantifies the degree of phylogenetic morphological integration between two or more sets of Procrustes shape variables. 
#' The approach is based on a Brownian motion model of evolution. It is 
#'   assumed that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].   
#' 
#' The function estimates the degree of  morphological covariation between two or sets of variables 
#' while accounting for phylogeny using partial least squares (Adams and Felice 2014), and under a Brownian
#' motion model of evolution. If more than two partitions are defined, the average pairwise PLS correlation is 
#' utilized as the test statistic. The observed value is statistically assessed using permutation, where data for 
#' one partition are permuted relative to the other partitions. In addition, a multivariate effect size 
#' describing the strength of the effect is estimated from the empirically-generated sampling distribution 
#' (see details in Adams and Collyer 2019). Note that this permutation is performed on phylogenetically-
#' transformed data, so that the probability of phylogenetic association of A vs. B is similar to that of B vs. A: 
#' i.e., prob(A,B|phy)~prob(B,A|phy); thus, shuffling the correct exchangeable units under the null 
#' hypothesis of no integration (Adams and Collyer 2018). 
#' 
#'  Input for the analysis can take one of two forms. First, one can input a single dataset (as a matrix or 3D array, along with 
#'  a vector describing which variables correspond to which partitions (for the case of a 3D array, which landmarks belong to which 
#'  partitions is specified). Alternatively, when evaluating the integration between two structures or partitions, two datasets may be provided.
#'
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{phylo.integration}}.  
#'  The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which produces an ordination plot.  
#'  An additional argument allows one to include a vector to label points.  Starting with version 3.1.0, warpgrids are no longer available with \code{\link{plot.pls}}
#'  but after making a plot, the function returns values that can be used with \code{\link{picknplot.shape}} or a combination of 
#' \code{\link{shape.predictor}} and \code{\link{plotRefToTarget}} to visualize shape changes in the plot (via warpgrids).
#'  
#' \subsection{Similarity to \code{\link{two.b.pls}} and \code{\link{compare.pls}} }{ 
#' Note that \code{phylo.integration} performed on two matrices or arrays returns the same results as a phylogenetic variation of
#'  \code{\link{two.b.pls}}.  It might be of interest with 3+ modules to perform separate phylogenetic integration tests
#' between all pairwise comparisons of modules.  This can be done, test by test, and the levels of integration can be compared with
#' \code{\link{compare.pls}}.  Such results are different than using the average amount of integration when more than two modules 
#' are input, as found with \code{phylo.integration}.
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
#' @param A A 2D array (n x [p1 x k1]) or 3D array (p1 x k1 x n) containing Procrustes shape variables for the first block
#' @param A2 An optional 2D array (n x [p2 x k2]) or 3D array (p2 x k2 x n) containing Procrustes shape variables for the second block 
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param partition.gp A list of which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C). This is required when only 1 dataset provided.
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
#'   \item{Effect.Size}{The multivariate effect size associated with sigma.d.ratio.}
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
#' @references  Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological 
#' integration and trait covariation in morphometric data using evolutionary covariance 
#' matrices. PLOS ONE. 9(4):e94335.
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @references Adams, D.C. and M.L. Collyer. 2016.  On the comparison of the strength of morphological integration across morphometric 
#' datasets. Evolution. 70:2623-2631.
#' @references Adams, D.C. and M.L. Collyer. 2018. Multivariate comparative methods: evaluations, comparisons, and
#' recommendations. Systematic Biology. 67:14-31.
#' @references Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
#' alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
#' Evolution. 73:2352-2367.
#' @seealso \code{\link{integration.test}}, \code{\link{modularity.test}}, and 
#' \code{\link{two.b.pls}}
#' @examples
#' 
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
#' IT<- phylo.integration(Y.gpa$coords,partition.gp=land.gps,
#'   phy=plethspecies$phy,iter=999)
#' summary(IT) # Test summary
#' plot(IT) # PLS plot
#' 
#'  ### Visualize shape variation using picknplot.shape Because picknplot  
#'  ### requires user decisions, the following example
#'  ### is not run (but can be with removal of #).
#'  ### For detailed options, see the picknplot help file
#'  # picknplot.shape(plot(IT))
#'
phylo.integration <-function(A, A2 = NULL, phy, 
                             partition.gp = NULL,iter = 999, 
                             seed = NULL, print.progress = TRUE){ 
  if(any(is.na(A)))
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').",
         call. = FALSE)  
  
  if (!inherits(phy, "phylo"))
    stop("phy must be of class 'phylo.'") 
  
  x <- try(two.d.array(A), silent = TRUE)
  if(inherits(x, "try-error")) x <- try(as.matrix(A), silent = TRUE)
  if(inherits(x, "try-error"))
    stop("\nA is not a suitable data array for analysis. ", call. = FALSE)
  
  namesX <- rownames(x)
  
  if (is.null(namesX))
    stop("\nNo specimen names in data matrix. Please assign specimen names.", call. = FALSE)  
  
  num.taxa <- length(phy$tip.label)
  num.obs <- length(namesX)
  
  if(num.obs < num.taxa)
    stop("\nTree contains some taxa not present in present in the data matrix", call. = FALSE)  
  
  if(num.obs > num.taxa)
    stop("\nTree is missing some taxa present in the data matrix", call. = FALSE) 
  
  if(length(unique(c(namesX, phy$tip.label))) > num.taxa)
    stop("\n The data names and taxa names do not match exactly.  Check for discrepancies.",
         call. = FALSE)
  
  if(!is.null(A2)) {
    if(any(is.na(A2)))
      stop("\nData matrix 2 contains missing values. Estimate these first (see 'estimate.missing').",
           call. = FALSE) 
    
    y <- try(two.d.array(A2), silent = TRUE)
    if(inherits(y, "try-error")) y <- try(as.matrix(A), silent = TRUE)
    if(inherits(y, "try-error"))
      stop("\nA2 is not a suitable data array for analysis. ", call. = FALSE)
    
    namesY <- rownames(y)
    
    if(is.null(namesY))
      stop("\nNo specimen names in data matrix 2. Please assign specimen names",
           call. = FALSE)  
    
    if(length(unique(c(namesX, namesY))) > num.taxa)
      stop("\n The data names in the two matrices do not match.  Check for discrepancies.",
           call. = FALSE)
    
    if(length(unique(c(namesY, phy$tip.label))) > num.taxa)
      stop("\n The data names in the second matrix and taxa names do not match exactly.  
           Check for discrepancies.",
           call. = FALSE)
    
    A1.new <- A
    A2.new <- A2
    
  }
  
  #PhyloPrep  
  phy.parts<-phylo.mat(x, phy)
  invC <- phy.parts$invC
  D.mat<-phy.parts$D.mat
  
  #Analysis  
  one <- matrix(1, nrow(x))
  I <- diag(nrow(x)) 
  Ptrans <- D.mat %*% (I-one %*% crossprod(one, invC)/sum(invC))
  
  
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
    
    ngps <- nlevels(gps)
    
    if(ngps == 2){
      y <- x[, which(gps == levels(gps)[2])]
      x <- x[, which(gps == levels(gps)[1])]
    }
  }
  
  
  if(!is.null(A2)){
    ngps <- 2
    y <- as.matrix(y[match(namesX, namesY),])  
  }
  
  n <- NROW(x)
  
  if(!is.null(seed) && seed == "random") seed = sample(1:iter, 1)
  ind <- perm.index(n, iter, seed = seed)
  perms <- length(ind)
  
  if(print.progress){
    cat(paste("\nRandom PLS calculations:", perms, "permutations.\n"))
    pb <- txtProgressBar(min = 0, max = perms+1, initial = 0, style=3)
  }

  if(ngps == 2){
    pls.obs <- pls.phylo(x, y, Ptrans, verbose=TRUE)
    
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
    
    x <- center(Ptrans %*% x)
    y <- center(Ptrans %*% y)
    pls.rand <- sapply(1:perms, function(j) {
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      s <- ind[[j]]
      quick.pls(x[s,], y)
    })
    
    names(pls.rand) <- c("obs", paste("iter", 1:iter, sep = "."))
    p.vals <- NULL
    Zs <- NULL
    p.val <- pval(abs(pls.rand))
    Z <- effect.size(pls.rand, center=TRUE) 
    XScores <- pls.obs$XScores
    YScores <- pls.obs$YScores
  }
  
  if(ngps > 2){
    y <- center(Ptrans %*% x)
    g <- factor(as.numeric(gps))
    gps.combo <- combn(ngps, 2)
    pls.rand <- sapply(1:perms, function(j) {
      step <- j
      if(print.progress) setTxtProgressBar(pb,step)
      s <- ind[[j]]
      sapply(1:ncol(gps.combo), function(jj) {
        y1 <- as.matrix(y[s, g == gps.combo[1,jj]])
        y2 <- y[ , g == gps.combo[2,jj]]
        quick.pls(y1, y2)
      })
    })
    
    nms <- levels(gps)
    rnms <- apply(gps.combo, 2, function(x) paste(nms[x[1]], nms[x[2]], sep = "-"))
    cnms <- c("obs", paste("iter", 1:iter, sep = "."))
    dimnames(pls.rand) <- list(rnms, cnms)
    
    p.vals <- apply(abs(pls.rand), 1, pval)
    Zs <- apply(pls.rand, 1, effect.size)
    p.val <- pval(colMeans(abs(pls.rand)))
    Z <- effect.size(colMeans(pls.rand), center=TRUE)
    
    r.pls.mat <- matrix(0, length(nms), length(nms))
    dimnames(r.pls.mat) <- list(nms, nms)
    r.pls.mat <- as.dist(r.pls.mat)
    r.pls.mat[1:nrow(pls.rand)] <- pls.rand[, 1]
  }  
    
  step <- perms + 1
  if(print.progress) {
    setTxtProgressBar(pb,step)
    cat("\n")
    close(pb)
  }
  
  ####OUTPUT
  
  if(ngps == 2){
    out <- list(r.pls = pls.obs$r.pls, r.pls.mat = NULL, P.value = p.val, Z = Z,
                left.pls.vectors = pls.obs$left.vectors,
                right.pls.vectors = pls.obs$right.vectors,
                random.r = pls.rand, 
                XScores = pls.obs$XScores,
                YScores = pls.obs$YScores,
                svd = pls.obs$pls.svd,
                A1 = A1.new, A2 = A2.new,
                A1.matrix = x, A2.matrix =y,
                permutations = iter+1, call=match.call(),
                method = "PLS")
  }
  if(ngps>2){
    out <- list(r.pls = mean(r.pls.mat), r.pls.mat = r.pls.mat, 
                P.value = p.val, Z = Z,
                pairwise.P.values = p.vals, pairwise.Z = Zs,
                random.r = pls.rand, 
                permutations = iter+1, call=match.call(),
                method = "PLS")
  }
  
  class(out) <- "pls"
  out  
}

