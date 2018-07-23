#' Principal components analysis of shape data
#'
#' Function performs raw or weighted PCA on Procrustes shape coordinates  
#'
#' The function performs a principal components analysis of shape variation, with the possibility
#' of weighing the analysis using a phylogenetic tree, or a variance-covariance matrix, 
#' to incorporate the expected covariance among dependent observations 
#' (e.g. due to spatial autocorrelation). At present combined analysis allowing the use of both 
#' a phylogeny and another covariance matrix is not implemented.
#' Several complementary types of analyses are available when using a phylogeny, namely:
#'      1) phyloPCA: A phylogenetic PCA, where a phylogenetic transformation is applied to take expected
#'      covariance due to shared phylogenetic history into account (sensu Revell 2009).
#'      2) phylomorphospace: A non-weighted PCA of the raw data, followed by a projection to PCA space of 
#'      the estimated ancestral shape values for the nodes of the phylogeny.
#'      3) alldataPCA: A non-weighted PCA of both the tip and estimated ancestor shapes together. This is
#'      a legacy version of the operations previously included in plotGMPhyloMorphoSpace.
#' 
#' WARNING: Note that in the case of a phylogenetic PCA, the sum of eigenvalues does not correspond to the
#' total variance of the data, or the variance of PC scores. Also, note that due to the phylogenetic centering
#' and weighting applied, this method distorts shape space. For more details on phylogenetic PCA and its
#' properties, users are strongly encouraged to consult Polly et al. 2013.
#' 
#' PLOTTING: Contrary to previous geomorph implementations, gm.prcomp does not produce plots. 
#' For plotting options of gm.prcomp class objects combine \code{\link{plot.gm.prcomp}} and 
#' \code{\link{picknplot.shape}} following the examples below. Note that trying to plot the full result
#' of gm.prcomp will give an error. Choose a specific PCA method to be plotted, by pointing to one of the 
#' components of the list returned by gm.prcomp.
#' 
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of aligned specimens
#' @param phy An optional phylogenetic tree of class phylo - see \code{\link{read.tree}} in library ape
#' @param Cov An optional covariance matrix for weighting. See also details.
#' @param ... Other arguments passed to \code{\link{scale}}
#' @return An object of class "gm.prcomp" contains a list of results for each of the PCA approaches implemented.
#' Each of these lists includes the following components:
#' \item{x}{Principal component scores for all specimens.}
#' \item{anc.x}{Principal component scores for the ancestors on the phylogeny.}
#' \item{d}{The singular values of the decomposed VCV matrix.}
#' \item{rotation}{The matrix of variable loadings, i.e. the eigenvectors of the decomposed matrix.}
#' \item{shapes}{A list with the shape coordinates of the extreme ends of all PC axes.}
#' \item{ancestors}{The matrix of estimated ancestral shapes, if a phylogeny is used.}
#' @export
#' @keywords visualization
#' @seealso \code{\link{plot.gm.prcomp}}
#' @seealso \code{\link{picknplot.shape}}
#' @author Antigoni Kaliontzopoulou, Michael Collyer & Dean Adams
#' @references Revell, L.J. (2009) Size-correction and principal components for interspecific comparative studies. 
#' Evolution, 63, 3258-3268.
#' Polly, P.D.P. et al. (2013) Phylogenetic principal components analysis and geometric morphometrics, 
#' Hystrix, the Italian Journal of Mammalogy, 24(1), 33–41.
#' @examples
#'  data(plethspecies) 
#'  Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
#'  
#'  ### PCA on raw data
#'  pleth.raw <- gm.prcomp(Y.gpa$coords)
#'  summary(pleth.raw)
#'  
#'  ### PCA with phylogeny
#'  pleth.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#'  summary(pleth.phylo)
#'  summary(pleth.phylo$phylomorphospace)
#'  
#'  #### Plotting
#'  plot(pleth.raw$rawPCA)
#'  gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
#'  par(mar=c(2, 2, 2, 2))
#'  plot(pleth.raw$rawPCA, pch=22, cex = 1.5, bg = gps) # Modify options as desired
#'  #  Add things as desired using standard R plotting
#'  text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", pos = 4, font = 2)
#'  text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
#'  legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))
#'  
#'  ### Plotting including a phylogeny
#'  plot(pleth.phylo$phyloPCA, phylo = TRUE, cex = 1.5, pch = 22, bg = gps, xlim = c(-0.05, 0.015),
#'       phylo.par = list(edge.color = "grey", edge.width = 2,
#'                        node.bg = "black", node.pch = 22, node.cex = 0.5))
#'  text(pleth.phylo$phyloPCA$x, labels = labels(pleth.phylo$phyloPCA$x)[[1]],
#'       pos = 2, font = 4) 
#'  text(pleth.phylo$phyloPCA$anc.x, labels = 1:nrow(pleth.phylo$phyloPCA$anc.x),
#'       adj = c(-0.1, -0.1), font = 2) 
#'  
#'  ### Visualize shape variation using picknplot.shape 
#'  # For detailed options, see the picknplot help file
#'  picknplot.shape(plot(pleth.phylo$phyloPCA))
#'  
#'  ### All methods 
#'  layout(matrix(1:4, nrow=2, byrow=TRUE))
#'  plot(pleth.phylo$rawPCA, pch=21, bg=1:nrow(pleth.phylo$rawPCA$x)); title(main="rawPCA")
#'  plot(pleth.phylo$phyloPCA, pch=21, bg=1:nrow(pleth.phylo$rawPCA$x), phylo = TRUE,
#'       phylo.par = list(edge.color="grey", node.cex=0)); title(main="phyloPCA")
#'  plot(pleth.phylo$phylomorphospace, pch=21, bg=1:nrow(pleth.phylo$rawPCA$x), phylo = TRUE,
#'       phylo.par = list(edge.color="grey", node.cex=0)); title(main="phylomorphospace")
#'  plot(pleth.phylo$alldataPCA, pch=21, bg=1:nrow(pleth.phylo$rawPCA$x), phylo = TRUE,
#'       phylo.par = list(edge.color="grey", node.cex=0)); title(main="alldataPCA")
#' 
#' 

gm.prcomp <- function (A, phy = NULL, Cov = NULL, ...) {
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }
  if(!is.null(phy) & !is.null(Cov)){
    stop("Method not implemented for weighting with BOTH a phylogeny and another covariance matrix.
         Only the phylogeny will be considered for weighting the PCA.")
  }
  
  p <- dim(A)[1]; k <- dim(A)[2]; n <- dim(A)[3]
  dots <- list(...)
  scale. <- dots$scale.
  if(is.null(scale.)) scale. <- FALSE
  center <- dots$center
  if(is.null(center)) center <- TRUE
  Y <- scale(two.d.array(A), center = center, scale = scale.)
  vcv.y <- crossprod(Y) / (n-1)
  svd.y <- svd(vcv.y)

  if(!is.null(phy)){
    if (!inherits(phy, "phylo"))
      stop("Tree must be of class 'phylo.'")
    N <- length(phy$tip.label); Nnode <- phy$Nnode
    if(N!=n)
      stop("Number of taxa in data matrix and tree are not equal.")
    if(is.null(rownames(Y))) {
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")
    } else {
      A <- A[,,phy$tip.label]
      Y <- Y[phy$tip.label, ]
    }
    anc.raw <- shape.ace(two.d.array(A), phy)
    ancY <- shape.ace(Y, phy)
    all.dt <- rbind(Y, ancY)
    
    C  <- vcv.phylo(phy)
    ones <- matrix(1, ncol(C))
    invC <- fast.solve(C)
    rootY <- fast.solve(t(ones)%*%invC%*%ones)%*%t(ones)%*%invC%*%Y
    
    Yc.phylo <- scale(Y, center = rootY, scale = F)
    anc.phylo <- scale(ancY, center = rootY, scale = F)
    vcv.p <- (t(Yc.phylo)%*%invC%*%Yc.phylo) / (nrow(Yc.phylo) - 1)
    svd.p <- svd(vcv.p)
    
    vcv.a <- crossprod(all.dt) / (nrow(all.dt)-1)
    svd.a <- svd(vcv.a)
  }
  
  if(!is.null(Cov) & is.null(phy)){
    C <- Cov
    C <- C[rownames(Y), rownames(Y)]
    ones <- matrix(1, ncol(C))
    invC <- fast.solve(C)
    w.meanY <- fast.solve(t(ones)%*%invC%*%ones)%*%t(ones)%*%invC%*%Y
    Ywc <- scale(Y, center = w.meanY, scale = F)
    vcv.w <- (t(Ywc)%*%invC%*%Ywc) / (N-1)
    svd.w <- svd(vcv.w)
  }
  
  ### Output
  rawPCA <- phyloPCA <- phylomorphospace <- alldataPCA <- wPCA <- NULL
  
  # Raw data PCA
  nPC <- max(which(zapsmall(svd.y$d) > 0))
  rawPCA$x <- Y%*%svd.y$v[,1:nPC]
  colnames(rawPCA$x) <- paste("PC", 1:nPC, sep = "")
  rawPCA$d <- svd.y$d[1:nPC]
  rawPCA$rotation <- svd.y$v
  rawPCA$shapes <- lapply(1:ncol(rawPCA$x),  
                         function(x){shape.predictor(A, rawPCA$x[,x], 
                                                     min = min(rawPCA$x[,x]),
                                                     max = max(rawPCA$x[,x]))})
  names(rawPCA$shapes) <- paste("shapes.PC", 1:nPC, sep = "")
  class(rawPCA) <- c("gm.prcomp", "rawPCA")
  
  if(!is.null(phy)){
    # Phylomorphospace
    phylomorphospace$x <- (all.dt%*%svd.y$v[, 1:nPC])[1:N,]
    phylomorphospace$anc.x <- (all.dt%*%svd.y$v[, 1:nPC])[(N+1):(N+Nnode),]
    colnames(phylomorphospace$x) <- colnames(phylomorphospace$anc.x) <- paste("PC", 1:nPC, sep = "")
    phylomorphospace$d <- svd.y$d[1:nPC]
    phylomorphospace$rotation <- svd.y$v
    phylomorphospace$shapes <- lapply(1:ncol(phylomorphospace$x), 
                                      function(x){shape.predictor(A, phylomorphospace$x[,x],
                                                                  min = min(phylomorphospace$x[,x]),
                                                                  max = max(phylomorphospace$x[,x]))})
    names(phylomorphospace$shapes) <- paste("shapes.PC", 1:nPC, sep = "")
    phylomorphospace$ancestors <- anc.raw
    class(phylomorphospace) <- c("gm.prcomp", "phylomorphospace")
    
    # Phylogenetic PCA
    nPC <- max(which(zapsmall(svd.p$d) > 0))
    phyloPCA$x <- (Yc.phylo%*%svd.p$v[, 1:nPC])
    phyloPCA$anc.x <- (anc.phylo%*%svd.p$v[, 1:nPC])
    phyloPCA$d <- svd.p$d[1:nPC]
    phyloPCA$rotation <- svd.p$v
    colnames(phyloPCA$x) <- colnames(phyloPCA$anc.x) <- paste("PC", 1:nPC, sep = "")
    phyloPCA$shapes <- lapply(1:ncol(phyloPCA$x),
                              function(x){shape.predictor(A, phyloPCA$x[,x],
                                                          min = min(phyloPCA$x[,x]),
                                                          max = max(phyloPCA$x[,x]))})
    names(phyloPCA$shapes) <- paste("shapes.PC", 1:nPC, sep = "")
    phyloPCA$ancestors <- anc.raw
    class(phyloPCA) <- c("gm.prcomp", "phyloPCA")
    
    # All data PCA
    nPC <- max(which(zapsmall(svd.a$d) > 0))
    alldataPCA$x <- (all.dt%*%svd.a$v[, 1:nPC])[1:N,]
    alldataPCA$anc.x <- (all.dt%*%svd.a$v[, 1:nPC])[(N+1):(N+Nnode),]
    alldataPCA$x <- scale(alldataPCA$x, center = alldataPCA$anc.x[1,], scale = F)
    alldataPCA$anc.x <- scale(alldataPCA$anc.x, center = alldataPCA$anc.x[1,], scale = F)
    alldataPCA$d <- svd.a$d[1:nPC]
    alldataPCA$rotation <- svd.a$v
    colnames(alldataPCA$x) <- colnames(alldataPCA$anc.x) <- paste("PC", 1:nPC, sep = "")
    alldataPCA$shapes <- lapply(1:ncol(alldataPCA$x),
                             function(x){shape.predictor(A, alldataPCA$x[,x],
                                                         min = min(alldataPCA$x[,x]),
                                                         max = max(alldataPCA$x[,x]))})
    names(alldataPCA$shapes) <- paste("shapes.PC", 1:nPC, sep = "")
    alldataPCA$ancestors <- anc.raw
    class(alldataPCA) <- c("gm.prcomp", "alldataPCA")
  }
  
  if(!is.null(Cov) & is.null(phy)){
    # Weighted PCA
    nPC <- max(which(zapsmall(svd.w$d) > 0))
    wPCA$x <- Ywc%*%svd.w$v[,1:nPC]
    wPCA$d <- svd.w$d[1:nPC]
    wPCA$rotation <- svd.w$v
    colnames(wPCA$x) <- paste("PC", 1:nPC, sep = "")
    wPCA$shapes <- lapply(1:ncol(wPCA$x),
                          function(x){shape.predictor(A, wPCA$x[,x],
                                                      min = min(wPCA$x[,x]),
                                                      max = max(wPCA$x[,x]))})
    names(wPCA$shapes) <- paste("shapes.PC", 1:nPC, sep = "")    
    class(wPCA) <- "gm.prcomp"
  }

  out <- list(rawPCA = rawPCA, phyloPCA = phyloPCA, phylomorphospace = phylomorphospace,
              alldataPCA = alldataPCA, wPCA = wPCA)
  
  class(out) = c("gm.prcomp", "list")
  attributes(out)$phy <- phy
  attributes(out)$Cov <- Cov
  attributes(out)$A <- A
  return(out)
}
