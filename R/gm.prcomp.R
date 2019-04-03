#' Principal components analysis of shape data
#'
#' Function performs PCA on Procrustes shape coordinates  
#'
#' The function performs a principal components analysis of shape variation.
#' It also allows a phylomorphospace approach, where the user can provide a phylogeny
#' and obtain a graph with estimated ancestral shape values and the tree projected
#' into PCA space.
#' 
#' PLOTTING: Contrary to previous geomorph implementations, gm.prcomp does not produce plots. 
#' For plotting options of gm.prcomp class objects combine \code{\link{plot.gm.prcomp}} and 
#' \code{\link{picknplot.shape}} following the examples below. Note that trying to plot the full result
#' of gm.prcomp will give an error. Choose a specific PCA method to be plotted, by pointing to one of the 
#' components of the list returned by gm.prcomp.
#' 
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of aligned specimens
#' @param phy An optional phylogenetic tree of class phylo 
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
#' @examples
#'  data(plethspecies) 
#'  Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
#'  
#'  ### PCA 
#'  pleth.raw <- gm.prcomp(Y.gpa$coords)
#'  summary(pleth.raw)
#'  
#'  ### PCA with phylogeny (result is same as above, but with additional components)
#'  pleth.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#'  summary(pleth.phylo)
#'  
#'  #### Plotting
#'  plot(pleth.raw)
#'  gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
#'  par(mar=c(2, 2, 2, 2))
#'  plot(pleth.raw, pch=22, cex = 1.5, bg = gps) # Modify options as desired
#'  #  Add things as desired using standard R plotting
#'  text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", pos = 4, font = 2)
#'  text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
#'  legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))
#'  
#'  ### Phylomorphospace plot
#'  plot(pleth.phylo, pch=21, bg=1:nrow(pleth.phylo$x), phylo = TRUE,
#'       phylo.par = list(edge.color="grey", node.cex=0)); title(main="phylomorphospace")
#'  text(pleth.phylo$x, labels = labels(pleth.phylo$x)[[1]],
#'       pos = 2, font = 4) 
#'  text(pleth.phylo$anc.x, labels = 1:nrow(pleth.phylo$anc.x),
#'       adj = c(-0.1, -0.1), font = 2) 
#'  
#'  ### Visualize shape variation using picknplot.shape Because picknplot requires 
#'  ### user decisions, the following example
#'  ### is not run (but can be with removal of #).
#'  ### For detailed options, see the picknplot help file
#'  # picknplot.shape(plot(pleth.phylo))
#'  
#' 

gm.prcomp <- function (A, phy = NULL, ...) {
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }

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
    if(N != n)
      stop("Number of taxa in data matrix and tree are not equal.")
    if(is.null(rownames(Y))) {
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")
    } else {
      A <- A[,, phy$tip.label]
      Y <- Y[phy$tip.label, ]
    }
    anc.raw <- anc.BM(phy, two.d.array(A))
    ancY <- anc.BM(phy, Y)
    all.dt <- rbind(Y, ancY)
    vcv.a <- crossprod(all.dt) / (nrow(all.dt)-1)
    svd.a <- svd(vcv.a)
  }
  
  ### Output
  rawPCA <- phylomorphospace <- NULL
  
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
    phylomorphospace$phy <- phy
    }

  if(is.null(phy)) {
    out <- rawPCA
  } else {
    out <- phylomorphospace
  }
  
  out$A <- A

  class(out) = c("gm.prcomp")
  
  out
}
