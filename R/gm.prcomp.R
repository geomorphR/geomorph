#' Principal components analysis of shape data
#'
#' Function performs raw or weighted PCA on superimposed shape coordinates  
#'
#' The function performs a principal components analysis of shape variation, with the possibility
#' of weighing the analysis using a phylogenetic tree, or a variance-covariance matrix, 
#' to incorporate the expected covariance among dependent observations 
#' (e.g. due to spatial autocorrelation). At present combined analysis allowing the use of both 
#' a phylogeny and another covariance matrix is not implemented.
#' If a phylogeny is provided and phylo.pca = TRUE, observations are phylogenetically adjusted to 
#' calculate a phylogenetic PCA (Revell, 2009). 
#' If phylo.pca = FALSE, ancestral shapes are first estimated for the phylogeny nodes and PCA is then calculated
#' on the combined matrix of raw and ancestor data, producing the data necessary for a phylomorphospace plot.
#' 
#' PLOTTING: Contrary to previous geomorph implementations, gm.prcomp does not produce plots. 
#' For plotting options of gm.prcomp class objects see \code{\link{plot.gm.prcomp}} and the 
#' examples below.
#' 
#'
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of aligned specimens
#' @param phy An optional phylogenetic tree of class phylo - see \code{\link{read.tree}} in library ape
#' @param phylo.pca A logical value indicating whether phylogenetic PCA (TRUE) 
#' or the phylomorphospace approach (FALSE) should be used. See also details.
#' @param Cov An optional covariance matrix for weighting. See also details.
#' @param ... Arguments passed on to \code{\link{prcomp}}.  By default, \code{\link{gm.prcomp}}
#' will attempt to remove redundant axes (eigenvalues effectively 0).  To override this, adjust the 
#' argument, tol, from \code{\link{prcomp}}.
#' @return An object of class "gm.prcomp" is a list with the following components:
#' \item{pc.summary}{A table summarizing the percent variation explained by each pc axis, equivalent to summary of \code{\link{prcomp}}.}
#' \item{pc.scores}{The set of principal component scores for all specimens.}
#' \item{pc.shapes}{A list with the shape coordinates of the extreme ends of all PC axes, e.g. $PC1min}
#' \item{sdev}{The standard deviations of the principal components (i.e., the square roots of the eigenvalues of the 
#' covariance/correlation matrix, as per \code{\link{prcomp}}.}
#' \item{rotation}{The matrix of variable loadings, as per \code{\link{prcomp}}.}
#' \item{anc.states}{The matrix of estimated ancestral shapes, if a phylogeny is used.}
#' \item{anc.pcscores}{The matrix of principal component scores for the phylogeny nodes.}
#' @export
#' @keywords visualization
#' @author Antigoni Kaliontzopoulou
#' @references Revell, L. J. (2009) Size-correction and principal components for interspecific comparative studies. 
#' Evolution, 63, 3258-3268.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' 
#' ### Raw data PCA
#' pleth.raw <- gm.prcomp(Y.gpa$coords)
#' summary(pleth.raw)
#' 
#' # Plotting
#' gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
#' par(mgp = c(2.5, 0.5, 0))
#' plot(pleth.raw, pch=22, cex = 1.5, xlab = "My PCA - axis 1", bg = gps,
#'     font.lab = 2, cex.lab = 2) # Modify options as desired
#  Add things as desired using standard R plotting
#' segments(0.95*par()$usr[1], 0, 0.95*par()$usr[2], 0, lty = 2, lwd = 1)
#' segments(0, 0.95*par()$usr[3], 0, 0.95*par()$usr[4], lty = 2, lwd = 1)
#' legend("topright", pch=22, pt.bg = unique(gps), legend = levels(gps), cex = 2)
#' 
#' ### Phylogenetic PCA
#' pleth.ppca <- gm.prcomp(Y.gpa$coords, plethspecies$phy, phylo.pca = TRUE)
#' summary(pleth.ppca) 
#' 
#' # Plotting
#' par(mgp = c(2, 0.5, 0))
#' plot(pleth.ppca, phylo = TRUE, cex = 1.5, pch = 22, bg = gps, cex.lab = 2, 
#'      font.lab = 2, xlim = c(-0.007, 0.017),
#'      phylo.par = list(edge.color = "grey", edge.width = 2,
#'      node.bg = "black", node.pch = 22, node.cex = 0.5))
#' text(pleth.ppca$pc.scores, labels = labels(pleth.ppca$pc.scores)[[1]],
#'      pos = 2, font = 4) 
#' text(pleth.ppca$anc.pcscores, labels = labels(pleth.ppca$anc.pcscores)[[1]],
#'      adj = c(-0.1, -0.1), font = 2) 
#'      
#' ### Phylomorphospace
#' pleth.phylomorpho <- gm.prcomp(Y.gpa$coords, plethspecies$phy)
#' summary(pleth.phylomorpho)
#' 
#' # Plotting
#' plot(pleth.phylomorpho, phylo = TRUE, cex = 2, pch = 22, bg = gps, 
#'      phylo.par = list(edge.color = "blue", edge.width = 2, edge.lty = 2,
#'      node.cex = 0)) # Supresses plotting of nodes
#' text(pleth.phylomorpho$pc.scores, labels = labels(pleth.phylomorpho$pc.scores)[[1]],
#'      pos = 3, font = 4) 
#'
#' # same plot but change the PC axes type
#' 
#' plot(pleth.phylomorpho, phylo = TRUE, cex = 2, pch = 22, bg = gps, 
#'      phylo.par = list(edge.color = "blue", edge.width = 2, edge.lty = 2,
#'      node.cex = 0), axes = FALSE) # Supresses plotting of nodes and axes
#'  abline(h = 0, lty = 3, cex = 0.7) # subtle axes
#'  abline(v = 0, lty = 3, cex = 0.7)
#' text(pleth.phylomorpho$pc.scores, labels = labels(pleth.phylomorpho$pc.scores)[[1]],
#'      pos = 3, font = 4) 

gm.prcomp <- function (A, phy = NULL, phylo.pca = FALSE, Cov = NULL, ...){
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').") }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }

  dots <- list(...)
  retx <- dots$retx
  if(is.null(retx)) retx <- TRUE
  scale. <- dots$scale.
  if(is.null(scale.)) scale. <- FALSE
  center <- dots$center
  if(is.null(center)) center <- TRUE
  tol <- dots$tol
  p <- dim(A)[1]; k <- dim(A)[2]; n <- dim(A)[3]
  ref <- mshape(A)
  x <- scale(two.d.array(A), center = center, scale = scale.)

  if(is.null(phy) & phylo.pca == T){
    stop("To perform phylogenetic pca, please provide a phylogeny.")
  }
  
  if(!is.null(phy)){
    if (!inherits(phy, "phylo"))
      stop("Tree must be of class 'phylo.'")
    if (!is.binary.tree(phy)) 
      stop("Tree is not fully bifurcating (consider 'multi2di' in ape).")
    if(!is.null(Cov)){
      stop("Method not implemented for weighting with BOTH a phylogeny and another covariance matrix.")
    }
    N <- length(phy$tip.label); Nnode <- phy$Nnode
    if(N!=n)
      stop("Number of taxa in data matrix and tree are not equal.")
    if(is.null(rownames(x))) {
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")
    } else x <- x[phy$tip.label, ]
    anc <- shape.ace(x, phy)
    
    if(phylo.pca == T){
      phy.parts <- phylo.mat(x, phy)
      invC <- phy.parts$invC; D.mat <- phy.parts$D.mat
      one <- matrix(1, nrow(x)); I <- diag(1, nrow(x)) 
      Ptrans <- D.mat%*%(I-one%*%crossprod(one, invC)/sum(invC))
      x <- Ptrans%*%x
      center = apply(x, 2, mean)
      anc <- scale(anc, center = anc[1,], scale = F)
    } else {
      x <- rbind(x, anc)
      center <- anc[1,]
    }
  } else {
    if(!is.null(Cov)){
      cov.parts <- cov.mat(x, Cov)
      invC <- cov.parts$invC; D.mat <- cov.parts$D.mat
      one <- matrix(1, nrow(x)); I <- diag(1, nrow(x)) 
      Ptrans <- D.mat%*%(I-one%*%crossprod(one, invC)/sum(invC))
      x <- Ptrans%*%x
      center = apply(x, 2, mean)
    }
  }
  
  if(is.null(tol)){
    d <- prcomp(x)$sdev^2
    cd <-cumsum(d)/sum(d)
    cd <- length(which(cd < 1)) 
    if(length(cd) < length(d)) cd <- cd + 1
    tol <- max(c(d[cd]/d[1], 0.005))
  }

  pc.res <- prcomp(x, retx = retx, tol = tol, center = center)
  pcdata <- pc.res$x
  shapes <- shape.names <- NULL
  for(i in 1:ncol(pcdata)){
    pcaxis.min <- min(pcdata[, i]) ; pcaxis.max <- max(pcdata[, i])
    pc.min <- pc.max <- rep(0, dim(pcdata)[2])
    pc.min[i] <- pcaxis.min ; pc.max[i] <- pcaxis.max
    pc.min <- as.matrix(pc.min %*% (t(pc.res$rotation))) + as.vector(t(ref))
    pc.max <- as.matrix(pc.max %*% (t(pc.res$rotation))) + as.vector(t(ref))
    shapes <- rbind(shapes,pc.min, pc.max)
    shape.names <- c(shape.names,paste("PC",i,"min", sep=""),paste("PC",i,"max", sep=""))
  }
  shapes <- arrayspecs(shapes,p, k)
  shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[,,x])
  names(shapes) <- shape.names
  
  # OUTPUT
  if(is.null(phy) & is.null(Cov)) meth <- "Raw data PCA"
  if(!is.null(phy) & phylo.pca == FALSE) meth <- "Phylomorphospace"
  if(!is.null(phy) & phylo.pca == TRUE) meth <- "Phylogenetic PCA"
  if(!is.null(Cov)) meth <- "Cov-weighted PCA"
  
  out <- list(pc.summary = summary(pc.res), pc.scores = pcdata[1:n, ], pc.shapes = shapes, 
              sdev = pc.res$sdev, rotation = pc.res$rotation, anc.states = NULL, anc.pcscores = NULL)
  
  if(meth == "Phylomorphospace") {
    out$anc.states <- anc
    out$anc.pcscores <- pcdata[(N+1):nrow(pcdata),]
  }
  
  if(meth == "Phylogenetic PCA"){
    out$anc.states <- anc
    out$anc.pcscores <- anc%*%pc.res$rotation
  }
  
  class(out) = "gm.prcomp"
  attributes(out)$method <- meth
  attributes(out)$phy <- phy
  attributes(out)$Adata <- A
  return(out)
}
