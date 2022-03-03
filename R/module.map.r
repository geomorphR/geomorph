#' Make a map of module covariances
#'
#' Function makes a plot of a covariance matrix arrange by modules, with colors that
#' measure the intensity of covariances within and between modules.
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, 
#' or a matrix (n x variables), or an already estimated covariance or correlation matrix.  
#' @param partition.gp A vector for an a priori, modular hypothesis, explaining
#' which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C). This is the same as partition.gp in, e.g., \code{\link{module.eigen}}.  
#' @param mincol The plot color for the minimum, non-0 covariance found.
#' @param maxcol The plot color for the maximum covariance found.
#' @param col0 The plot color for 0 (lacking) covariance.
#' @param bins The number  of color bins to us for covariance gradient.
#' @param near0 A logical value for whether covariances near 0 should be considered 0.  This is especially helpful
#' for geometric morphometric data that tend to be correlated because of generalized Procrustes analysis.  If TRUE,
#' more of the module map will be colored the same as col0, except the more intense covariances.
#' @param phy Optional argument to include a class \code{phylo} phylogenetic tree.  Tip labels must match data names.
#' This argument instructs the function to estimate a phylogenetic covariance matrix based on a Brownian motion model
#' of evolutionary divergence.  The Cov argument allows a user to define a hypothetical covariance matrix, if different 
#' than a BM model.
#' @param Cov Optional argument to include a hypothetical covariance matrix used for non-independence of observations.  
#' Row and column names must match data names.  If both a phy and Cov are provided, Cov will override phy.
#' @param transform. A logical argument for whether to use transformed residuals, if a phylogeny is provided.  If TRUE,
#' a GLS covariance matrix will be estimated; if FALSE, data will be centered on GLS mean but an OLS covariance matrix 
#' will be estimated.  The former is more representative of covariances independent of phylogeny;  the latter 
#' is more representative of dispersion in the tangent space.  See \code{\link{module.eigen}}.
#' @param label.vars A logical variable for whether to include variable names as row and column labels.
#' @param label.mods A logical variable for whether to numerically label modules.
#' @param label.cex Magnification of the variable labels in the plot.
#' @param borders A logical variable for whether to add borders (lines) around modules.
#' @export
#' @keywords analysis
#' @author Michael Collyer
#' @references Collyer et al. In review.
#' @seealso \code{\link{module.eigen}}, 
#' \code{\link{summary.module.eigen}}, \code{\link{plot.module.eigen}}
#' \code{\link{two.b.pls}}, \code{\link{modularity.test}}, 
#' \code{\link{phylo.integration}}, \code{\link{phylo.modularity}},
#' and \code{\link{compare.pls}}
#' @examples
#' 
#' # OLS vs GLS example
#' 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
#' # OLS approach
#' 
#' module.map(Y.gpa$coords, partition.gp = land.gps)
#' 
#' # GLS-centered approach.  This approach find the GLS mean but does not
#' # transform residuals.  Thus, it is produces OLS-like covariance matrices
#' # that are centered on the GLS mean rather than the OLS mean
#' 
#' module.map(Y.gpa$coords, partition.gp = land.gps,
#' phy = plethspecies$phy, transform. = FALSE)
#' 
#' # GLS-transformed approach.  This approach find the GLS mean and
#' # transforms residuals.  Thus, it produces GLS covariance matrices,
#' # from residuals rendered independent of phylogenetic covariances.
#' 
#' module.map(Y.gpa$coords, partition.gp = land.gps,
#' phy = plethspecies$phy, transform. = TRUE)
#' 
#' # Multi-module-hypothesis example
#' data(pupfish)
#' 
#' # Two-module hypothesis
#' hyp2 <- c(rep(1, 3), 2, rep(1, 5),
#' rep(2, 9), rep(1, 20), rep(2, 18))
#' 
#' # Three-module hypothesis
#' hyp3 <- c(rep(1, 3), 2, 1, rep(3, 4),
#' rep(2, 9), rep(1, 12), rep(3, 8), rep(2, 18))
#' 
#' par(mfrow = c(1, 2))
#' plot(pupfish$coords[,,1], pch = 21, bg = hyp2, asp = 1)
#' plot(pupfish$coords[,,1], pch = 21, bg = hyp3, asp = 1)
#' 
#' module.map(pupfish$coords, partition.gp = hyp2)
#' module.map(pupfish$coords, partition.gp = hyp3)
#' 
#' par(mfrow = c(1, 1))

module.map <- function(A, partition.gp, 
                       mincol = "gray", maxcol = "dark red", 
                       col0 = "white",
                       bins = 10, 
                       near0 = TRUE,
                       phy = NULL, Cov = NULL,
                       transform. = TRUE,
                       label.vars = FALSE,
                       label.mods = TRUE,
                       label.cex = 0.5,
                       borders = TRUE){
  
  
  if(any(is.na(A)))
    stop("\nData matrix contains missing values. Estimate these first (see 'estimate.missing').",
         call. = FALSE)  
  
  partition.gp <- as.factor(partition.gp)
  if (length(dim(A)) == 3){ 
    dims <- dim(A)
    p <- dims[1]
    k <- dims[2]
    
    if(length(partition.gp) != p) 
      stop("\nNot all landmarks are assigned to a partition.", call. = FALSE)
    
    gps <- as.factor(rep(partition.gp, k, each = k, length = p * k)) 
  }
  
  if (length(dim(A)) == 2){ 
    
    if(length(partition.gp) != ncol(A))
      stop("\nNot all variables are assigned to a partition.", call. = FALSE)
    
    gps <- as.factor(partition.gp) 
    
  }
  
  dims <- dim(A)
  n <- dims[1]
  p <- dims[2]
  
  if(n == p){
    if(round(sum(A - t(A)), 12) == 0)
      V <- A
  } else V <- get.VCV(A, phy, Cov, transform.)
  
  or.ord <- 1:nrow(V)
  df <- data.frame(gps = gps, or.ord = or.ord)
  df <- df[, order(c(1,2))]
  gps.ord <- df[, 1]
  
  V <- V[order(gps.ord, gps), order(gps.ord)] 
  
  ngps <- nlevels(gps.ord)
  ind.levels <- levels(gps.ord)
  
  plot(1:nrow(V), 1:nrow(V), type = "n", xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", bty = "n")
  colpal <- colorRampPalette(colors = c(mincol, maxcol))
  
  V <- abs(V)
  VV <- round(V, 6)
  diag(VV) <- 0
  if(all(VV == 0)) VV <- V
 
  col.check <- seq.int(min(VV), max(VV), length.out = bins)
  
  col.vec <- c(col0, colpal(bins - 1))
  
  if(near0) {
    l <- round(bins/4)
    col.vec[1:l] <- col0
  }
  
  nRow <- nrow(V)
  Cex <- 100 / nRow
  
  colors.id <-  lapply(1:length(V), function(x) which.min(abs(V[x] - col.check)))
  colors <- sapply(colors.id, function(x) col.vec[x])
  
  xid <- nRow:1
  xy <- cbind(rep(xid, nRow), rep(1:nRow, each = nRow))
  points(xy, pch = 22, col = colors, bg = colors, cex = Cex)
  gpn <- by(gps.ord, gps.ord, length)
  
  spots <- cumsum(gpn) - gpn /2
  spots <- cbind(spots, nrow(V) - spots)
  
  if(label.mods) text(spots, ind.levels, cex = 4*label.cex)
  
  if(borders) {
    for(i in 1:(length(gpn) - 1)) 
      abline(v = Cex + cumsum(gpn)[i])
    for(i in 1:(length(gpn) - 1)) 
      abline(h = nrow(V) - Cex - cumsum(gpn)[i])
  }

  hlab <- rownames(V)
  vlab <- rev(hlab)
  
  if(label.vars) {
    text(cbind(1:nrow(V), rep(-0.5, nrow(V))), hlab, cex = label.cex)
    text(cbind(1:nrow(V), rep(nrow(V) + 1, nrow(V))), hlab, cex = label.cex)
    text(cbind(rep(-0.3, nrow(V)), 1:nrow(V)), vlab, cex = label.cex)
    text(cbind(rep(nrow(V) + 1, nrow(V)), 1:nrow(V)), vlab, cex = label.cex)
    
  }
  
}