#' Principal and phylogenetically-aligned components analysis of shape data
#'
#' Function performs principal components analysis (PCA) or 
#' phylogenetically-aligned components (PaCA) on Procrustes shape coordinates.
#'
#' The function performs a series of ordinations, taking into account, phylogeny, if desired.
#' There are two main types of ordinations: principal components analysis (PCA) and phylogenetically-
#' aligned components analysis (PaCA).  Both of these have two variants: centering and projection
#' via ordinary least-squares (OLS) or via generalized least-squares (GLS).  The name,
#' "gm.prcomp", references that this function performs much like \code{\link{prcomp}}, in terms
#' of arguments and output, but this function is quite a bit more diverse.  This function has
#' the capability of performing analyses genereally referred to as:
#' 
#' \itemize{
#' \item{\bold{PCA}}{  Traditional PCA based on OLS-centering and projection of data.}
#' \item{\bold{  Phylomorphospace}}{  Traditional PCA with estimated ancestral states and phylogenetic branches
#' projected into ordination plots.}
#' \item{\bold{phyloPCA}}{  PCA based on GLS-centering and projection of data.  Also possible to 
#' project ancestral states into plots.}
#' \item{\bold{PaCA}}{  Alignment of components to maximum phylogenetic signal rather than
#' maximum variation.  This analysis can use either OLS- or GLS-centering and projection.  Phylogenetic
#' signal is strongest in the first few components.  See Collyer and Adams (in review) for more details.}
#' }
#' 
#' The examples provided will illustrate each of these methods, but the following arguments are essential for directing
#' the function:
#' 
#' \itemize{
#' \item{\bold{phy}}{  Whether a phylogeny and estimated ancestral states are considered in plots.}
#' \item{\bold{align.to.phy}}{  Whether components are aligned to phylogenetic signal (rather than principal axes).}
#' \item{\bold{GLS}}{  Whether to use PGLS residuals.  Note that doing so creates an oblique projection rather than 
#' an orthogonal projection of data, but the origin of the plot will be the tree root rather than the center of gravity.}
#' }
#' 
#' PLOTTING: Contrary to previous geomorph implementations, gm.prcomp does not produce plots. 
#' For plotting options of gm.prcomp class objects combine \code{\link{plot.gm.prcomp}} and 
#' \code{\link{picknplot.shape}} following the examples below. 
#' 
#' NOTE: The \code{\link{plot.gm.prcomp}} function performs the same plotting that was previously 
#' posible with \code{\link{plotTangentSpace}} and \code{\link{plotGMPhyloMorphoSpace}}, which have now been 
#' deprecated.
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of aligned specimens.  
#' Alternatively, this can be an n x p matrix of any data, but output will not conatin information about shapes.
#' @param phy An optional phylogenetic tree of class phylo 
#' @param align.to.phy An optional argument for whether \bold{PaCA} (if TRUE) should be performed
#' @param GLS An optional argument for whether GLS should be used for centering and projecting data.
#' @param ... Other arguments passed to \code{\link{ordinate}} and \code{\link{scale}}.  The most common
#' arguments are scale., tol, and rank.
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
#' @seealso \code{\link{ordinate}}{ A more bare-bones ordination function on which gm.prcomp depends.}
#' @author Antigoni Kaliontzopoulou, Michael Collyer, & Dean Adams
#' @examples
#'  data(plethspecies) 
#'  Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
#'  
#'  ###  Traditional PCA 
#'  PCA <- gm.prcomp(Y.gpa$coords)
#'  summary(PCA)
#'  plot(PCA, main = "PCA")
#'  
#'  ### Phylomorphospace - PCA with phylogeny (result is same as above, 
#'  ### but with estimated ancestral states projected into plot)
#'  PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#'  summary(PCA.w.phylo)
#'  plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo")
#'  
#'  ### Phylogenetic PCA - PCA based on GLS-centering and projection
#'  phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE)
#'  summary(phylo.PCA)
#'  plot(phylo.PCA, phylo = TRUE, main = "phylo PCA")
#'  
#'  ### PaCA - Alignment of data to physlogenetic signal rather than axis of 
#'  ### greatest variation, like in PCA
#'  
#'  # OLS method (rotation of PCA)
#'  PaCA.ols <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, align.to.phy = TRUE)
#'  summary(PaCA.ols)
#'  plot(PaCA.ols, phylo = TRUE, main = "PaCA using OLS")
#'  
#'  # GLS method (rotation of Phylogenetic PCA)
#'  PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
#'  align.to.phy = TRUE, GLS = TRUE)
#'  summary(PaCA.gls)
#'  plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS")
#'  
#'  ### Advanced Plotting
#'  gps <- as.factor(c(rep("gp1", 5), rep("gp2", 4))) # Two random groups
#'  par(mar=c(2, 2, 2, 2))
#'  plot(PaCA.ols, pch=22, cex = 1.5, bg = gps, phylo = TRUE) # Modify options as desired
#'  #  Add things as desired using standard R plotting
#'  text(par()$usr[1], 0.1*par()$usr[3], labels = "PC1 - 45.64%", pos = 4, font = 2)
#'  text(0, 0.95*par()$usr[4], labels = "PC2 - 18.80%", pos = 4, font = 2)
#'  legend("topleft", pch=22, pt.bg = unique(gps), legend = levels(gps))
#'  


gm.prcomp <- function (A, phy = NULL, align.to.phy = FALSE,
                       GLS = FALSE, ...) {
  
  if(is.array(A)) {
    
    dims <- dim(A)
    if(length(dims) == 3) { 
      
      if(any(is.na(A)))
        stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').\n",
             call. = FALSE) 
      p <- dims[1]; k <- dims[2]; n <- dims[3]
      
      } else if(length(dims) == 2) {
      
      n <- dims[1]; k <- NULL; p <- dims[2]
    }
    
  } else {
      A <- try(as.matrix(A), silent = TRUE)
      if(inherits(A, "try-error"))
        stop("Data not of a form coercible to matrix or array.\n", call. = FALSE)
      dims <- dim(A)
      n <- dims[1]; k <- NULL; p <- dims[2]
    }
  
  ord.args <- list(...)
  if(!is.null(k)) {
    Y <- two.d.array(A)
    ord.args$Y <- Y
  } else ord.args$Y <- Y <- A
  

  if(!is.null(phy)){
    if (!inherits(phy, "phylo"))
      stop("Tree must be of class 'phylo.'\n", call. = FALSE)
    N <- length(phy$tip.label); Nnode <- phy$Nnode
    if(N != n)
      stop("Number of taxa in data matrix and tree are not equal.\n", call. = FALSE)
    if(is.null(rownames(Y))) 
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")

    ancY <- anc.BM(phy, Y)
    C <- fast.phy.vcv(phy)
    if(!is.null(rownames(Y))) C <- C[rownames(Y), rownames(Y)]
    if(align.to.phy) ord.args$A <- C
    if(GLS) ord.args$Cov <- C
    ord.args$newdata <- ancY

}

  out <- do.call(ordinate, ord.args)
  if(!is.null(phy)) 
    names(out)[[which(names(out) == "xn")]] <- "anc.x"
  
  names(out)[[which(names(out) == "rot")]] <- "rotation"
  if(!is.null(k)) {
    
    out$A <- A
    out$shapes <- lapply(1:ncol(out$x),  
                         function(x){shape.predictor(A, out$x[,x], 
                                                     min = min(out$x[,x]),
                                                     max = max(out$x[,x]))})
    names(out$shapes) <- paste("shapes.comp", 1:length(out$d), sep = "")
  }

  if(!is.null(phy)) {
    out$ancestors <- ancY
    out$phy <- phy
  }
  class(out) <- c("gm.prcomp", class(out))
  out
}
