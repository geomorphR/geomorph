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
#' the capability of performing analyses generally referred to as:
#' 
#' \itemize{
#' \item{\bold{PCA}}{  Traditional PCA based on OLS-centering and projection of data.}
#' \item{\bold{  Phylomorphospace}}{  Traditional PCA with estimated ancestral states and phylogenetic branches
#' projected into ordination plots.}
#' \item{\bold{phyloPCA}}{  PCA based on GLS-centering and projection of data.  Also possible to 
#' project ancestral states into plots. Note that if transformed GLS-residuals are used for projection, the ancestral states
#' might not appear logical, as the projection is independent of phylogeny.  With OLS-centering, a phyloPCA as described by 
#' Revell (2009) is produced.}
#' \item{\bold{PaCA}}{  Phylogenetically-aligned component analysis.  Data are aligned to an axis of greatest phylogenetic
#' signal rather than axis of greatest dispersion.  This analysis can use either OLS- or GLS-centering and projection.  
#' Phylogenetic signal is strongest in the first few components of the OLS approach.  This analysis will make little sense with 
#' GLS-centering and projection of transformed residuals, since phylogenetic signal is removed the transformed data.  
#' See Collyer and Adams (2020) for more details. For greater flexibility for type of residuals and projection of trees, 
#' use \code{\link{ordinate}}.  See Collyer and Adams (2020) for details.}
#' }
#' 
#' \itemize{
#' \item{\bold{phy}}{  Whether a phylogeny and estimated ancestral states are considered in plots.}
#' \item{\bold{align.to.phy}}{  Whether components are aligned to phylogenetic signal (rather than principal axes).}
#' \item{\bold{GLS}}{  Whether to use GLS-centering and estimation of covariance matrix.}
#' \item{\bold{transform}} {Whether to transform GLS-residuals (making them independent of phylogeny and an orthogonal
#' projection from the transformed data space, as opposed to an oblique projection from the untransformed data space).}
#' }
#' 
#' PLOTTING: Contrary to previous geomorph implementations, gm.prcomp does not produce plots. 
#' For plotting options of gm.prcomp class objects combine \code{\link{plot.gm.prcomp}} and 
#' \code{\link{picknplot.shape}} following the examples below. 
#' 
#' SUMMARY STATISTICS: For principal component plots, the traditional statistics to summarize the analysis include
#' eigenvalues (variance by component), proportion of variance by component, and cumulative proportion of variance. 
#' When data are aligned to a phylogenetic covariance matrix, the statistics are less straightforward.  A summary of
#' of such an analysis (performed with \code{\link{summary.gm.prcomp}}) will produce these additional statistics:
#'
#' \itemize{
#' \item{\bold{Singular Value}}{  Rather than eigenvalues, the singular values from singular value decomposition of the 
#' cross-product of the scaled phylogenetic covariance matrix and the data.}
#' \item{\bold{Proportion of Covariance}}{  Each component's singular value divided by the sum of singular values.  The cumulative
#' proportion is also returned.  Note that these values do not explain the amount of covariance between phylogeny and data, but
#' explain the distribution of the covariance.  Large proportions can be misleading.}
#' \item{\bold{RV by Component}}{  The partial RV statistic by component.  Cumulative values are also returned.  The sum of partial
#' RVs is Escoffier's RV statistic, which measures the amount of covariation between phylogeny and data.  Caution should
#' be used in interpreting these values, which can vary with the number of observations and number of variables.  However,
#' the RV is more reliable than proportion of singular value for interpretation of the strength of linear association for 
#' phylogenetically-aligned components.}
#' \item{\bold{Tips or Ancestors Dispersion}}{  The variances of points by component for tip data and estimated ancestral
#' character states, after projection.  These values will differ from variances from PCA with GLS estimation, as the "Importance
#' of Components" weights variances by phylogenetic covariances.  Dispersion statistics correspond to the amount of scatter in
#' plots of component scores.}
#' }
#' 
#' 
#' NOTE: The \code{\link{plot.gm.prcomp}} function performs the same plotting that was previously 
#' possible with \code{plotTangentSpace} and \code{plotGMPhyloMorphoSpace}, which
#' have now been deprecated.
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of aligned specimens.  
#' Alternatively, this can be an n x p matrix of any data, but output will not contain information about shapes.
#' @param phy An optional phylogenetic tree of class phylo 
#' @param align.to.phy An optional argument for whether \bold{PaCA} (if TRUE) should be performed
#' @param GLS Whether GLS-centering and covariance estimation should be used (rather than OLS).
#' @param transform A logical value to indicate if transformed residuals should be projected.  This is only applicable if 
#' GLS = TRUE.  If TRUE, an orthogonal projection of transformed data is made; if FALSE an oblique projection of untransformed 
#' data is made.
#' @param ... Other arguments passed to \code{\link{ordinate}} and \code{\link{scale}}.  The most common
#' arguments are scale., tol, and rank.
#' @return An object of class "gm.prcomp" contains a list of results for each of the PCA approaches implemented.
#' Each of these lists includes the following components:
#' \item{x}{Component scores for all specimens.}
#' \item{anc.x}{Component scores for the ancestors on the phylogeny.}
#' \item{d}{The singular values of the decomposed VCV matrix.}
#' \item{rotation}{The matrix of variable loadings, i.e. the eigenvectors of the decomposed matrix.}
#' \item{shapes}{A list with the shape coordinates of the extreme ends of all PC axes.}
#' \item{ancestors}{The matrix of estimated ancestral shapes, if a phylogeny is used.}
#' \item{anc.var}{The variances among ancestor scores, by component, if a phylogeny is used.}
#' @export
#' @keywords visualization
#' @seealso \code{\link{plot.gm.prcomp}}
#' @seealso \code{\link{picknplot.shape}}
#' @seealso \code{\link{ordinate}}{ A more bare-bones ordination function on which gm.prcomp depends.}
#' @references Collyer, M.L and D.C. Adams, 2020. Phylogenetically-aligned component analysis. 
#' Methods in Ecology and Evolution (in review).
#' @references Revell, L. J. (2009). Size-correction and principal components for interspecific 
#' comparative studies. Evolution, 63: 3258-3268.
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
#'  # This is the same as the method described by Revell (2009)
#'  phylo.PCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, GLS = TRUE)
#'  summary(phylo.PCA)
#'  plot(phylo.PCA, phylo = TRUE, main = "phylo PCA")
#'  
#'  ### Phylogenetic PCA - PCA based on GLS-centering and transformed projection
#'  # This produces a PCA independent of phylogeny
#'  phylo.tPCA <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
#'  GLS = TRUE, transform = TRUE)
#'  summary(phylo.tPCA)
#'  plot(phylo.tPCA, phylo = TRUE, main = "phylo PCA")
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
#'  # GLS method (rotation of Phylogenetic PCA with transformed data)
#'  PaCA.gls <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy, 
#'  align.to.phy = TRUE, GLS = TRUE, transform = TRUE)
#'  summary(PaCA.gls)
#'  plot(PaCA.gls, phylo = TRUE, main = "PaCA using GLS and transformed projection")
#'  
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
#'  ### 3D plot with a phylogeny and time on the z-axis
#'  plot(PCA.w.phylo, time.plot = TRUE)
#'  plot(PCA.w.phylo, time.plot = TRUE, bg = "red", phylo.par = list(tip.labels = TRUE, 
#'  tip.txt.cex = 2, edge.color = "blue", edge.width = 2))
#'  


gm.prcomp <- function (A, phy = NULL, align.to.phy = FALSE,
                       GLS = FALSE, transform = FALSE, ...) {
  
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
    if(is.null(rownames(Y))) {
      warning("Shape dataset does not include species names. Assuming the order of data matches phy$tip.label")
      rownames(Y) <- phy$tip.label
    }

    ancY <- anc.BM(phy, Y)
    if(is.null(rownames(Y))) {
      rownames(Y) <- phy$tip.label
      cat("Warning: Data are not labeled so it is assumed they are in the same order as tree names.\n")
    }
    phy.mat <- phylo.mat(Y, phy)
    C <- phy.mat$C
    if(align.to.phy) ord.args$A <- C
    if(GLS) ord.args$Cov <- C

  }
  
  ord.args$transform. = transform

  out <- do.call(ordinate, ord.args)
  if(out$alignment == "A") out$alignment <- "phy"
  
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
    out$x <- as.matrix(out$x)
    ancs <- as.matrix(anc.BM(phy, out$x))
    out$anc.x <- ancs
    out$anc.var <- apply(ancs, 2, var)
  }
  if(GLS) {
    Pcov <- phy.mat$D.mat
    Pcov <- Pcov[rownames(C), rownames(C)]
    out$Pcov <- Pcov
  } else out$Pcov <- NULL
  
  if(!is.null(rownames(Y))) out$x <- out$x[rownames(as.matrix(Y)), ]
  
  class(out) <- c("gm.prcomp", class(out))
  out
}
