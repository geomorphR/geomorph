#' Plot phylogenetic tree and specimens in tangent space
#'
#' Function plots a phylogenetic tree and a set of Procrustes shape variables in tangent space
#'
#' The function creates a plot of the principal dimensions of tangent space for a set of Procrustes shape variables. 
#' Default is a plot of PC axis 1 and 2. The phylogenetic tree for these specimens is superimposed in this plot revealing how shape 
#'   evolves (e.g., Rohlf 2002; Klingenberg and Gidaszewski 2010). The plot also displays the ancestral 
#'   states for each node of the phylogenetic tree (analogous to  fastAnc from phytools), whose values can optionally be returned. 
#'   If a tree with branch lengths scaled by time is used, with the option zaxis = "time", the function plots a 3D phylomorphospace, with internal nodes positioned along the Z-axis scaled 
#'   to time (a.k.a. Chronophylomorphospace, Sakamoto & Ruta 2012).
#'   
#' \subsection{Notes for geomorph 3.3.0 and subsequent versions}{
#'  This function is defunct and replaced with \code{\link{gm.prcomp}}.  The functions,
#'  \code{\link{summary.gm.prcomp}}, \code{\link{plot.gm.prcomp}}, and \code{\link{picknplot.shape}} 
#'  allow summaries, ammendable plots, and warpgrids, for points anywhere in a plot, respectively, can
#'  be used to accomplish any task of plotGMPhyloMorphospace, with greater flexibility.
#'  }
#'
#' @param phy A phylogenetic tree of {class phylo} 
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing Procrustes shape variables for a set of specimens
#' @param tip.labels A logical value indicating whether taxa labels (tips) should be included
#' @param node.labels A logical value indicating whether node labels (ancestors) should be included
#' @param xaxis A numeric value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param yaxis A numeric value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param zaxis Optional, a numeric value indicating which PC axis should be displayed as the Z-axis (e.g. PC3) or if zaxis="time", 
#' internal nodes are plotted along the Z-axis relative to time
#' @param ancStates Either a logical value indicating whether ancestral state values should be returned, or a matrix of ancestral states
#' @param plot.param A list of plotting parameters for the tips (t.bg, t.pch, t.cex), nodes (n.bg, n.pch, n.cex), 
#' branches (l.col, lwd), taxa labels (txt.cex, txt.adj, txt.col) and node labels (n.txt.cex, n.txt.adj, n.txt.col)
#' @param shadow A logical value indicating whether a 2D phylomorphospace should be plotted at the base when zaxis="time"
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @seealso \code{\link{gm.prcomp}}
#' @return Function returns estimated ancestral states if {ancStates=TRUE}
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic 
#'   signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @references Rohlf, F. J. 2002. Geometric morphometrics and phylogeny. Pp.175-193 in N. Macleod, and
#'   P. Forey, eds. Morphology, shape, and phylogeny. Taylor & Francis, London.
#' @references Sakamoto, M. and Ruta, M. 2012. Convergence and Divergence in the Evolution of Cat
#' Skulls: Temporal and Spatial Patterns of Morphological Diversity. PLoSONE 7(7): e39752.

plotGMPhyloMorphoSpace <- function(phy,A,tip.labels=TRUE,node.labels=TRUE,ancStates=TRUE, xaxis=1, yaxis=2, zaxis=NULL, plot.param = list(), shadow=FALSE){
  .Defunct("gm.prcomp", package = "geomorph", 
              msg = "plotTangentSpace has been removed from geomorph.
              A combination of gm.prcomp, summary.prcomp, and plot.gm.prcomp has much greater flexibility.
              One can also use picknplot.shape to generate warpgrids, anywhere in a plot.")
}
