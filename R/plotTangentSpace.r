#' Plot specimens in tangent space
#'
#' Function plots a set of Procrustes shape variables in tangent space along their principal component axes

#' \subsection{Notes for geomorph 3.3.0 and subsequent versions}{
#'  This function is defunct and replaced with \code{\link{gm.prcomp}}.  The functions,
#'  \code{\link{summary.gm.prcomp}}, \code{\link{plot.gm.prcomp}}, and \code{\link{picknplot.shape}} 
#'  allow summaries, ammendable plots, and warpgrids, for points anywhere in a plot, respectively, can
#'  be used to accomplish any task of plotTangentSpace, with greater flexibility.
#'  }
#' 
plotTangentSpace<-function (){
  .Defunct("gm.prcomp", package = "geomorph", 
              msg = "plotTangentSpace has been removed from geomorph.
              A combination of gm.prcomp, summary.prcomp, and plot.gm.prcomp has much greater flexibility.
              One can also use picknplot.shape to generate warpgrids, anywhere in a plot.")
}
