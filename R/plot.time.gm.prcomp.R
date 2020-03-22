#' Create a 3d phylomorphospace plot of a \code{\link{gm.prcomp}} object with a phylogeny associated to it,
#' with time displayed on the z-axis.
#' 
#' This function replaces the 3d plotting functionality of plotGMPhyloMorphoSpace, now deprecated.
#' 
#' @param x A \code{\link{gm.prcomp}} object
#' @param axis1 A numeric value indicating which PC axis should be displayed as the X-axis
#' @param axis2 A numeric value indicating which PC axis should be displayed as the Y-axis
#' @param tip.labels A logical value indicating whether taxa labels (tips) should be included
#' @param node.labels A logical value indicating whether node labels (ancestors) should be included
#' @param anc.states A logical value indicating whether points for the nodes should be included
#' @param plot.param A list of plotting parameters for the tips (tip.bg, tip.cex), nodes (node.bg, node.cex),
#' edges (edge.colour, edge.width), and labels of the tips (txt.cex, txt.col, txt.adj) and nodes 
#' (node.txt.cex, node.txt.col, node.txt.adj)
#' @keywords visualization
#' @author Antigoni Kaliontzopoulou
#' @seealso \code{\link{gm.prcomp}} 
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#' 
#' PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#' 
#' plot.time.gm.prcomp(x)
#' plot.time.gm.prcomp(x, node.labels = F, anc.states = F, 
#' plot.param = list(tip.bg = "red", tip.cex = 3, 
#' edge.colour = "blue", edge.width = 2, txt.cex = 2))

plot.time.gm.prcomp <- function(x, axis1 = 1, axis2 = 2, 
                                tip.labels = TRUE, node.labels = TRUE, anc.states = TRUE, 
                                plot.param = list(tip.bg = "black", tip.cex = 1,
                                                  node.bg = "grey", node.cex = 1,
                                                  edge.color = "black", edge.width = 1,
                                                  txt.cex = 1, txt.col = "black", txt.adj = c(-0.1,-0.1),
                                                  node.txt.cex = 1, node.txt.col = "grey", node.txt.adj = c(-0.1, -0.1))){
  
  if(!inherits(x, "gm.prcomp")) stop("x must be an object of class gm.prcomp")   
  
  if(is.null(x$phy)) stop("x must include a phylogeny for plotting")  
  
  phy <- x$phy
  phy.pcdata <- rbind(x$x[x$phy$tip.label,], x$anc.x)
  phy.pcdata <- as.matrix(phy.pcdata[, c(axis1, axis2)])
  
  zaxis <- getNodeDepth(phy)
  zaxis <- abs(zaxis - max(zaxis))
  
  limits = function(x,s){ 
    r = range(x)
    rc = scale(r, scale=F)
    l = mean(r) + s*rc 
    }
  
  p.p <- plot.param
  if(is.null(p.p$tip.bg)) p.p$tip.bg <- "black"
  if(is.null(p.p$tip.cex)) p.p$tip.cex <- 1
  if(is.null(p.p$node.bg)) p.p$node.bg <- "grey"
  if(is.null(p.p$node.cex)) p.p$node.cex <- 1
  if(is.null(p.p$edge.color)) p.p$edge.color <- "black"
  if(is.null(p.p$edge.width)) p.p$edge.width <- 1
  if(is.null(p.p$txt.cex)) p.p$txt.cex <- 1
  if(is.null(p.p$txt.col)) p.p$txt.col <- "black"
  if(is.null(p.p$txt.adj)) p.p$txt.adj <- c(-0.1, -0.1)
  if(is.null(p.p$node.txt.cex)) p.p$node.txt.cex <- 1
  if(is.null(p.p$node.txt.col)) p.p$node.txt.col <- "grey"
  if(is.null(p.p$node.txt.adj)) p.p$node.txt.adj <- c(-0.1, -0.1)
  
  view3d(phi=90, fov=30)
  plot3d(phy.pcdata[,1], phy.pcdata[,2], zaxis, type="n", 
         xlim = limits(phy.pcdata[,1], 1.5),
         ylim = limits(phy.pcdata[,2], 1.5),
         zlim = c(max(zaxis), min(zaxis)),
         asp = c(1,1,1),
         xlab= paste("PC", axis1), ylab = paste("PC", axis2), zlab = "Time")
  
  for (i in 1:nrow(phy$edge)) {
    lines3d(phy.pcdata[(phy$edge[i, ]), 1], phy.pcdata[(phy$edge[i, ]), 2], zaxis[(phy$edge[i, ])], 
            col = p.p$edge.color, lwd = p.p$edge.width)}
  
  points3d(phy.pcdata[1:Ntip(phy), 1], phy.pcdata[1:Ntip(phy), 2], zaxis[1:Ntip(phy)],
           col = p.p$tip.bg, size = p.p$tip.cex*4)
  
  if(anc.states){
    points3d(phy.pcdata[(Ntip(phy) + 1):nrow(phy.pcdata), 1], 
             phy.pcdata[(Ntip(phy) + 1):nrow(phy.pcdata), 2], 
             zaxis[(Ntip(phy) + 1):nrow(phy.pcdata)], 
             col = p.p$node.bg, size = p.p$node.cex*4)
  }
  
  if(tip.labels){
    text3d(phy.pcdata[1:Ntip(phy), 1], phy.pcdata[1:Ntip(phy), 2], zaxis[1:Ntip(phy)], 
           rownames(phy.pcdata)[1:Ntip(phy)],
           col = p.p$txt.col, cex = p.p$txt.cex, adj = p.p$txt.adj) }
  
  if(node.labels){
    text3d(phy.pcdata[(Ntip(phy) + 1):nrow(phy.pcdata), 1], phy.pcdata[(Ntip(phy) + 1):nrow(phy.pcdata), 2],
           zaxis[(Ntip(phy) + 1):nrow(phy.pcdata)], rownames(phy.pcdata)[(Ntip(phy) + 1):nrow(phy.pcdata)],
           col = p.p$node.txt.col, cex = p.p$node.txt.cex, adj = p.p$node.txt.adj)}
}
