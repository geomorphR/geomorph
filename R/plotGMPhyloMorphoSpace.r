#' Plot phylogenetic tree and specimens in tangent space
#'
#' Function plots a phylogenetic tree and a set of Procrustes-aligned specimens in tangent space
#'
#' The function creates a plot of the principal dimensions of tangent space for a set of Procrustes-aligned 
#'   specimens. Default is a plot of PC axis 1 and 2. The phylogenetic tree for these specimens is superimposed in this plot revealing how shape 
#'   evolves (e.g., Rohlf 2002; Klingenberg and Gidaszewski 2010). The plot also displays the ancestral 
#'   states for each node of the phylogenetic tree (analogous to from \code{\link[phytools]{fastAnc}} from phytools), whose values can optionally be returned. 
#'   If a tree with branch lengths scaled by time is used, with the option zaxis = "time", the function plots a 3D phylomorphospace, with internal nodes positioned along the Z-axis scaled 
#'   to time (a.k.a. Chronophylomorphospace, Sakamoto & Ruta 2012).
#'
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param tip.labels A logical value indicating whether taxa labels (tips) should be included
#' @param node.labels A logical value indicating whether node labels (ancestors) should be included
#' @param xaxis A numeric value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param yaxis A numeric value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param zaxis Optional, a numeric value indicating which PC axis should be displayed as the Z-axis (e.g. PC3) or if zaxis="time", 
#' internal nodes are plotted along the Z-axis relative to time
#' @param ancStates Either a logical value indicating whether ancestral state values should be returned, or a matrix of ancestral states (i.e. calculated with \code{\link[phytools]{fastAnc}} or \code{\link[ape]{ace}})
#' @param plot.param A list of plotting parameters for the tips (t.bg, t.pch, t.cex), nodes (n.bg, n.pch, n.cex), 
#' branches (l.col, lwd), taxa labels (txt.cex, txt.adj, txt.col) and node labels (n.txt.cex, n.txt.adj, n.txt.col)
#' @param shadow A logical value indicating whether a 2D phylomorphospace should be plotted at the base when zaxis="time"
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @return Function returns estimated ancestral states if {ancStates=TRUE}
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic 
#'   signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @references Rohlf, F. J. 2002. Geometric morphometrics and phylogeny. Pp.175-193 in N. Macleod, and
#'   P. Forey, eds. Morphology, shape, and phylogeny. Taylor & Francis, London.
#' @references Sakamoto, M. and Ruta, M. 2012. Convergence and Divergence in the Evolution of Cat
#' Skulls: Temporal and Spatial Patterns of Morphological Diversity. PLoSONE 7(7): e39752.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords)
#' plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords, 
#'                  plot.param=list(t.bg="blue",txt.col="red",n.cex=1))
#' #NOTE: 3D plot also available: plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords, zaxis= "time",
#' #                 plot.param=list(n.cex=2, n.bg="blue"), shadow=TRUE)
plotGMPhyloMorphoSpace<-function(phy,A,tip.labels=TRUE,node.labels=TRUE,ancStates=TRUE, xaxis=1, yaxis=2, zaxis=NULL, plot.param = list(), shadow=FALSE){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A))==3){ 
    if(is.null(dimnames(A)[[3]])){
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
    x<-two.d.array(A)}
  if (length(dim(A))==2){ 
    if(is.null(dimnames(A)[[1]])){
      stop("Data matrix does not include taxa names as dimnames for rows.")  }
    x<-A }
  if (!inherits(phy, "phylo"))
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating (consider 'multi2di' in ape.")
  N<-length(phy$tip.label)
  Nnode <- phy$Nnode
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  x<-x[phy$tip.label, ]  
  anc.states<-NULL   #follows fastAnc in phytools
  for (i in 1:ncol(x)){
    x1<-x[,i]
    tmp <- vector()
    for (j in 1:Nnode + N) {
      a <- multi2di(root(phy, node = j))
      tmp[j - N] <- ace(x1, a, method = "pic")$ace[1]
    }
    anc.states<-cbind(anc.states,tmp)   }
  colnames(anc.states)<-NULL
  row.names(anc.states)<-1:length(tmp)
  all.data<-rbind(x,anc.states)  
  pcdata<-prcomp(all.data)$x 
  pcdata<-pcdata-matrix(rep(pcdata[(N+1),],nrow(pcdata)), nrow=nrow(pcdata),byrow=T)  #phylogenetic mean adjustment
#plotting  
  p.p <- plot.param
  if(is.null(p.p$t.bg)) p.p$t.bg="black" ; if(is.null(p.p$t.pch)) p.p$t.pch=21
  if(is.null(p.p$t.cex)) p.p$t.cex=2 ; if(is.null(p.p$n.bg)) p.p$n.bg="white"
  if(is.null(p.p$n.pch)) p.p$n.pch=21 ; if(is.null(p.p$n.cex)) p.p$n.cex=1.25
  if(is.null(p.p$l.col)) p.p$l.col="black" ; if(is.null(p.p$lwd)) p.p$lwd=3
  if(is.null(p.p$txt.adj)) p.p$txt.adj=c(-.1,-.1) ; if(is.null(p.p$txt.col)) p.p$txt.col="black"
  if(is.null(p.p$txt.cex)) p.p$txt.cex=1 ; if(is.null(p.p$n.txt.adj)) p.p$n.txt.adj=c(-.1,-.1) 
  if(is.null(p.p$n.txt.col)) p.p$n.txt.col="black" ; if(is.null(p.p$n.txt.cex)) p.p$n.txt.cex=0.6
  limits = function(x,s){ 
    r = range(x)
    rc=scale(r,scale=F)
    l=mean(r)+s*rc}
  # regular 2D phylomorphospace
  if(is.null(zaxis)){
    if(tip.labels==TRUE){
      plot(pcdata[,xaxis],pcdata[,yaxis],type="n",xlim=limits(pcdata[,xaxis],1.5),ylim=limits(pcdata[,yaxis],1.5),asp=1,
           xlab = paste("PC", xaxis), ylab = paste("PC", yaxis)) }
    if(tip.labels==FALSE) {
      plot(pcdata[,xaxis],pcdata[,yaxis],type="n",asp=1, xlab = paste("PC", xaxis), ylab = paste("PC", yaxis)) }
    for (i in 1:nrow(phy$edge)){
      lines(pcdata[(phy$edge[i,]),xaxis],pcdata[(phy$edge[i,]),yaxis],type="l",col=p.p$l.col,lwd=p.p$lwd)
    }
    points(pcdata[1:N,xaxis], pcdata[1:N,yaxis],pch=p.p$t.pch, bg=p.p$t.bg, cex=p.p$t.cex)
    points(pcdata[(N+1):nrow(pcdata),xaxis], pcdata[(N+1):nrow(pcdata),yaxis],pch=p.p$n.pch, bg=p.p$n.bg, cex=p.p$n.cex)
    if(tip.labels==TRUE){
      text(pcdata[1:N,xaxis],pcdata[1:N,yaxis],rownames(pcdata)[1:N],
           col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj)}
    if(node.labels==TRUE){
      text(pcdata[(N + 1):nrow(pcdata),xaxis],pcdata[(N + 1):nrow(pcdata),yaxis],rownames(pcdata)[(N + 1):nrow(pcdata)],
           col=p.p$n.txt.col,cex=p.p$n.txt.cex,adj=p.p$n.txt.adj)}
  }
  # 3d phylomorphospace in rgl
  if(is.numeric(zaxis)){
    plot3d(pcdata[1:N,xaxis], pcdata[1:N,yaxis], pcdata[1:N,zaxis],type="s",xlim=limits(pcdata[,xaxis],1.5),
             ylim=limits(pcdata[,yaxis],1.5), zlim=limits(pcdata[,zaxis],1.5), asp=1,
             xlab= paste("PC",xaxis), ylab= paste("PC",yaxis), zlab=paste("PC",zaxis),
             col= p.p$t.bg, size=p.p$t.cex)
    if(p.p$n.bg == "white"){ p.p$n.bg <- "grey"}
    points3d(pcdata[(N + 1):nrow(pcdata),xaxis], pcdata[(N + 1):nrow(pcdata),yaxis], 
             pcdata[(N + 1):nrow(pcdata),zaxis], 
             col= p.p$n.bg, size=p.p$n.cex*4)
    for (i in 1:nrow(phy$edge)) {
      lines3d(pcdata[(phy$edge[i, ]), xaxis], pcdata[(phy$edge[i, ]), yaxis],pcdata[(phy$edge[i, ]), zaxis], 
              col=p.p$l.col, lwd=p.p$lwd)}
    if(tip.labels==TRUE){
      text3d(pcdata[1:N,xaxis],pcdata[1:N,yaxis],pcdata[1:N,zaxis],rownames(pcdata)[1:N],
             col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj) }
    if(node.labels==TRUE){
      text3d(pcdata[(N + 1):nrow(pcdata),xaxis],pcdata[(N + 1):nrow(pcdata),yaxis],
             pcdata[(N + 1):nrow(pcdata),zaxis],rownames(pcdata)[(N + 1):nrow(pcdata)],
             col=p.p$n.txt.col,cex=p.p$n.txt.cex,adj=p.p$n.txt.adj) }
  }
  # 3d phylomorphospace in rgl with time on Z-axis
  if(is.character(zaxis)){
    zaxis <- node.depth.edgelength(phy)
    zaxis <- abs(node.depth.edgelength(phy) - max(zaxis))
    view3d(phi=90, fov=30)
    plot3d(pcdata[,xaxis],pcdata[,yaxis],zaxis,type="n",xlim=limits(pcdata[,xaxis],1.5),
             ylim=limits(pcdata[,yaxis],1.5),
             zlim=c(max(zaxis), min(zaxis)),
             asp=c(1,1,1),
           xlab= paste("PC",xaxis), ylab= paste("PC",yaxis), zlab="Time")
    points3d(pcdata[1:N,xaxis], pcdata[1:N,yaxis], zaxis[1:N],
             col= p.p$t.bg, size=p.p$t.cex*4)
    if(p.p$n.bg == "white"){ p.p$n.bg <- "grey"}
    points3d(pcdata[(N + 1):nrow(pcdata),xaxis], pcdata[(N + 1):nrow(pcdata),yaxis], 
             zaxis[(N + 1):nrow(pcdata)], 
             col= p.p$n.bg, size=p.p$n.cex*4)
    for (i in 1:nrow(phy$edge)) {
      lines3d(pcdata[(phy$edge[i, ]), xaxis], pcdata[(phy$edge[i, ]), yaxis],zaxis[(phy$edge[i, ])], 
              col=p.p$l.col, lwd=p.p$lwd)}
    if(tip.labels==TRUE){
      text3d(pcdata[1:N,xaxis],pcdata[1:N,yaxis],zaxis[1:N],rownames(pcdata)[1:N],
             col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj) }
    if(node.labels==TRUE){
      text3d(pcdata[(N + 1):nrow(pcdata),xaxis],pcdata[(N + 1):nrow(pcdata),yaxis],
             zaxis[(N + 1):nrow(pcdata)],rownames(pcdata)[(N + 1):nrow(pcdata)],
             col=p.p$n.txt.col,cex=p.p$n.txt.cex,adj=p.p$n.txt.adj) }
      if(shadow==TRUE){
        #plot shadow version at base
        points3d(pcdata[1:N,xaxis], pcdata[1:N,yaxis], max(zaxis),
                 col= p.p$t.bg, size=p.p$t.cex*4, alpha = 0.5)
        points3d(pcdata[(N + 1):nrow(pcdata),xaxis], pcdata[(N + 1):nrow(pcdata),yaxis], 
                 max(zaxis), col= p.p$n.bg, size=p.p$n.cex*4,alpha = 0.5)
        for (i in 1:nrow(phy$edge)) {
          lines3d(pcdata[(phy$edge[i, ]), xaxis], pcdata[(phy$edge[i, ]), yaxis],max(zaxis), 
                  col=p.p$l.col, lwd=p.p$lwd, alpha = 0.5)}
    }
  }
  if(ancStates==TRUE){ return(anc.states)  }
}