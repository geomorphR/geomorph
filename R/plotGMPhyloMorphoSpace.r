#' Plot phylogenetic tree and specimens in tangent space
#'
#' Function plots a phylogenetic tree and a set of Procrustes-aligned specimens in tangent space
#'
#' The function creates a plot of the principal dimensions of tangent space for a set of Procrustes-aligned 
#'   specimens. Default is a plot of PC axis 1 and 2. The phylogenetic tree for these specimens is superimposed in this plot revealing how shape 
#'   evolves (e.g., Rohlf 2002; Klingenberg and Gidaszewski 2010). The plot also displays the ancestral 
#'   states for each node of the phylogenetic tree (obtained from \code{\link[phytools]{fastAnc}}), whose values can optionally be returned. 
#'
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param labels A logical value indicating whether taxa labels (tips and ancestors) should be included
#' @param axis1 A value indicating which PC axis should be displayed as the X-axis (default = PC1)
#' @param axis2 A value indicating which PC axis should be displayed as the Y-axis (default = PC2)
#' @param axis3 Optional, a value indicating which PC axis should be displayed as the Z-axis (default = PC3)
#' @param ancStates A logical value indicating whether ancestral state values should be returned
#' @param plot.param A list of plotting parameters for the tips (t.bg, t.pch, t.cex), nodes (n.bg, n.pch, n.cex), branches (l.col, lwd), and taxa labels (txt.cex, txt.adj, txt.col)
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @return Function returns estimated ancestral states if {ancStates=TRUE}
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic 
#'   signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @references Rohlf, F. J. 2002. Geometric morphometrics and phylogeny. Pp.175-193 in N. Macleod, and 
#'   P. Forey, eds. Morphology, shape, and phylogeny. Taylor & Francis, London.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords)
#' plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords, 
#'                  plot.param=list(t.bg="blue",txt.col="red",n.cex=1))
plotGMPhyloMorphoSpace<-function(phy,A,labels=TRUE,ancStates=TRUE, axis1=1, axis2=2, axis3=NULL, 
                                 plot.param = list()){
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
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  N<-length(phy$tip.label)
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  p.p <- plot.param
    if(is.null(p.p$t.bg)) p.p$t.bg="black" ; if(is.null(p.p$t.pch)) p.p$t.pch=21
    if(is.null(p.p$t.cex)) p.p$t.cex=2 ; if(is.null(p.p$n.bg)) p.p$n.bg="white"
    if(is.null(p.p$n.pch)) p.p$n.pch=21 ; if(is.null(p.p$n.cex)) p.p$n.cex=1.25
    if(is.null(p.p$l.col)) p.p$l.col="black" ; if(is.null(p.p$lwd)) p.p$lwd=3
    if(is.null(p.p$txt.adj)) p.p$txt.adj=c(-.1,-.1) ; if(is.null(p.p$txt.col)) p.p$txt.col="black"
    if(is.null(p.p$txt.cex)) p.p$txt.cex=1
  x<-x[phy$tip.label, ]  
  names<-row.names(x)
  anc.states<-NULL
  for (i in 1:ncol(x)){
    options(warn=-1)  
    tmp <- as.vector(fastAnc(phy, x[, i]))
    anc.states<-cbind(anc.states,tmp)   }
  colnames(anc.states)<-NULL
  ## add labels to anc.states
  row.names(anc.states)<-1:length(tmp)
  all.data<-rbind(x,anc.states)  
  pcdata<-prcomp(all.data)$x 
  limits = function(x,s){ 
    r = range(x)
    rc=scale(r,scale=F)
    l=mean(r)+s*rc}
  # regular 2D phylomorphospace
  if(is.null(axis3)){
    if(labels==TRUE){
      plot(pcdata,type="n",xlim=limits(pcdata[,axis1],1.5),ylim=limits(pcdata[,axis2],1.5),asp=1) }
    if(labels==FALSE) {
      plot(pcdata,type="n",asp=1) }
    for (i in 1:nrow(phy$edge)){
      lines(pcdata[(phy$edge[i,]),axis1],pcdata[(phy$edge[i,]),axis2],type="l",col=p.p$l.col,lwd=p.p$lwd)
    }
    points(pcdata[1:N,], pch=p.p$t.pch, bg=p.p$t.bg, cex=p.p$t.cex)
    points(pcdata[(N+1):nrow(pcdata),], pch=p.p$n.pch, bg=p.p$n.bg, cex=p.p$n.cex)
    if(labels==TRUE){
      text(pcdata[,axis1],pcdata[,axis2],rownames(pcdata),col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj)
      }
  }
  # 3d phylomorphospace in rgl
  if(!is.null(axis3)){
    if(labels==TRUE){
      plot3d(pcdata,type="n",xlim=limits(pcdata[,axis1],1.5),
             ylim=limits(pcdata[,axis2],1.5),
             zlim=limits(pcdata[,axis3],1.5),asp=1,
             xlab= "PC 1", ylab= "PC 2", zlab="PC 3") }
    if(labels==FALSE) {
      plot3d(pcdata,type="n",asp=1) }
    points3d(pcdata[1:N,1], pcdata[1:N,2], pcdata[1:N,3],
             col= p.p$t.bg, size=p.p$t.cex*4)
    if(p.p$n.bg == "white"){ p.p$n.bg <- "grey"}
    points3d(pcdata[(N + 1):nrow(pcdata),1], pcdata[(N + 1):nrow(pcdata),2], 
             pcdata[(N + 1):nrow(pcdata),3], 
             col= p.p$n.bg, size=p.p$n.cex*4)
    for (i in 1:nrow(phy$edge)) {
      lines3d(pcdata[(phy$edge[i, ]), 1], pcdata[(phy$edge[i, ]), 2],pcdata[(phy$edge[i, ]), 3], 
              lwd=2)}
    
    if(labels==TRUE){
      text3d(pcdata[,axis1],pcdata[,axis2],pcdata[,axis3],rownames(pcdata),col=p.p$txt.col,cex=p.p$txt.cex,adj=p.p$txt.adj)
    }
  }
  if(ancStates==TRUE){ return(anc.states)  }
}