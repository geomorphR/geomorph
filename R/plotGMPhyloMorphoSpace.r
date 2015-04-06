#' Plot phylogenetic tree and specimens in tangent space
#'
#' Function plots a phylogenetic tree and a set of Procrustes-aligned specimens in tangent space
#'
#' The function creates a plot of the first two dimensions of tangent space for a set of Procrustes-aligned 
#'   specimens. The phylogenetic tree for these specimens is superimposed in this plot revealing how shape 
#'   evolves (e.g., Rohlf 2002; Klingenberg and Gidaszewski 2010). The plot also displays the ancestral 
#'   states for each node of the phylogenetic tree (obtained from \code{\link[phytools]{fastAnc}}), whose values can optionally be returned. 
#'
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param labels A logical value indicating whether taxa labels (tips and ancestors) should be included
#' @param ancStates A logical value indicating whether ancestral state values should be returned
#' @export
#' @keywords visualization
#' @author Dean Adams & Emma Sherratt
#' @return Function returns estimated ancestral states if {ancStates=TRUE}
#' @references Klingenberg, C. P., and N. A. Gidaszewski. 2010. Testing and quantifying phylogenetic 
#'   signals and homoplasy in morphometric data. Syst. Biol. 59:245-261.
#' @references Rohlf, F. J. 2002. Geometric morphometrics and phylogeny. Pp. 175'193 in N. Macleod, and 
#'   P. Forey, eds. Morphology, shape, and phylogeny. Taylor & Francis, London.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#' plotGMPhyloMorphoSpace(plethspecies$phy,Y.gpa$coords)
plotGMPhyloMorphoSpace<-function(phy,A,labels=TRUE,ancStates=TRUE){
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
  if(labels==TRUE){
    plot(pcdata,type="n",xlim=limits(pcdata[,1],1.5),ylim=limits(pcdata[,2],1.5),asp=1) }
  if(labels==FALSE) {
    plot(pcdata,type="n",asp=1) }
  for (i in 1:nrow(phy$edge)){
    lines(pcdata[(phy$edge[i,]),1],pcdata[(phy$edge[i,]),2],type="l",pch=21,col="black",lwd=3)
  }
  points(pcdata[1:N,],pch=21,bg="black",cex=2)
  points(pcdata[(N+1):nrow(pcdata),],pch=21,bg="white",cex=1.25)
  if(labels==TRUE){
    text(pcdata[,1],pcdata[,2],rownames(pcdata),cex=0.75,adj=c(-.1,-.1))}
  if(ancStates==TRUE){ return(anc.states)  }
}