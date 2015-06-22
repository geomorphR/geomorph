#' Comparing rates of shape evolution on phylogenies
#'
#' Function calculates rates of shape evolution for two or more groups of species on a phylogeny from a set of Procrustes-aligned specimens
#'
#' The function compares rates of morphological evolution for two or more groups of species on a phylogeny, under a 
#'  Brownian motion model of evolution. It is assumed that the landmarks have previously been aligned 
#'  using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The approach is based on the distances
#'  between species in morphospace after phylogenetic transformation (Adams 2014). From the data the rate of shape evolution
#'  for each group is calculated, and a ratio of rates is obtained. If three or more groups of species are used, the ratio of 
#'  the maximum to minimum rate is used as a test statistic (see Adams 2014). Significance testing 
#'  is accomplished by phylogenetic simulation in which tips data are obtained under Brownian motion using a common 
#'  evolutionary rate pattern for all species on the phylogeny (i.e., a common evolutionary rate matrix). This procedure is
#'  more general, and retains the desirable statistical properties of earlier methods, but for a wider array of data types.  
#'  If three or more groups of species are used, pairwise p-values are also returned. A histogram of evolutionary rate ratios obtained 
#'  via phylogenetic simulation is presented, 
#'  with the observed value designated by an arrow in the plot. The function can be used to obtain a rate for the whole
#'  dataset of species by using a dummy group factor assigning all species to one group.
#'  
#'  This function can be used with univariate data (i.e. centroid size) if imported as matrix with rownames
#'  giving the taxa names.
#'
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param gp A factor array designating group membership
#' @param ShowPlot A logical value indicating whether or not the plot should be returned
#' @param iter Number of iterations for significance testing
#' @keywords analysis
#' @author Dean Adams & Emma Sherratt
#' @export
#' @return Function returns a list with the following components: 
#'   \item{sigma.d}{The phylogenetic evolutionary rate for all species on the phylogeny}
#'   \item{sigmad.all}{The phylogenetic evolutionary rate for each group of species on the phylogeny}
#'   \item{sigmad.ratio}{The ratio of maximum to minimum evolutionary rates}
#'   \item{pvalue}{The significance level of the observed ratio}
#'   \item{pairwise.pvalue}{Matrix of pairwise significance levels comparing each pair of rates}
#'   
#' @references Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
#'  shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'
#'  gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
#'  names(gp.end)<-plethspecies$phy$tip
#' 
#' #Calculate rates of shape
#' compare.evol.rates(plethspecies$phy,Y.gpa$coords,gp=gp.end,iter=49)
#' 
#' #Calculate rates of size
#' Csize <- matrix(Y.gpa$Csize, dimnames=list(names(Y.gpa$Csize))) # make matrix Csize with names
#' compare.evol.rates(plethspecies$phy,Csize,gp=gp.end,iter=49)
compare.evol.rates<-function(phy,A,gp,ShowPlot=TRUE,iter=999 ){
  if (length(dim(A))==3){ 
      if(is.null(dimnames(A)[[3]])){
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
      x<-two.d.array(A)}
  if (length(dim(A))==2){ 
      if(is.null(rownames(A))){
      stop("Data matrix does not include taxa names as dimnames for rows.")  }
      x<-A }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(!is.factor(gp)){
    stop("gp is not a factor.")}
  if (is.null(names(gp))){
    stop("Factor contains no names. Use names() to assign specimen names to group factor.")}
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  ntaxa<-length(phy$tip.label)
  N<-nrow(x)  
  if(N!=dim(x)[1]){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x), phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  p<-ncol(x)           
  sigma.d<-function(phy,x,N,gp){
    x<-prcomp(x)$x
    gp<-gp[rownames(x)]
    ngps<-nlevels(gp)
    gpsz<-table(gp)
    ones<-array(1,N)
    C<-vcv.phylo(phy)
    C<-C[rownames(x),rownames(x)]
    a.obs<-colSums(solve(C))%*%x/sum(solve(C))  
    eigC <- eigen(C)
    lambda <- zapsmall(eigC$values)
    if(any(lambda == 0)){
      warning("Singular phylogenetic covariance matrix. Proceed with caution")
      lambda = lambda[lambda > 0]
    }
    eigC.vect = eigC$vectors[,1:(length(lambda))]
    D.mat <- solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect)) 
    dist.adj<-as.matrix(dist(rbind((D.mat%*%(x-(ones%*%a.obs))),0))) 
    vec.d2<-dist.adj[N+1,1:N]^2
    sigma.d<-tapply(vec.d2,gp,sum)/gpsz/p   
    sigma.d.all<-sum(vec.d2)/N/p   
    if (ngps==1){ return(sigma.d.all) }
    if (ngps==2){sigma.d.rat<-max(sigma.d)/min(sigma.d)
                 return(list(sigma.all = sigma.d.all, ratio = sigma.d.rat, sigma.d.all=sigma.d)) }
    if(ngps>2){
      sigma.d.rat.gp<- array(0, dim = c(ngps, ngps))
      for (i in 1:(ngps - 1)) {
        for (j in 2:ngps) {
          tmp<-c(sigma.d[i],sigma.d[j])
          sigma.d.rat.gp[i, j] <-max(tmp)/min(tmp)
          diag(sigma.d.rat.gp) <- 0
          sigma.d.rat.gp[lower.tri(sigma.d.rat.gp)] <- 0 }
      }
        sigma.d.rat <- max(sigma.d.rat.gp)
        return(list(sigma.all = sigma.d.all, ratio = sigma.d.rat, sigma.d.all=sigma.d,
                    sigma.gp = sigma.d.rat.gp)) 
    }
  }
  if (nlevels(gp) == 1) {
    print("Single group. Sigma calculated for all specimens together.")
    sigmad.obs<-sigma.d(phy,x,ntaxa,gp) 
    return(list(sigma.d = sigmad.obs))
  }
  if (nlevels(gp) > 1) {
    sigmad.obs<-sigma.d(phy,x,ntaxa,gp) 
    ones<-array(1,N)
    C<-vcv.phylo(phy); C<-C[rownames(x),rownames(x)]
    a.obs<-colSums(solve(C))%*%x/sum(solve(C))  
    rate.mat<-t(x-ones%*%a.obs)%*%solve(C)%*%(x-ones%*%a.obs)/N
    x.sim<-sim.char(phy,rate.mat,nsim=iter) 
    sig.sim<-1
    if (nlevels(gp) > 2) {
      gp.sig.sim <- array(1, dim = c(dim(sigmad.obs$sigma.gp)[1], 
                                     dim(sigmad.obs$sigma.gp)[1]))}
    rate.val<-rep(0,iter)
    for(ii in 1:iter){
      sigmad.sim<-sigma.d(phy,x.sim[,,ii],ntaxa,gp)
      sig.sim<-ifelse(sigmad.sim$ratio>=sigmad.obs$ratio, sig.sim+1,sig.sim)
      rate.val[ii]<-sigmad.sim$ratio
      if (nlevels(gp) > 2) {
        gp.sig.sim <- ifelse(sigmad.sim$sigma.gp >= sigmad.obs$sigma.gp, 
                             gp.sig.sim + 1, gp.sig.sim)}
    }
    sig.sim<-sig.sim/(iter+1)
    if (nlevels(gp) > 2) {
      gp.sig.sim <- gp.sig.sim/(iter + 1)
      rownames(gp.sig.sim) <- colnames(gp.sig.sim) <- levels(gp)
      gp.sig.sim[lower.tri(gp.sig.sim)]<-NA
    }
    rate.val[iter+1]=sigmad.obs$ratio
    if(ShowPlot==TRUE){ 
      hist(rate.val,30,freq=TRUE,col="gray",xlab="SigmaD ratio")
      arrows(sigmad.obs$ratio,50,sigmad.obs$ratio,5,length=0.1,lwd=2)
      
    }
    if (nlevels(gp) > 2) {
      return(list(sigma.d = sigmad.obs$sigma.all, sigmad.all = sigmad.obs$sigma.d.all, 
                  sigmad.ratio = sigmad.obs$ratio, pvalue = sig.sim, pairwise.pvalue = gp.sig.sim))
    } else if (nlevels(gp) == 2) {
      return(list(sigma.d = sigmad.obs$sigma.all, sigmad.all = sigmad.obs$sigma.d.all, 
                  sigmad.ratio = sigmad.obs$ratio, pvalue = sig.sim))
    }
  }
}