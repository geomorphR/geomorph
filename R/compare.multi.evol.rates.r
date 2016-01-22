#' Comparing rates of shape evolution among traits on phylogenies
#'
#' Function calculates rates of shape evolution for two or more multi-dimensional traits on a 
#' phylogeny from a set of Procrustes-aligned specimens
#'
#' The function compares rates of morphological evolution for two or more multi-dimensional traits
#' on a phylogeny, under a Brownian motion model of evolution following the procedure of Denton and 
#' Adams (2015). It is assumed that the landmarks for all traits have previously been aligned using
#' Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The approach calculates
#' multivariate evolutionary rates found from the distances between species in morphospace after 
#' phylogenetic transformation (sensu Adams 2014). From the data the rate of shape evolution for 
#' each multi-dimensional trait is calculated, and a ratio of rates is obtained. If three or more 
#' traits are used, the ratio of the maximum to minimum rate is used as a test statistic (see 
#' Denton and Adams 2015). Significance testing is accomplished by phylogenetic simulation in 
#' which tips data are obtained under Brownian motion using a an evolutionary rate matrix 
#' for all traits, which contains a common rate for all trait dimensions (Denton and Adams 2015).
#' If three or more traits are used, pairwise p-values are 
#' also returned. A histogram of evolutionary rate ratios obtained via phylogenetic simulation 
#' is presented, with the observed value designated by an arrow in the plot. 
#' 
#' The shape data may be input as either a 3D array (p x k x n) containing GPA-aligned coordinates 
#' for a set of species, or as a matrix (n x [p x k]) whose rows correspond to each species. In 
#' both cases, species names must be provided as rownames (for a matrix) or as the names of the 
#' third dimension of the array. Landmark  groups for each trait are then specified by a factor
#' array designating which landmark belongs to which trait. Additionally, if the method is to be 
#' used with other data (i.e., a set of length measurements), the input A should be a matrix 
#' of n rows of species and p columns of variables. In this case, the grouping factor should 
#' have each variable assigned to a trait group. 
#' 
#' Comparisons of evolutionary rates between traits may be accomplished in one of two ways. First, 
#' if the traits are are part of a single shape that was subjected to a single Procrustes 
#' superimposition (i.e., they are subsets of landmarks in the configuration), then the procedure
#' is performed without alteration as described above. However, if the shapes are derived from 
#' different structures (shapes) that were superimposed separately, then the estimates of the rates must 
#' take the difference in the number of trait dimensions into account (see discussion in Denton and
#' Adams 2015). This option is identified by selecting Subset = FALSE.
#'  
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param gp A factor array designating group membership
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param Subset A logical value indicating whether or not the traits are subsets from a single 
#' landmark configuration (default is TRUE)
#' @param ShowPlot A logical value indicating whether or not the plot should be returned
#' @param iter Number of iterations for significance testing
#' @keywords analysis
#' @author Dean Adams
#' @export
#' @return Function returns a list with the following components: 
#'   \item{rates.all}{The phylogenetic evolutionary rates for each trait}
#'   \item{rate.ratio}{The ratio of maximum to minimum evolutionary rates}
#'   \item{pvalue}{The significance level of the observed rate ratio}
#'   \item{pvalue.gps}{Matrix of pairwise significance levels comparing each pair of rates}
#'   
#' @references Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
#'  shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.
#' @references Denton, J.S.S., and D.C. Adams. 2015. A new phylogenetic test for comparing 
#' multiple high-dimensional evolutionary rates suggests interplay of evolutionary rates and 
#' modularity in lanternfishes (Myctophiformes; Myctophidae). Evolution. 69: doi:10.1111/evo.12743
#' @examples
#' 
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#' land.gp<-c("A","A","A","A","A","B","B","B","B","B","B")  #mandible and cranium subsets
#'
#' compare.multi.evol.rates(Y.gpa$coords,land.gp,plethspecies$phy,iter=99)
compare.multi.evol.rates<-function(A,gp,phy,Subset=TRUE,ShowPlot=TRUE,iter=999){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")}
  gp<-as.factor(gp)
  if (length(dim(A))==3){ x<-two.d.array(A)
                          p<-dim(A)[1]; k<-dim(A)[2]
                          if(length(gp)!=p){stop("Not all landmarks are assigned to a partition.")}
                            gps<-as.factor(rep(gp,k,each = k, length=p*k))  }
  if (length(dim(A))==2){ x<-A
                          if(length(gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
                          gps<-gp  }
  ngps<-nlevels(gp)
  if(ngps==1){stop("Only one shape assigned.")}
  ntaxa<-length(phy$tip.label)
  N<-nrow(x) 
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if(is.null(rownames(x))){
    stop("Data matrix does not include taxa names.")  }
  if(N!=ntaxa){
    stop("Number of taxa in data matrix and tree are not not equal.")  }
  if(length(match(rownames(x),phy$tip.label))!=N) 
    stop("Data matrix missing some taxa present on the tree.")
  if(length(match(phy$tip.label,rownames(x)))!=N) 
    stop("Tree missing some taxa in the data matrix.")
  sigma.d<-function(phy,x,Subset){
    ones<-array(1,N)
    p<-ncol(x)
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
    sigma<-sum(vec.d2)/N/p  
    if(Subset==FALSE){sigma<-sum(vec.d2)/N}
    return(sigma = sigma)
  }
  rate.global<-sigma.d(phy,x,Subset)
  rate.gps<-array(NA,ngps)
  for (i in 1:nlevels(gps)){rate.gps[i]<-sigma.d(phy,x[,which(gps==levels(gps)[i])],Subset)}
  rate.ratio<-max(rate.gps)/min(rate.gps)
  if(ngps>2){
    sig.rate.gps<-array(1,dim=c(ngps,ngps))
    rate.ratio.gps<- array(0, dim = c(ngps, ngps))
    for (i in 1:(ngps - 1)) {
      for (j in 2:ngps) { tmp<-c(rate.gps[i],rate.gps[j])
        rate.ratio.gps[i, j]<-rate.ratio.gps[j,i]  <-max(tmp)/min(tmp)
        diag(rate.ratio.gps) <- 0
      }
    }
  }
  ones<-array(1,N)
  p<-ncol(x)
  C<-vcv.phylo(phy)
  C<-C[rownames(x),rownames(x)]
  a.obs<-colSums(solve(C))%*%x/sum(solve(C))  
  rate.mat<-t(x-ones%*%a.obs)%*%solve(C)%*%(x-ones%*%a.obs)/ntaxa
  diag(rate.mat)<-rate.global
  rate.mat<-matrix(nearPD(rate.mat,corr=FALSE)$mat,nrow=p,ncol=p)
  x.sim<-sim.char(phy,rate.mat,nsim=iter) 
  sig.rate<-1
  rate.val<-rep(0,iter)
  for(ii in 1:iter){
    rate.gps.r<-array(NA,ngps)
    for (i in 1:nlevels(gps)){rate.gps.r[i]<-sigma.d(phy,x.sim[,which(gps==levels(gps)[i]),ii],Subset)}
    rate.ratio.r<-max(rate.gps.r)/min(rate.gps.r)
    if(ngps>2){
      rate.ratio.gps.r<- array(0, dim = c(ngps, ngps))
      for (i in 1:(ngps - 1)) {
        for (j in 2:ngps) {    tmp<-c(rate.gps.r[i],rate.gps.r[j])
          rate.ratio.gps.r[i, j] <-rate.ratio.gps.r[j, i]<-max(tmp)/min(tmp)
          diag(rate.ratio.gps.r) <- 0
        }
      }
    }
    sig.rate<-ifelse(rate.ratio.r>=rate.ratio, sig.rate+1,sig.rate)
    if(ngps>2){ sig.rate.gps<-ifelse(rate.ratio.gps.r>=rate.ratio.gps, sig.rate.gps+1,sig.rate.gps) }
    rate.val[ii]<-rate.ratio.r
  }
  sig.rate<-sig.rate/(iter+1)
  if(ngps>2){  sig.rate.gps<-sig.rate.gps/(iter+1)}
  rate.val[iter+1]=rate.ratio
  if(ShowPlot==TRUE){ 
  hist(rate.val,30,freq=TRUE,col="gray",xlab="SigmaD ratio")
  arrows(rate.ratio,50,rate.ratio,5,length=0.1,lwd=2) }
  if(ngps==2){
    return(list(rates.all = rate.gps, rate.ratio = rate.ratio,pvalue=sig.rate))}
  if(ngps>2){
    return(list(rates.all = rate.gps, rate.ratio = rate.ratio,pvalue=sig.rate,pvalue.gps=sig.rate.gps))}
}
