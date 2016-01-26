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
#'  evolutionary rate pattern for all species on the phylogeny. Specifically, the common evolutionary rate matrix for all
#'  species is used, with the multi-dimensional rate used along the diagonal elements (see Denton and Adams 2015). This procedure is
#'  more general than the original simulation procedure, and retains the desirable statistical properties of earlier methods, 
#'  and under a wider array of data types.  
#'  If three or more groups of species are used, pairwise p-values are also calculated. The function can be used to obtain a 
#'  rate for the whole dataset of species by using a dummy group factor assigning all species to one group.
#'  
#'  This function can be used with univariate data (i.e. centroid size) if imported as matrix with rownames
#'  giving the taxa names.
#'  
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{compare.evol.rates}}.
#'  The generic function, \code{\link{plot}}, produces a histogram of random rate-ratios associated with
#'  the resampling procedure.
#'
#'  \subsection{Notes for geomorph 3.0}{ 
#' Compared to older versions of geomorph, the order of input variables has changed, so that it is consistent with other functions
#' in the program.  Additionally, for 3 or more groups, the pairwise p-values are found in the output object.}
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param gp A factor array designating group membership
#' @param iter Number of iterations for significance testing
#' @keywords analysis
#' @author Dean Adams & Emma Sherratt
#' @export
#' @return An object of class "evolrate" returns a list with the following components: 
#'   \item{sigma.d.ratio}{The ratio of maximum to minimum evolutionary rates.}
#'   \item{P.value}{The significance level of the observed ratio.}
#'   \item{sigma.d.gp}{The phylogenetic evolutionary rate for each group of species on the phylogeny.}
#'   \item{random.sigma}{The sigma values found in random permutations of the resampling procedure.}
#'   \item{permutations}{The number of random permutations used.}
#'   
#' @references Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
#'  shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.
#' @references Denton, J.S.S., and D.C. Adams. 2015. A new phylogenetic test for comparing 
#' multiple high-dimensional evolutionary rates suggests interplay of evolutionary rates and 
#' modularity in lanternfishes (Myctophiformes; Myctophidae). Evolution. 69:2425-2440.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'  gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
#'  names(gp.end)<-plethspecies$phy$tip
#' 
#' ER<-compare.evol.rates(A=Y.gpa$coords, phy=plethspecies$phy,gp=gp.end,iter=999)
#' summary(ER)
#' plot(ER)
compare.evol.rates<-function(A,phy,gp,iter=999 ){
  gp<-as.factor(gp)
  if (length(dim(A))==3){ 
      if(is.null(dimnames(A)[[3]])){
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")  }
      x<-two.d.array(A)}
  if (length(dim(A))==2){ 
      if(is.null(rownames(A))){
      stop("Data matrix does not include taxa names as dimnames for rows.")  }
      x<-A }
  if (is.vector(A)== TRUE){ 
    if(is.null(names(A))){
      stop("Data vector does not include taxa names as names.")  }
    x<-as.matrix(A) }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
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
  gp<-gp[rownames(x)]
  phy.parts<-phylo.mat(x,phy)
  invC<-phy.parts$invC; D.mat<-phy.parts$D.mat;C = phy.parts$C
  sigma.obs<-sigma.d(x,invC,D.mat,gp)
  rate.mat<-sigma.obs$R
  diag(rate.mat)<-sigma.obs$sigma.d.all
  rate.mat<-matrix(nearPD(rate.mat,corr=FALSE)$mat,nrow=ncol(rate.mat),ncol=ncol(rate.mat))
  x.sim<-sim.char(phy=phy,par=rate.mat,nsim=iter,model="BM") 
  sigma.rand <- sapply(1:(iter), function(j) sigma.d(as.matrix(x.sim[,,j]),invC,D.mat,gp))
  random.sigma<-c(sigma.obs$sigma.d.ratio,as.vector(unlist(sigma.rand[1,])))
  if(nlevels(gp)>1){
    p.val <- pval(random.sigma)
    p.val.mat<-NULL
    if(nlevels(gp)==2) p.val.mat<-p.val
    if(nlevels(gp)>2){
      ratio.vals<-matrix(NA,nrow=(iter+1),ncol=length(unlist(sigma.obs[4])))
      ratio.vals[1,]<-as.vector(sigma.obs$sigma.d.gp.ratio)
      for(i in 1:iter) ratio.vals[i+1,]<-as.vector(unlist(sigma.rand[4,][[i]]))
      tmp.p.val.mat <- sapply(1:ncol(ratio.vals), function(j){ pval(ratio.vals[,j])})
      p.val.mat<-dist(matrix(0,length(tmp.p.val.mat)))
      for(i in 1:length(p.val.mat)) p.val.mat[[i]] <- tmp.p.val.mat[i]
    }    
  }
  if(nlevels(gp)==1){ 
    out <- list(sigma.d.all = sigma.obs$sigma.d.all,
                Ngroups = nlevels(gp))
    
    class(out) <- "evolrate1"
    }
  if(nlevels(gp)>1){
    out <- list(sigma.d.ratio = sigma.obs$sigma.d.ratio, P.value=p.val,
                sigma.d.all = sigma.obs$sigma.d.all,
                sigma.d.gp = sigma.obs$sigma.d.gp,
                sigma.d.gp.ratio = sigma.obs$sigma.d.gp.ratio,
                pairwise.pvalue = p.val.mat, Ngroups = nlevels(gp),
                groups = levels(gp),
                random.sigma = random.sigma, permutations=iter+1, 
                call = match.call())
    
    class(out) <- "evolrate"
  }
  out 
}