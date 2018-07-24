#' Comparing net rates of shape evolution on phylogenies
#'
#' Function calculates net rates of shape evolution for two or more groups of species on a phylogeny from a set of Procrustes shape variables
#'
#' The function compares net rates of morphological evolution for two or more groups of species on a phylogeny, under a 
#'  Brownian motion model of evolution. It is assumed that the landmarks have previously been aligned 
#'  using Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The approach is based on the outer-product matrix of  
#'  between species differences in morphospace after phylogenetic transformation (Adams 2014). From the data the net rate of shape evolution
#'  for each group in the multi-dimensional space is calculated, and a ratio of rates is obtained. If three or more groups of species are used, the ratio of 
#'  the maximum to minimum rate is used as a test statistic (see Adams 2014). The function can be used with univariate data (i.e. 
#'  centroid size) if imported as matrix with rownames giving the taxa names.
#'  
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{compare.evol.rates}}.
#'  The generic function, \code{\link{plot}}, produces a histogram of random rate-ratios associated with
#'  the resampling procedure.
#'
#' \subsection{Notes for geomorph 3.0.4 and subsequent versions}{ 
#' Significance testing is now accomplished in one of two ways. First, phylogenetic simulation may be used, in which tips data are 
#' obtained under Brownian motion using a common evolutionary rate pattern for all species on the phylogeny. Specifically, the 
#' common evolutionary rate matrix for all species is used, with the multi-dimensional rate used along the diagonal elements (see 
#' Denton and Adams 2015). This procedure is more general than the original simulation procedure, and retains the desirable 
#' statistical properties of earlier methods, and under a wider array of data types.  Second, significance may be accomplished via 
#' permutation, where data values at the tips are permuted relative to the (see Adams and Collyer 2018). This procedure is shown to 
#' retain all appropriate statistical properties, including rotation-invariance of significance levels (see results of Adams and Collyer 2018).
#' }
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param gp A factor array designating group membership
#' @param method One of "simulation" or "permutation", to choose which approach should be used to assess significance. 
#' @param iter Number of iterations for significance testing
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @keywords analysis
#' @author Dean Adams
#' @seealso  \code{\link[ape]{vcv.phylo}}, \code{\link[geiger]{sim.char}}, \code{\link[Matrix]{nearPD}}
#'  (used in some internal computations)
#' @export
#' @return An object of class "evolrate" returns a list with the following components: 
#'   \item{sigma.d.ratio}{The ratio of maximum to minimum net evolutionary rates.}
#'   \item{P.value}{The significance level of the observed ratio.}
#'   \item{sigma.d.gp}{The phylogenetic net evolutionary rate for each group of species on the phylogeny.}
#'   \item{random.sigma}{The sigma values found in random permutations of the resampling procedure.}
#'   \item{permutations}{The number of random permutations used.}
#'   
#' @references Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
#'  shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.
#' @references Denton, J.S.S., and D.C. Adams. 2015. A new phylogenetic test for comparing 
#' multiple high-dimensional evolutionary rates suggests interplay of evolutionary rates and 
#' modularity in lanternfishes (Myctophiformes; Myctophidae). Evolution. 69:2425-2440.
#' @references Adams, D.C. and M.L. Collyer. 2018. Multivariate comparative methods: evaluations, comparisons, and
#' recommendations. Systematic Biology. 67:14-31.
#' @examples
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#'  gp.end<-factor(c(0,0,1,0,0,1,1,0,0))  #endangered species vs. rest
#'  names(gp.end)<-plethspecies$phy$tip
#' 
#' ER<-compare.evol.rates(A=Y.gpa$coords, phy=plethspecies$phy,method="simulation",gp=gp.end,iter=999)
#' summary(ER)
#' plot(ER)
compare.evol.rates<-function(A,phy,gp,iter=999,method=c("simulation","permutation"),print.progress=TRUE ){
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
  if (!inherits(phy, "phylo")){
    stop("tree must be of class 'phylo.'")}
  ntaxa<-length(phy$tip.label)
  N<-nrow(x)  
  method <- match.arg(method)
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
  if(method == "permutation"){
    ind<-perm.index(N,iter)
    if(ncol(x)==1){
      x.r <-simplify2array(lapply(1:iter, function(i) x[ind[[i]]]))
      dim(x.r) <- c(nrow(x.r), 1, iter); dimnames(x.r)[[1]]<-names(x)
    }else
    x.r <-simplify2array(lapply(1:iter, function(i) x[ind[[i]],]))
  }
  if(method != "permutation") {x.r<-sim.char(phy=phy,par=rate.mat,nsim=iter,model="BM") }
  ones <- matrix(1,N,N); I <- diag(1,N)
  Xadj <- I -crossprod(ones,invC)/sum(invC) 
  Ptrans <- D.mat%*%Xadj
  g<-factor(as.numeric(gp))
  ngps<-nlevels(g)
  gps.combo <- combn(ngps, 2)
  if(nlevels(gp) > 1){
    if(print.progress){
      pb <- txtProgressBar(min = 0, max = iter, initial = 0, style=3) 
      sigma.rand <- sapply(1:iter, function(j) {
        setTxtProgressBar(pb,j)
        fast.sigma.d(as.matrix(x.r[,,j]),Ptrans,g, ngps, gps.combo, N,p )
      })
      close(pb)
    } else sigma.rand <- sapply(1:(iter), 
                                function(j) fast.sigma.d(as.matrix(x.r[,,j]),Ptrans,g, ngps, gps.combo,N,p))
    if(nlevels(gp) == 2) 
      sigma.rand <- random.sigma <- c(sigma.obs$sigma.d.gp.ratio, sigma.rand) else {
        sigma.rand <- cbind(as.vector(sigma.obs$sigma.d.gp.ratio), sigma.rand)
        random.sigma<- sapply(1:(iter+1), function(j) {x <- sigma.rand[,j]; max(x)/min(x)})
      }
    p.val <- pval(random.sigma)
    if(nlevels(gp) > 2) {
      p.val.mat <- dist(sigma.obs$sigma.d.gp)
      p.val.mat[1:length(p.val.mat)] <- apply(sigma.rand, 1, pval)
    } else p.val.mat <- p.val
    if(nlevels(gp)==2) p.val.mat<-p.val
    if(nlevels(gp)>2){
      ratio.vals<-sigma.rand
      tmp.p.val.mat <- sapply(1:ncol(ratio.vals), function(j){ pval(ratio.vals[,j])})
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