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
#' also returned. 
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
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with 
#'  \code{\link{compare.multi.evol.rates}}.
#'  The generic function, \code{\link{plot}}, produces a histogram of random rate-ratios associated with
#'  the resampling procedure.
#'
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param gp A factor array designating group membership
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param Subset A logical value indicating whether or not the traits are subsets from a single 
#' landmark configuration (default is TRUE)
#' @param iter Number of iterations for significance testing
#' @keywords analysis
#' @author Dean Adams
#' @export
#' @return An object of class "evolrate" returns a list of the following: 
#'   \item{rates.all}{The phylogenetic evolutionary rates for each trait.}
#'   \item{rate.ratio}{The ratio of maximum to minimum evolutionary rates.}
#'   \item{pvalue}{The significance level of the observed rate ratio.}
#'   \item{pvalue.gps}{Matrix of pairwise significance levels comparing each pair of rates.}
#'   \item{call}{The matched call.}
#'   
#' @references Adams, D.C. 2014. Quantifying and comparing phylogenetic evolutionary rates for 
#'  shape and other high-dimensional phenotypic data. Syst. Biol. 63:166-177.
#' @references Denton, J.S.S., and D.C. Adams. 2015. A new phylogenetic test for comparing 
#' multiple high-dimensional evolutionary rates suggests interplay of evolutionary rates and 
#' modularity in lanternfishes (Myctophiformes; Myctophidae). Evolution. 69:2425-2440.
#' @examples
#' 
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment    
#' land.gp<-c("A","A","A","A","A","B","B","B","B","B","B")  #mandible and cranium subsets
#'
#' EMR<-compare.multi.evol.rates(A=Y.gpa$coords,gp=land.gp, 
#'     Subset=TRUE, phy= plethspecies$phy,iter=999)
#' summary(EMR)
#' plot(EMR)
compare.multi.evol.rates<-function(A,gp,phy,Subset=TRUE,iter=999){
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
  ngps<-nlevels(gps)
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
  x<-as.matrix(x)
  phy.parts<-phylo.mat(x,phy)
  invC<-phy.parts$invC; D.mat<-phy.parts$D.mat;C = phy.parts$C
  sigma.obs<-sigma.d.multi(x,invC,D.mat,gps,Subset)
  R<-sigma.obs$R; diag(R)<-sigma.obs$rate.global
  R<-matrix(nearPD(R,corr=FALSE)$mat,nrow=ncol(R),ncol=ncol(R))
  x.sim<-sim.char(phy,R,nsim=iter) 
  sigma.rand <- sapply(1:(iter), function(j) sigma.d.multi(as.matrix(x.sim[,,j]),invC,D.mat,gps,Subset))
  random.sigma<-c(sigma.obs$sigma.d.ratio,as.vector(unlist(sigma.rand[1,])))
  p.val <- pval(random.sigma)
  ratio.vals<-matrix(NA,nrow=(iter+1),ncol=length(unlist(sigma.obs[4])))
  ratio.vals[1,]<-as.vector(sigma.obs$sigma.d.gp.ratio)
  for(i in 1:iter) ratio.vals[i+1,]<-as.vector(unlist(sigma.rand[4,][[i]]))
  tmp.p.val.mat <- sapply(1:ncol(ratio.vals), function(j){ pval(ratio.vals[,j])})
  p.val.mat<-dist(matrix(0,length(tmp.p.val.mat)))
  if(ngps==2) p.val.mat<-tmp.p.val.mat
  if(ngps>2){
    for(i in 1:length(p.val.mat)) p.val.mat[[i]] <- tmp.p.val.mat[i]
  }
  out <- list(sigma.d.ratio = sigma.obs$sigma.d.ratio, P.value=p.val,
              sigma.d.gp = sigma.obs$rate.gps,
              sigma.d.gp.ratio = sigma.obs$sigma.d.gp.ratio,
              pairwise.pvalue = p.val.mat, groups = levels(gp), 
              random.sigma = random.sigma, permutations=iter+1, call= match.call())
  
  class(out) <- "evolrate"
  out 

}
