#' Evaluate the degree of phylogenetic modular signal in morphometric datasets
#'
#' Function quantifies the degree of modularity between two or more hypothesized modules of Procrustes-aligned 
#'   landmark coordinates in a phylogenetic context and compares this to patterns found by randomly assigning landmarks into subsets
#'
#' The function quantifies the degree of phylogenetic modularity in two or more hypothesized modules of shape data as 
#' defined by landmark coordinates, under a Brownian motion model of evolution. The degree of modularity 
#' is characterized by the covariance ratio covariance ratio (CR: see Adams 2016). The phylogenetic version of the approach 
#' procedure utilizes the evolutionary covariance matrix among traits found under a Brownian
#' motion model of evolution as the basis of the analysis. This is the same matrix used to evaluate patterns of phylogenetic 
#' morphological integration as described in Adams and Felice (2014). 
#' 
#' Input may be either a 2D matrix of phenotypic values, or a 3D array of aligned Procrustes coordinates. It 
#' is assumed that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA) [e.g., 
#' with \code{\link{gpagen}}]. The degree of modularity is quantified using the CR coefficient (Adams 2016). If more than 
#' two modules are defined, the average pairwise CR coefficient is utilized. The CR coefficient for the observed modular 
#' hypothesis is then compared to a distribution of values obtained by randomly assigning landmarks into subsets, with the 
#' restriction that the number of landmarks in each subset is identical to that observed in each of the original partitions. 
#' A significant modular signal is found when the observed CR coefficient is small relative to this distribution (see Adams 2016). 
#' Such a result implies that there is significantly greater independence among modules than is expected under the null 
#' hypothesis of random associations of variables (neither modular nor integrated structure). This  
#' result is consistent with the identification of significant modular structure in the data. For landmark data, the CR coefficient 
#' found from the average CR across a 90 degree rotation of the data is used as the test statistic (see Adams 2016). 
#'
#' @param A A 3D array (p x k x n) containing GPA-aligned coordinates for all specimens, or a matrix (n x variables)
#' @param partition.gp A list of which landmarks (or variables) belong in which partition (e.g. A,A,A,B,B,B,C,C,C)
#' @param CI A logical argument indicating whether bootstrapping should be used for estimating confidence intervals
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Objects of class "CR" from modularity.test return a list of the following:
#'    \item{CR}{Covariance ratio: The estimate of the observed modular signal.}
#'    \item{CInterval}{The bootstrapped 95 percent confidence intervals of the CR, if CI = TRUE.}
#'    \item{CR.boot}{The bootstrapped CR values, if CI = TRUE
#'    (For more than two partitions, this is the mean CR of pairwise CRs.)}
#'    \item{P.value}{The empirically calculated P-value from the resampling procedure.}
#'    \item{CR.mat}{For more than two partitions, the pairwise CRs among partitions.}
#'    \item{random.CR}{The CR calculated in each of the random permutations of the resampling procedure.}
#'    \item{permutations}{The number of random permutations used in the resampling procedure.}
#'    \item{call}{The match call.}
#'    
#' @references Adams, D.C. 2016.Evaluating modularity in morphometric data: Challenges with the RV coefficient and a 
#' new test measure. Methods in Ecology and Evolution. (Accepted). 
#' @references  Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological 
#' integration and trait covariation in morphometric data using evolutionary covariance 
#' matrices. PLOS ONE. 9(4):e94335.
#' @examples
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
#' MT <- phylo.modularity(Y.gpa$coords, partition.gp=land.gps, phy=plethspecies$phy, 
#' CI = FALSE, iter=999)
#' summary(MT) # Test summary
#' plot(MT) # Histogram of CR sampling distribution 
phylo.modularity<-function(A,partition.gp,phy, CI=FALSE, iter=999, seed=NULL){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (class(phy) != "phylo")
    stop("phy must be of class 'phylo.'") 
  if (length(dim(A))==3){ x<-two.d.array(A)
           p<-dim(A)[1]; k<-dim(A)[2];n<-dim(A)[3]
           if(length(partition.gp)!=p){stop("Not all landmarks are assigned to a partition.")}
           }
  if (length(dim(A))==2){ x<-A; k <-1; p <- ncol(A)
           if(length(partition.gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
           }
  gps<-factor(as.numeric(as.factor(partition.gp)))
  gps.obs <- as.factor(rep(gps,k,each = k, length=p*k))
  ngps<-nlevels(gps)
  Nspec<-num.taxa.X<-nrow(x)
  namesX<-rownames(x)
  if (is.null(namesX)){
    stop("No specimen names in data matrix. Please assign specimen names.")  } 
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && length(phy$tip.label) < num.taxa.X)
    stop("Tree is missing some taxa present in the data matrix") 
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && num.taxa.X < length(phy$tip.label)) 
    stop("Tree contains some taxa not present in present in the data matrix")  
  phy.parts<-phylo.mat(x,phy)
  invC<-phy.parts$invC; D.mat<-phy.parts$D.mat
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  if (length(dim(A))==2){
    CR.obs<-CR.phylo(x,invC,gps.obs)
    if(ngps > 2) CR.mat <- CR.obs$CR.mat else CR.mat <- NULL
    CR.obs <- CR.obs$CR
    CR.rand <- apply.phylo.CR(x,invC, gps, k, iter=iter, seed=seed)
    p.val <- 1-pval(CR.rand)  #b/c smaller values more significant
    if (p.val==0){p.val<-1/(iter+1)}
    if(CI=="TRUE"){
      CR.boot<- boot.CR(x, gps, k,iter=iter, seed=seed)
      CR.CI<-quantile(CR.boot, c(.025, .975)) 
    }
    if(CI=="FALSE"){
      CR.boot <- NULL
      CR.CI <- NULL
    }
  }
  if (length(dim(A))==3){
    angle <- seq(0,89.95,.05)
    if(k==2){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
              sin(angle[i]*pi/180),-sin(angle[i]*pi/180),cos(angle[i]*pi/180)),ncol=2))      
    }
    if(k==3){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                sin(angle[i]*pi/180),0,-sin(angle[i]*pi/180),cos(angle[i]*pi/180), 0,0,0,1),ncol=3))      
    }
    Alist <-lapply(1:n,function(j) A[,,j]) # convert array to list
    rotatedCRs <-sapply(1:length(rot.mat), function(j) {
      r <- rot.mat[[j]]
      rotA <- t(mapply(function(a) matrix(t(a%*%r)), Alist))
      CR.phylo(rotA,invC,gps)$CR
    })
    avgCR <- mean(rotatedCRs)
    angCheck <- abs(rotatedCRs-avgCR)
    optAngle <- angle[angCheck==min(angCheck)]
    # Optimal rotation 
    if(k==2) optRot <- matrix(c(cos(optAngle*pi/180),
             sin(optAngle*pi/180),-sin(optAngle*pi/180),cos(optAngle*pi/180)),ncol=2) else
              optRot <- matrix(c(cos(optAngle*pi/180),
               sin(optAngle*pi/180),0,-sin(optAngle*pi/180),cos(optAngle*pi/180), 0,0,0,1),ncol=3)
    x <- t(mapply(function(a) matrix(t(a%*%optRot)), Alist))
    CR.rand <- apply.phylo.CR(x, invC, gps, k, iter=iter, seed=seed)
    CR.rand[1] <- CR.obs <- avgCR
    if(ngps > 2) CR.mat <- CR(x,gps)$CR.mat else CR.mat <- NULL
    p.val <- pval(1/CR.rand)  #b/c smaller values more significant
    if(CI=="TRUE"){
      CR.boot<- boot.CR(x, gps, k,iter=iter, seed=seed)
      CR.CI<-quantile(CR.boot, c(.025, .975)) 
    }
    if(CI=="FALSE"){
      CR.boot <- NULL
      CR.CI <- NULL
    }
  }
  
  out <- list(CR=CR.obs, CInterval=CR.CI, CR.boot = CR.boot, P.value=p.val,
              CR.mat = CR.mat, random.CR = CR.rand,
              permutations=iter+1, call=match.call())
  class(out) <- "CR"
  out  
}
