#' Evaluate the degree of phylogenetic modular signal in Procrustes shape variables
#'
#' Function quantifies the degree of modularity between two or more hypothesized modules of Procrustes shape variables
#'  in a phylogenetic context and compares this to patterns found by randomly assigning landmarks into subsets
#'
#' The function quantifies the degree of phylogenetic modularity in two or more hypothesized modules of Procrustes shape variables as 
#' defined by landmark coordinates, under a Brownian motion model of evolution. The degree of modularity 
#' is characterized by the covariance ratio covariance ratio (CR: see Adams 2016). The phylogenetic version of the approach 
#' procedure utilizes the evolutionary covariance matrix among traits found under a Brownian
#' motion model of evolution as the basis of the analysis. This is the same matrix used to evaluate patterns of phylogenetic 
#' morphological integration as described in Adams and Felice (2014). 
#' 
#' Input may be either a 2D matrix of phenotypic values, or a 3D array of Procrustes shape variables. It 
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
#' In addition, a multivariate effect size describing the strength of the effect is 
#'   estimated from the empirically-generated sampling distribution (see details in Adams and Collyer 2019).
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param partition.gp A list of which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C)
#' @param CI A logical argument indicating whether bootstrapping should be used for estimating confidence intervals
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Objects of class "CR" from modularity.test return a list of the following:
#'    \item{CR}{Covariance ratio: The estimate of the observed modular signal.}
#'    \item{CInterval}{The bootstrapped 95 percent confidence intervals of the CR, if CI = TRUE.}
#'    \item{CR.boot}{The bootstrapped CR values, if CI = TRUE. 
#'    For more than two partitions, this is the mean CR of pairwise CRs.}
#'    \item{P.value}{The empirically calculated P-value from the resampling procedure.}
#'   \item{Effect.Size}{The multivariate effect size associated with sigma.d.ratio.}
#'    \item{CR.mat}{For more than two partitions, the pairwise CRs among partitions.}
#'    \item{random.CR}{The CR calculated in each of the random permutations of the resampling procedure.}
#'    \item{permutations}{The number of random permutations used in the resampling procedure.}
#'    \item{call}{The match call.}
#'    
#' @references Adams, D.C. 2016.Evaluating modularity in morphometric data: Challenges with the RV coefficient and a 
#' new test measure. Methods in Ecology and Evolution. 7:565-572.
#' @references  Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological 
#' integration and trait covariation in morphometric data using evolutionary covariance 
#' matrices. PLOS ONE. 9(4):e94335.
#' @references Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
#' alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
#' Evolution. 73:2352-2367.
#' @examples
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B") 
#' 
#' MT <- phylo.modularity(Y.gpa$coords, partition.gp=land.gps, 
#' phy=plethspecies$phy, 
#' CI = FALSE, iter=499)
#' summary(MT) # Test summary
#' plot(MT) # Histogram of CR sampling distribution 
phylo.modularity<-function(A, partition.gp, phy, CI = FALSE, 
                           iter = 999, seed = NULL, print.progress = TRUE){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (!inherits(phy, "phylo"))
    stop("phy must be of class 'phylo.'") 
  if (length(dim(A))==3){ 
    p<-dim(A)[1]; k<-dim(A)[2];n<-dim(A)[3]
    gps<-as.factor(partition.gp)
    gps.obs <- as.factor(rep(gps,k,each = k, length=p*k))
    angle <- seq(0,89.95,0.05)
    if(k==2){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                                                              sin(angle[i]*pi/180),-sin(angle[i]*pi/180),cos(angle[i]*pi/180)),ncol=2))      
    }
    if(k==3){
      rot.mat<-lapply(1:(length(angle)), function(i) matrix(c(cos(angle[i]*pi/180),
                                                              sin(angle[i]*pi/180),0,-sin(angle[i]*pi/180),cos(angle[i]*pi/180), 0,0,0,1),ncol=3))      
    }
    Alist <-lapply(1:n,function(j) A[,,j]) # convert array to list
    if(print.progress){
      cat("\nFinding the optimal rotation for CR\n")
      pb <- txtProgressBar(min = 0, max = length(rot.mat), initial = 0, style=3) 
      rotatedCRs <- sapply(1:length(rot.mat), function(j) {
        r <- rot.mat[[j]]
        rotA <- t(mapply(function(a) matrix(t(a%*%r)), Alist))
        setTxtProgressBar(pb,j)
        quick.CR(rotA, gps=gps.obs)
      })
      close(pb)
    } else {
      rotatedCRs <-sapply(1:length(rot.mat), function(j) {
        r <- rot.mat[[j]]
        rotA <- t(mapply(function(a) matrix(t(a%*%r)), Alist))
        quick.CR(rotA, gps=gps.obs)
      })
    }
    avgCR <- mean(rotatedCRs)
    angCheck <- abs(rotatedCRs-avgCR)
    optAngle <- angle[angCheck==min(angCheck)]; optAngle<-optAngle[1]
    # Optimal rotation 
    if(k==2) optRot <- matrix(c(cos(optAngle*pi/180),
                                sin(optAngle*pi/180),-sin(optAngle*pi/180),cos(optAngle*pi/180)),ncol=2) else
                                  optRot <- matrix(c(cos(optAngle*pi/180),
                                                     sin(optAngle*pi/180),0,-sin(optAngle*pi/180),cos(optAngle*pi/180), 0,0,0,1),ncol=3)
    x <- t(mapply(function(a) matrix(t(a%*%optRot)), Alist)); rownames(x) <- dimnames(A)[[3]]
    p<-dim(A)[1]; k<-dim(A)[2];n<-dim(A)[3]
    if(length(partition.gp)!=p){stop("Not all landmarks are assigned to a partition.")}
  }
  if (length(dim(A))==2){ x<-A; k <-1; p<-ncol(A)
  if(length(partition.gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
  }
  phy.parts<-phylo.mat(x,phy)
  invC<-phy.parts$invC; D.mat<-phy.parts$D.mat
  one<-matrix(1,nrow(x)); I = diag(1,nrow(x),) 
  Ptrans<-D.mat%*%(I-one%*%crossprod(one,invC)/sum(invC))
  x <- Ptrans%*%x
  if (length(dim(A))==3){x <- arrayspecs(x,p=p,k=k)}
  res <- modularity.test(x,partition.gp,CI=CI, iter=iter, seed=seed, opt.rot = FALSE,
                         print.progress=print.progress)
  out <- list(CR=res$CR, CInterval=res$CInterval, CR.boot = res$CR.boot, 
              P.value=res$P.value, Z = res$Z,
              CR.mat = res$CR.mat, random.CR = res$random.CR,
              permutations=iter+1, call=match.call())
  class(out) <- "CR"
  out  
}
