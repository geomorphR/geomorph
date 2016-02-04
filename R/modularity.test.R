#' Evaluate the degree of modular signal in morphometric datasets
#'
#' Function quantifies the degree of modularity between two or more hypothesized modules of Procrustes-aligned 
#' landmark coordinates and compares this to patterns found by randomly assigning landmarks into subsets
#'
#' The function quantifies the degree of modularity in two or more hypothesized modules of shape data as 
#' defined by landmark coordinates, and compares this to what is expected under the null hypothesis of random assignment
#' of variables to partitions (i.e., neither modular nor integrated structure). 
#' Input may be either a 2D matrix of phenotypic values, or a 3D array of aligned Procrustes coordinates. It 
#' is assumed that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA); e.g., 
#' with \code{\link{gpagen}}. The degree of modularity is quantified using the CR coefficient (Adams 2016). If more than 
#' two modules are defined, the average pairwise CR coefficient is utilized. The CR coefficient for the observed modular 
#' hypothesis is then compared to a distribution of values obtained by randomly assigning landmarks into subsets, with the 
#' restriction that the number of landmarks in each subset is identical to that observed in each of the original partitions. 
#' A significant modular signal is found when the observed CR coefficient is small relative to this distribution (see Adams 2016). 
#' Such a result implies that there is significantly greater independence among modules than is expected under the null 
#' hypothesis of random associations of variables (neither modular nor integrated structure). This  
#' result is consistent with the identification of significant modular structure in the data. A histogram of coefficients obtained via 
#' resampling is presented, with the observed value designated by an arrow in the plot. For landmark data, the CR coefficient 
#' found from the average CR across a 90 degree rotation of the data is used as the test statistic (see Adams 2016). For all
#' data, the CR coefficient is returned, and (optionally)  its 95% confidence intervals (based on bootstrapping) may be requested. NOTE that 
#' for landmark data, estimation of the CI can take some time to compute.
#' 
#' Landmark groups can be defined using \code{\link{define.modules}}, or made by hand (see example below).
#' To use this method with other data (i.e., a set of length measurements), the input A should be a matrix 
#' of n rows of specimens and variables arranged in columns. 
#' In this case, the partition.gp input should have each variable assigned to a partition. 
#' 
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{modularity.test}}.
#'  The generic function, \code{\link{plot}}, produces a two-block.pls plot.  This function calls \code{\link{plot.pls}}, which has two additional
#'  arguments (with defaults): label = NULL, warpgrids = TRUE.  These arguments allow one to include a vector to label points and a logical statement to
#'  include warpgrids, respectively.  Warpgrids can only be included for 3D arrays of Procrustes residuals. The plot is a plot of PLS scores from 
#'  Block1 versus Block2 performed for the first set of PLS axes. 
#'  
#' @param A A 3D array (p x k x n) containing GPA-aligned coordinates for all specimens, or a matrix (n x variables)
#' @param partition.gp A list of which landmarks (or variables) belong in which partition (e.g. A,A,A,B,B,B,C,C,C)
#' @param CI A logical argument indicating whether bootstrapping should be used for estimating confidence intervals
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @keywords analysis
#' @export
#' @author Dean Adams
#' @return An object of class "CR" is a list containing the following
#' \item{CR}{Covariance ratio: The estimate of the observed modular signal.}
#' \item{CInterval}{The bootstrapped 95 percent confidence intervals of the CR, if CI = TRUE.}
#'    \item{CR.boot}{The bootstrapped CR values, if CI = TRUE}
#' \item{P.value}{The empirically calculated P-value from the resampling procedure.}
#'    \item{CR.mat}{For more than two partitions, the pairwise CRs among partitions.}
#'    \item{random.CR}{The CR calculated in each of the random permutations of the resampling procedure.}
#'    \item{permutations}{The number of random permutations used in the resampling procdure.}
#'    \item{call}{The match call.}
#' @references Adams, D.C. 2016.Evaluating modularity in morphometric data: Challenges with the RV coefficient and a 
#' new test measure. Methods in Ecology and Evolution. (Accepted). 
#' @seealso \code{\link{two.b.pls}}, \code{\link{integration.test}}, \code{\link{phylo.modularity}}, and 
#' \code{\link{phylo.integration}}
#' @examples
#' data(pupfish) 
#' Y.gpa<-gpagen(pupfish$coords)    #GPA-alignment    
#'  #landmarks on the body and operculum
#' land.gps<-rep('a',56); land.gps[39:48]<-'b'
#'
#' MT <- modularity.test(Y.gpa$coords,land.gps,CI=FALSE,iter=499)
#' summary(MT) # Test summary
#' plot(MT) # Histogram of CR sampling distribution 
#' # Result implies modularity present

modularity.test<-function(A,partition.gp,iter=999, CI=FALSE,seed=NULL){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A))==3){ x<-two.d.array(A)
           p<-dim(A)[1]; k<-dim(A)[2];n<-dim(A)[3]
           if(length(partition.gp)!=p){stop("Not all landmarks are assigned to a partition.")}
            }
  if (length(dim(A))==2){ x<-A; k<-1; p <- ncol(A); n <- nrow(A)
           if(length(partition.gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
            }
  gps<-factor(as.numeric(as.factor(partition.gp)))
  gps.obs <- as.factor(rep(gps,k,each = k, length=p*k))
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  ngps<-nlevels(gps)
  if (length(dim(A))==2){
    CR.obs<-CR(x,gps=gps.obs)
    if(ngps > 2) CR.mat <- CR.obs$CR.mat else CR.mat <- NULL
    CR.obs <- CR.obs$CR
    CR.rand <- apply.CR(x, gps, k=k, iter=iter, seed=seed)
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
    rotatedCRs <-sapply(1:length(rot.mat), function(j) {
      r <- rot.mat[[j]]
      rotA <- t(mapply(function(a) matrix(t(a%*%r)), Alist))
      CR(rotA, gps=gps.obs)$CR
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
    CR.rand <- apply.CR(x, g=gps, k=k, iter=iter, seed=seed)
    CR.rand[1] <- CR.obs <- avgCR
    if(ngps > 2) CR.mat <- CR(x,gps.obs)$CR.mat else CR.mat <- NULL
    p.val <- pval(1/CR.rand)  #b/c smaller values more significant
    if(CI=="TRUE"){
      CR.boot<- boot.CR(x, gps.obs,k, iter=iter, seed=seed)
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
