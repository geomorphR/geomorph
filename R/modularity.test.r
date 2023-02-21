#' Evaluate the degree of modular signal in shape data
#'
#' Function quantifies the degree of modularity between two or more hypothesized modules of Procrustes shape variables and 
#' compares this to patterns found by randomly assigning landmarks into subsets
#'
#' The function quantifies the degree of modularity in two or more hypothesized modules of Procrustes shape variables, 
#' and compares this to what is expected under the null hypothesis of random assignment
#' of variables to partitions (i.e., neither modular nor integrated structure). 
#' Input may be either a 2D matrix of phenotypic values, or a 3D array of Procrustes shape variables. It 
#' is assumed that the landmarks have previously been aligned using Generalized Procrustes Analysis (GPA); e.g., 
#' with \code{\link{gpagen}}. The degree of modularity is quantified using the CR coefficient (Adams 2016). If more than 
#' two modules are defined, the average pairwise CR coefficient is utilized. The CR coefficient for the observed modular 
#' hypothesis is then compared to a distribution of values obtained by randomly assigning landmarks into subsets, with the 
#' restriction that the number of landmarks in each subset is identical to that observed in each of the original partitions. 
#' A significant modular signal is found when the observed CR coefficient is small relative to this distribution (see Adams 2016). 
#' Such a result implies that there is significantly greater independence among modules than is expected under the null 
#' hypothesis of random associations of variables (neither modular nor integrated structure). This  
#' result is consistent with the identification of significant modular structure in the data. 
#' In addition, a multivariate effect size describing the strength of the effect is 
#'   estimated from the empirically-generated sampling distribution (see details in Adams and Collyer 2019).
#' A histogram of coefficients obtained via 
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
#'  The generic function, \code{\link{plot}}, produces a histogram of random CR values associated with the resampling procedure.
#'  
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for all specimens, or a matrix (n x variables)
#' @param partition.gp A list of which landmarks (or variables) belong in which partition: 
#' (e.g. A, A, A, B, B, B, C, C, C)
#' @param CI A logical argument indicating whether bootstrapping should be used for estimating confidence intervals
#' @param opt.rot A logical argument for whether the optimal rotation for CR should be used for landmark data (default = TRUE)
#' @param iter Number of iterations for significance testing
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#' If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @keywords analysis
#' @export
#' @author Dean Adams
#' @return An object of class "CR" is a list containing the following
#' \item{CR}{Covariance ratio: The estimate of the observed modular signal.}
#' \item{CInterval}{The bootstrapped 95 percent confidence intervals of the CR, if CI = TRUE.}
#'    \item{CR.boot}{The bootstrapped CR values, if CI = TRUE}
#' \item{P.value}{The empirically calculated P-value from the resampling procedure.}
#'   \item{Effect.Size}{The multivariate effect size associated with sigma.d.ratio.}
#'    \item{CR.mat}{For more than two partitions, the pairwise CRs among partitions.}
#'    \item{random.CR}{The CR calculated in each of the random permutations of the resampling procedure.}
#'    \item{permutations}{The number of random permutations used in the resampling procedure.}
#'    \item{call}{The match call.}
#' @references Adams, D.C. 2016.Evaluating modularity in morphometric data: Challenges with the RV coefficient and a 
#' new test measure. Methods in Ecology and Evolution 7:565-572.
#' @references Adams, D.C. and M.L. Collyer. 2019. Comparing the strength of modular signal, and evaluating 
#' alternative modular hypotheses, using covariance ratio effect sizes with morphometric data. 
#' Evolution. 73:2352-2367.
#' @seealso \code{\link{two.b.pls}}, \code{\link{integration.test}}, \code{\link{phylo.modularity}}, and 
#' \code{\link{phylo.integration}}
#' @examples
#' # Not Run
#' # data(pupfish) 
#' # Y.gpa<-gpagen(pupfish$coords, print.progress = FALSE)    #GPA-alignment    
#'  #landmarks on the body and operculum
#' # land.gps<-rep('a',56); land.gps[39:48]<-'b'
#'
#' # MT <- modularity.test(Y.gpa$coords,land.gps,CI=FALSE,iter=99)
#' # summary(MT) # Test summary
#' # plot(MT) # Histogram of CR sampling distribution 
#' # Result implies modularity present

modularity.test<-function(A, partition.gp, iter = 999, CI = FALSE, seed = NULL, 
                          opt.rot = TRUE, print.progress = TRUE){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if (length(dim(A))==3){ x<-two.d.array(A)
           p<-dim(A)[1]; k<-dim(A)[2];n<-dim(A)[3]
           if(length(partition.gp)!=p){stop("Not all landmarks are assigned to a partition.")}
           if(any(table(partition.gp)==1)){stop("Must have at least two landmarks per partition.")}
            }
  if (length(dim(A))==2){ x<-A; k<-1; p <- ncol(A); n <- nrow(A)
           if(length(partition.gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
            }
  gps<-as.factor(partition.gp)
  gps.obs <- as.factor(rep(gps,k,each = k, length=p*k))
  if(any(table(gps.obs)==1)){stop("Must have at least two variables per partition.")}
  if(!is.null(seed) && seed=="random") seed = sample(1:iter, 1)
  ngps<-nlevels(gps)
  if (length(dim(A))==2){
    CR.obs<-CR(x,gps=gps.obs)
    if(ngps > 2) CR.mat <- CR.obs$CR.mat else CR.mat <- NULL
    CR.obs <- CR.obs$CR
    if(print.progress) CR.rand <- apply.CR(x, gps, k=k, iter=iter, seed=seed) else
      CR.rand <- .apply.CR(x, gps, k=k, iter=iter, seed=seed)
    p.val <- 1-pval(CR.rand)  #b/c smaller values more significant 
    if (p.val==0){p.val<-1/(iter+1)}
    Z <- effect.size(CR.rand, center=TRUE) 
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
    if (opt.rot==TRUE){
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
    x <- t(mapply(function(a) matrix(t(a%*%optRot)), Alist))
    } else x <- two.d.array(A) #new
    if(print.progress) {
      cat("\nPerforming permutations\n")
      CR.rand <- apply.CR(x, gps, k=k, iter=iter, seed=seed)
      } else CR.rand <- .apply.CR(x, gps, k=k, iter=iter, seed=seed)
#    CR.rand[1] <- CR.obs <- avgCR
    if(ngps > 2) CR.mat <- CR(x,gps.obs)$CR.mat else CR.mat <- NULL
    p.val <- pval(1/CR.rand)  #b/c smaller values more significant
    Z <- effect.size(CR.rand, center=TRUE) 
    if(CI=="TRUE"){
      CR.boot<- boot.CR(x, gps.obs,k, iter=iter, seed=seed)
      CR.CI<-quantile(CR.boot, c(.025, .975))
    }
    if(CI=="FALSE"){
      CR.boot <- NULL
      CR.CI <- NULL
    }
  }
  
  out <- list(CR=CR.rand[1], CInterval=CR.CI, CR.boot = CR.boot, P.value=p.val, Z = Z,
              CR.mat = CR.mat, random.CR = CR.rand,
              permutations=iter+1, call=match.call())
  class(out) <- "CR"
  out  
}
