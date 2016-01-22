#' Compare modular signal to alternative landmark subsets
#'
#' Function quantifies the degree of modularity between two or more hypothesized modules of Procrustes-aligned 
#'   landmark coordinates and compares this to patterns found by randomly assigning landmarks into subsets
#'
#' The function quantifies the degree of modularity in two or more hypothesized modules of shape data as 
#'   defined by landmark coordinates, and compares this to the degree of modular signal found in random assignment of
#'   landmarks to modules. It is assumed that the landmarks have previously been aligned using Generalized 
#'   Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}]. The degree of modularity  
#'   is quantified using the RV coefficient (Klingenberg 2009). If more than two modules are defined, the average
#'   RV coefficient is utilized (see Klingenberg 2009). The RV coefficient for the observed modular 
#'   hypothesis is then compared to a distribution of values obtained by randomly assigning landmarks into 
#'   subsets, with the restriction that the number of landmarks in each subset is identical to that observed 
#'   in each of the original partitions. A significant modular signal is found when the observed RV coefficient 
#'   is small relative to this distribution (see Klingenberg 2009). A histogram of coefficients obtained via 
#'   resampling is presented, with the observed value designated by an arrow in the plot. 
#'   
#'   Landmark groups can be defined using \code{\link{define.modules}}, or made by hand (see example below).
#'   To use this method with other data (i.e., a set of length measurements), the input A should be a matrix 
#'   of n rows of specimens and variables arranged in columns. 
#'   In this case, the partition.gp input should have each variable assigned to a partition. 
#'   
#' @param A A 3D array (p x k x n) containing GPA-aligned coordinates for all specimens, or a matrix (n x variables)
#' @param partition.gp A list of which landmarks (or variables) belong in which partition (e.g. A,A,A,B,B,B,C,C,C)
#' @param ShowPlot A logical value indicating whether or not the plot should be returned
#' @param iter Number of iterations for significance testing
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{RV}{The estimate of the observed modular signal}
#'   \item{pvalue}{The significance level of the observed signal}
#'   \item{RV.min}{The minimal RV coefficient found via landmark permutation}
#'   \item{RV.min.partitions}{A list of landmarks assigned to partitions that yields the minimal RV coefficient}
#' @references Klingenberg, C. P. 2009. Morphometric integration and modularity in configurations of 
#'   landmarks: tools for evaluating a priori hypotheses. Evol. Develop. 11:405-421.
#' @seealso  \code{\link{define.modules}}
#' @examples
#' data(plethodon) 
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment    
#'  #landmarks on the skull and mandible assigned to partitions
#' land.gps<-c("A","A","A","A","A","B","B","B","B","B","B","B") 
#'
#' compare.modular.partitions(Y.gpa$coords,land.gps,iter=99)
#' #Result implies that the skull and mandible are not independent modules
compare.modular.partitions<-function(A,partition.gp,ShowPlot=TRUE,iter=999){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")}
  partition.gp<-as.factor(partition.gp)
  if (length(dim(A))==3){ x<-two.d.array(A)
                          p<-dim(A)[1]; k<-dim(A)[2]
                          if(length(partition.gp)!=p){stop("Not all landmarks are assigned to a partition.")}
                          gps<-as.factor(rep(partition.gp,k,each = k, length=p*k))  }
  if (length(dim(A))==2){ x<-A
                          if(length(partition.gp)!=ncol(x)){stop("Not all variables are assigned to a partition.")}
                          gps<-partition.gp  }
  ngps<-nlevels(partition.gp)
  S<-cov(x)
  RV.gp<-array(0,dim=c(ngps,ngps))
  for (i in 1:(ngps-1)){
    for (j in 2:ngps){
      S11<-S[which(gps==levels(gps)[i]),which(gps==levels(gps)[i])]
      S22<-S[which(gps==levels(gps)[j]),which(gps==levels(gps)[j])]
      S12<-S[which(gps==levels(gps)[i]),which(gps==levels(gps)[j])]
      S21<-t(S12)
      RV.gp[i,j]<- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22)))
      diag(RV.gp)<-0
    }
  }
  RV.obs<-sum(RV.gp)/(ngps/2*(ngps-1))
  RV.min<-RV.obs; partition.min<-partition.gp
  P.val<-1
  RV.val<-rep(0,iter)
  for(ii in 1:iter){
    partition.gp.r<-sample(partition.gp)
    if (length(dim(A))==3){ gps.r<-as.factor(rep(partition.gp.r,k,each = k, length=p*k)) }  
    if (length(dim(A))==2){ gps.r<-as.factor(partition.gp.r) }
    RV.gp.r<-array(0,dim=c(ngps,ngps))
    for (i in 1:(ngps-1)){
      for (j in 2:ngps){
        S11.r<-S[which(gps.r==levels(gps.r)[i]),which(gps.r==levels(gps.r)[i])]
        S22.r<-S[which(gps.r==levels(gps.r)[j]),which(gps.r==levels(gps.r)[j])]
        S12.r<-S[which(gps.r==levels(gps.r)[i]),which(gps.r==levels(gps.r)[j])]
        S21.r<-t(S12.r)
        RV.gp.r[i,j]<- sum(diag(S12.r%*%S21.r))/sqrt(sum(diag(S11.r%*%S11.r))*sum(diag(S22.r%*%S22.r)))
        diag(RV.gp.r)<-0
      }
    }
    RV.r<-sum(RV.gp.r)/(ngps/2*(ngps-1))
    RV.val[ii]<-RV.r
    if (RV.r< RV.min) {RV.min<-RV.r; partition.min<-partition.gp.r}
    P.val<-ifelse(RV.r<=RV.obs, P.val+1,P.val) 
  }
  RV.val[iter+1]=RV.obs
  P.val<-P.val/(iter+1)
  if(ShowPlot==TRUE){ 
  hist(RV.val,30,freq=TRUE,col="gray",xlab="RV Coefficient")
  arrows(RV.obs,50,RV.obs,5,length=0.1,lwd=2) }
  return(list(RV=RV.obs,pvalue=P.val,RV.min=RV.min,RV.min.partitions=partition.min))
}