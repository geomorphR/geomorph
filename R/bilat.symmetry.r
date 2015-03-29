#' Analysis of bilateral symmetry
#'
#' Function performs an analysis of directional and fluctuating asymmetry for bilaterally symmetric objects 
#'
#' The function quantifies components of shape variation for a set of specimens as described by their patterns of symmetry
#'  and asymmetry. Here, shape variation is decomposed into variation among individuals, variation among sides (directional 
#'  asymmetry), and variation due to an individual x side interaction (fluctuating symmetry). These components are then 
#'  statistically evaluated using Procrustes ANOVA and Goodall's F tests (i.e. an isotropic model of shape variation). Methods for both 
#'  matching symmetry and object symmetry can be implemented. Matching symmetry is when each object contains mirrored 
#'  pairs of structures (e.g., right and left hands) while object symmetry is when a single object is symmetric 
#'  about a midline (e.g., right and left sides of human faces). Analytical and computational details concerning the 
#'  analysis of symmetry in geometric morphometrics can be found in Mardia et al. 2000; Klingenberg et al. 2002.
#'
#' Analyses of symmetry for matched pairs of objects is implemented when {object.sym=FALSE}. Here, a 3D array [p x k x 2n] 
#'  contains the landmark coordinates for all pairs of structures (2 structures for each of n specimens). Because the two sets of 
#'  structures are on opposite sides, they represent mirror images, and one set must be reflected prior to the analysis to 
#'  allow landmark correspondence. IT IS ASSUMED THAT THE USER HAS DONE THIS PRIOR TO PERFORMING THE SYMMETRY ANALYSIS. 
#'  Reflecting a set of specimens may be accomplished by multiplying one coordinate dimension 
#'  by '-1' for these structures (either the x-, the y-, or the z-dimension). A vector containing information on individuals 
#'  and sides must also be supplied. Replicates of each specimen may also be included in the dataset, and when specified will be 
#'  used as measurement error (see Klingenberg and McIntyre 1998). 
#' 
#' Analyses of object symmetry is implemented when {object.sym=TRUE}. Here, a 3D array [p x k x n] contains the landmark 
#'  coordinates for all n specimens. To obtain information about asymmetry, the function generates a second set of objects 
#'  by reflecting them about one of their coordinate axes. The landmarks across the line of symmetry are then relabeled to obtain
#'  landmark correspondence. The user must supply a list of landmark pairs. A vector containing information on individuals 
#'  must also be supplied. Replicates of each specimen may also be included in the dataset, and when specified will be 
#'  used as measurement error. 
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens [for "object.sym=FALSE, A is of dimension (n x k x 2n)]
#' @param ind A vector containing labels for each individual. For matching symmetry, the matched pairs receive the same 
#' label (replicates also receive the same label).
#' @param side An optional vector (for matching symmetry) designating which object belongs to which 'side-group'
#' @param replicate An optional vector designating which objects belong to which group of replicates
#' @param object.sym A logical value specifying whether the analysis should proceed based on object symmetry {=TRUE} or matching symmetry {=FALSE}
#' @param land.pairs An optional matrix (for object symmetry) containing numbers for matched pairs of landmarks across the line of symmetry 
#' @param warpgrids A logical value indicating whether deformation grids for directional and fluctuating components
#' of asymmetry
#' @param mesh A mesh3d object to be warped to represent shape deformation of the directional and fluctuating components
#' of asymmetry if {warpgrids= TRUE} (see \code{\link{warpRefMesh}}).
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @keywords analysis
#' @export
#' @author Dean Adams & Emma Sherratt
#' @return Function returns a list with the following components:
#'   \item{ANOVA.shape}{Procrustes ANOVA table assessing patterns of shape asymmetry}
#'   \item{ANOVA.size}{Procrustes ANOVA table assessing patterns of shape asymmetry (when {object.sym=FALSE})} 
#'   \item{symm.shape}{The symmetric component of shape variation of the aligned specimens ({when verbose=TRUE})}
#'   \item{asymm.shape}{The asymmetric component of shape variation of the aligned specimens ({when verbose=TRUE})}
#' 
#' @references Klingenberg, C.P. and G.S. McIntyre. 1998. Quantitative genetics of geometric shape in the mouse mandible. Evolution. 55:2342-2352.
#' @references Mardia, K.V., F.L. Bookstein, and I.J. Moreton. 2000. Statistical assessment of bilateral symmetry of shapes. Biometrika. 87:285-300.
#' @references Klingenberg, C.P., M. Barluenga, and A. Meyer. 2002. Shape analysis of symmetric structures: quantifying variation among
#' individuals and asymmetry. Evolution. 56:1909-1920.
#' @examples
#' #Example of matching symmetry
#'
#' data(mosquito)
#' bilat.symmetry(mosquito$wingshape,ind=mosquito$ind,side=mosquito$side,
#' replicate=mosquito$replicate,object.sym=FALSE)
#'
#' #Example of object symmetry
#'
#' data(scallops)
#' bilat.symmetry(scallops$coorddata,ind=scallops$ind,object.sym=TRUE,land.pairs=scallops$land.pairs)
bilat.symmetry<-function(A,ind=NULL,side=NULL,replicate=NULL,object.sym=FALSE,land.pairs=NULL,
      warpgrids = TRUE, mesh=NULL, verbose =FALSE){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(is.null(ind)){stop("Individuals not specified.")}
  ind<-as.factor(ind)
  n<-dim(A)[3];   k<-dim(A)[2];  p<-dim(A)[1]; shpsp<-k*p-k-k*(k-1)/2-1; nind<-nlevels(ind) 
    spec.names<-dimnames(A)[[3]]
  if(!is.null(replicate)){replicate<-as.factor(replicate); nrep<-nlevels(replicate) }
  if(object.sym==FALSE){
    if(is.null(side)){stop("Sides not specified.")}  
    side<-as.factor(side)
    gpa.res<-gpagen(A,ShowPlot=FALSE)
    shape<-two.d.array(gpa.res$coords)
    f1<-"shape~ind*side"; if(!is.null(replicate)){f1<-paste(f1,"ind:side:replicate",sep="+")}
    f1<-as.formula(f1)
    f2<-"gpa.res$Csize~ind*side"; if(!is.null(replicate)){f2<-paste(f2,"ind:side:replicate",sep="+")}
    f2<-as.formula(f2)
    res.shape<-anova.parts(f1,Yalt="observed")$table    
    res.shape<-res.shape[1:(dim(res.shape)[1]-2),1:(dim(res.shape)[2]-2)]
    res.shape[1,1]<-(nind-1)*shpsp;res.shape[2,1]<-shpsp; res.shape[3,1]<-(nind-1)*shpsp
    if(!is.null(replicate)){res.shape[4,1]<-(nrep-1)*nind*2*shpsp}
    res.shape[,3]<-res.shape[,2]/res.shape[,1] 
    Goodall.F<-array(NA,nrow(res.shape));Goodall.F[1]<-res.shape[1,3]/res.shape[3,3]
    Goodall.F[2]<-res.shape[2,3]/res.shape[3,3]; if(!is.null(replicate)){Goodall.F[3]<-res.shape[3,3]/res.shape[4,3]}
    P.param<-array(NA,nrow(res.shape)); P.param[1]<-1-pf(Goodall.F[1],res.shape[1,1],res.shape[3,1])
    P.param[2]<-1-pf(Goodall.F[2],res.shape[2,1],res.shape[3,1]); if(!is.null(replicate)){
      P.param[3]<-1-pf(Goodall.F[3],res.shape[3,1],res.shape[4,1])}
    P.param[zapsmall(P.param)==0]=0.00001
    res.shape<-cbind(res.shape,Goodall.F,P.param)
    res.size<-summary(aov(f2))[[1]]; res.size<-res.size[,(1:3)] 
    colnames(res.size) <- c("df", "SS", "MS")
    F<-array(NA,nrow(res.size)); F[1]<-res.size[1,3]/res.size[3,3]
    F[2]<-res.size[2,3]/res.size[3,3]; if(!is.null(replicate)){F[3]<-res.size[3,3]/res.size[4,3]}
    P<-array(NA,nrow(res.size)); P[1]<-1-pf(F[1],res.size[1,1],res.size[3,1])
    P[2]<-1-pf(F[2],res.size[2,1],res.size[3,1]); if(!is.null(replicate)){
      P[3]<-1-pf(F[3],res.size[3,1],res.size[4,1])}
    P[zapsmall(P)==0]=0.00001
    res.size<-cbind(res.size,F,P);  if(!is.null(replicate)){res.size<-res.size[(1:4),];
                                                            rownames(res.size)[4]<-"replicate"}
    class(res.shape) <- c("anova", class(res.shape))
    class(res.size) <- c("anova", class(res.size))
    DA.mns <- arrayspecs((rowsum(predict(lm(shape~side)), 
                                 side)/as.vector(table(side))),p,k)
    MSCP.FA<-summary(manova(lm(f1)))$SS[[3]]/res.shape[3,1]
    eig.FA<-eigen(MSCP.FA); PC1.eigval<-eig.FA$values[1]/sum(eig.FA$values)
    PC1<-shape%*%eig.FA$vec[,1]
    FA.mns<-arrayspecs((rbind(shape[which.min(PC1),],shape[which.max(PC1),])),p,k) 
    mn.shape<-mshape(gpa.res$coords)
    symm.component<-arrayspecs(predict(lm(shape~ind)),p,k)
    asymm.component<-array(data=NA,dim=c(p,k,n))     
      dimnames(symm.component)[[3]] <- dimnames(asymm.component)[[3]]<-spec.names
    for (i in 1:n){ asymm.component[,,i]<-(gpa.res$coords[,,i]-symm.component[,,i]) + mn.shape}
      if(warpgrids == TRUE){
        if(k==2){  
          par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
          plotAllSpecimens(symm.component)
          plotAllSpecimens(asymm.component)
          plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method="TPS",main="Directional Asymmetry")
          plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method="TPS",main="Fluctuating Asymmetry")
          mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
          mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
          par(mfrow=c(1,1))
        }
        if (k==3){
          if (is.null(mesh)==TRUE){
            open3d()
            plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method="points",main="Directional Asymmetry")
            open3d()
            plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method="points",main="Fluctuating Asymmetry")
            } 
          if(is.null(mesh)==FALSE){
            plotRefToTarget(DA.mns[,,1],DA.mns[,,2],mesh,method="surface")
            title3d(main="Directional Asymmetry")
            plotRefToTarget(FA.mns[,,1],FA.mns[,,2],mesh,method="surface")
            title3d(main="Fluctuating Asymmetry")
            }
        }
      layout(1) 
      }
    if(verbose==TRUE){
    	class(res.size) <- c("anova", class(res.size))
    	class(res.shape) <- c("anova", class(res.shape))
    	return(list(symm.shape=symm.component,asymm.shape=asymm.component, 
                                     ANOVA.size=res.size,ANOVA.Shape=res.shape)) }
    if(verbose==FALSE){
    	class(res.size) <- c("anova", class(res.size))
    	class(res.shape) <- c("anova", class(res.shape))
    	return(list(ANOVA.size=res.size,ANOVA.Shape=res.shape)) }
    }
  
  if(object.sym==TRUE){
    if(is.null(land.pairs)){stop("Landmark pairs not specified.")} 
    npairs<-nrow(land.pairs); nl<-p-2*npairs
    A2<-A; 
    for (i in 1:n){
      for (j in 1:nrow(land.pairs)){
        A2[land.pairs[j,1],,i]<-A[land.pairs[j,2],,i]
        A2[land.pairs[j,2],,i]<-A[land.pairs[j,1],,i]            
      }
    }
    A2[,1,]<-A2[,1,]*-1
    A<-array(c(A,A2), c(p,k, 2*n))
    ind<-rep(ind,2);side<-gl(2,n); if(!is.null(replicate)){replicate<-rep(replicate,2)}
    gpa.res<-gpagen(A,ShowPlot = FALSE)
    shape<-two.d.array(gpa.res$coords)    
    f1<-"shape~ind*side"; if(!is.null(replicate)){f1<-paste(f1,"ind:side:replicate",sep="+")}
    f1<-as.formula(f1)
    res.shape<-anova.parts(f1,Yalt="observed")$table    
    res.shape<-res.shape[1:(dim(res.shape)[1]-2),1:(dim(res.shape)[2]-2)]
    res.shape[,2]<-res.shape[,2]/2 
    res.shape[,2]<-res.shape[,2]/2 
    res.shape[2,1]<-ifelse(k==2,((2*npairs+nl-2)),((3*npairs+nl-3)))
    res.shape[1,1]<-res.shape[3,1]<-(nind-1)*res.shape[2,1]
    if(k==3){res.shape[1,1]=res.shape[1,1]+((nind-1)*(nl-1))}
    if(!is.null(replicate)){res.shape[4,1]<-(nrep-1)*nind*shpsp}
    res.shape[,3]<-res.shape[,2]/res.shape[,1] 
    Goodall.F<-array(NA,nrow(res.shape));Goodall.F[1]<-res.shape[1,3]/res.shape[3,3]
    Goodall.F[2]<-res.shape[2,3]/res.shape[3,3]; if(!is.null(replicate)){Goodall.F[3]<-res.shape[3,3]/res.shape[4,3]}
    P.param<-array(NA,nrow(res.shape)); P.param[1]<-1-pf(Goodall.F[1],res.shape[1,1],res.shape[3,1])
    P.param[2]<-1-pf(Goodall.F[2],res.shape[2,1],res.shape[3,1]); if(!is.null(replicate)){
      P.param[3]<-1-pf(Goodall.F[3],res.shape[3,1],res.shape[4,1])}
    P.param[zapsmall(P.param)==0]=0.00001
    res.shape<-cbind(res.shape,Goodall.F,P.param)
    class(res.shape) <- c("anova", class(res.shape))
    
    DA.mns <- arrayspecs((rowsum(predict(lm(shape~side)), 
                                 side)/as.vector(table(side))),p,k)
    MSCP.FA<-summary(manova(lm(f1)))$SS[[3]]/res.shape[3,1]
    eig.FA<-eigen(MSCP.FA); PC1.eigval<-eig.FA$values[1]/sum(eig.FA$values)
    PC1<-shape%*%eig.FA$vec[,1]
    FA.mns<-arrayspecs((rbind(shape[which.min(PC1),],shape[which.max(PC1),])),p,k) 
    symm.component<-arrayspecs(predict(lm(shape~ind)),p,k)
    symm.component<-symm.component[,,(1:n)]
    asymm.component<-array(data=NA,dim=c(p,k,n)) 
      dimnames(symm.component)[[3]] <- dimnames(asymm.component)[[3]]<-spec.names
    mn.shape<-mshape(gpa.res$coords)
    for (i in 1:n){ asymm.component[,,i]<-(gpa.res$coords[,,i]-symm.component[,,i]) + mn.shape}
    if(warpgrids==TRUE){
        if(k==2){  
          par(mfrow=c(2,2),oma=c(1.5,0,1.5,0))
          plotAllSpecimens(symm.component)
          plotAllSpecimens(asymm.component)
          plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method="TPS",main="Directional Asymmetry")
          plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method="TPS",main="Fluctuating Asymmetry")
          mtext("Symmetric Shape Component (left) and Asymmetric Shape Component (right)",outer = TRUE,side=3)
          mtext("Mean directional (left) and fluctuating (right) asymmetry",side = 1, outer = TRUE)
        }
        if (k==3){
          if(is.null(mesh)==TRUE){
            open3d()
            plotRefToTarget(DA.mns[,,1],DA.mns[,,2],method="points",main="Directional Asymmetry")
            open3d()
            plotRefToTarget(FA.mns[,,1],FA.mns[,,2],method="points",main="Fluctuating Asymmetry")
          } 
          if(is.null(mesh)==FALSE){
            plotRefToTarget(DA.mns[,,1],DA.mns[,,2],mesh,method="surface")
            title3d(main="Directional Asymmetry")
            plotRefToTarget(FA.mns[,,1],FA.mns[,,2],mesh,method="surface")
            title3d(main="Fluctuating Asymmetry")
          }  
        }
      layout(1) 
      } 
    if(verbose==TRUE){
    	class(res.shape) <- c("anova", class(res.shape))
    	return(list(symm.shape=symm.component,asymm.shape=asymm.component, 
                                ANOVA.Shape=res.shape)) }
    if(verbose==FALSE){
    	class(res.shape) <- c("anova", class(res.shape))
    	return(ANOVA.Shape=res.shape) }
  }
}
