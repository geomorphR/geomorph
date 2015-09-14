#'  Quantify and compare shape change trajectories
#'
#'  Function estimates attributes of shape change trajectories or motion trajectories for a set of 
#'  Procrustes-aligned specimens and compares them statistically
#'
#'  The function quantifies phenotypic shape change trajectories from a set of specimens, and assesses variation 
#'  in these parameters via permutation. A shape change trajectory is defined by a sequence 
#'  of shapes in tangent space. These trajectories can be quantified for various attributes (their size, orientation, 
#'  and shape), and comparisons of these attribute enable the statistical comparison of shape change 
#'  trajectories (see Collyer and Adams 2013; Collyer and Adams 2007; Adams and Collyer 2007; Adams and Collyer 2009). 
#'
#'  Data input is specified by a formula (e.g., Y~X), where 'Y' specifies the response variables (trajectory data), 
#'  and 'X' contains one or more independent variables (discrete or continuous). The response matrix 'Y' must be 
#'  in the form of a two-dimensional data matrix of dimension (n x [p x k]), rather than a 3D array. The function
#'  \code{\link{two.d.array}} can be used to obtain a two-dimensional data matrix from a 3D array of landmark
#'  coordinates. It is assumed that the order of the specimens 'Y' matches the order of specimens in 'X'. 
#' 
#'  There are two primary modes of analysis through this function. If "estimate.traj=TRUE" the function 
#'  estimates shape trajectories using the least-squares means for groups, based on a two-factor model
#'  (e.g., Y~A+B+A:B). Under this implementation, the last factor in 'X' must be the interaction term, and
#'  the preceding two factors must be the effects of interest. Covariates may be included in 'X', and must
#'  precede the factors of interest (e.g., Y~cov+A*B). In this implementation, 'Y' contains a matrix of landmark
#'  coordinates. It is assumed that the landmarks have previously been aligned using Generalized Procrustes 
#'  Analysis (GPA) [e.g., with \code{\link{gpagen}}]. 
#'
#'  If "estimate.traj=FALSE" the trajectories are assembled directly from the set of shapes provided in 'Y'. 
#'  With this implementation, the user must specify the number of shapes that comprise each trajectory. This 
#'  approach is useful when the set of shapes forming each trajectory have been quantified directly 
#'  (e.g., when motion paths are compared: see Adams and Cerney 2007). With this implementation, variation in 
#'  trajectory size, shape, and orientation are evaluated for each term in 'X'.(see Adams and Cerney 2007). 
#'  Currently only single-factor analyses are supportd with this method.
#'
#'  Once the function has performed the analysis, it generates a plot of the trajectories as visualized in the 
#'  space of principal components (PC1 vs. PC2). The first point in each trajectory is displayed as white, the 
#' last point is black, and any middle points on the trajectories are in gray.  The colors of trajectories follow
#'  the order in which they are found in the dataset as a default, using R's standard color palette: black, red, green3,
#'  blue, cyan, magenta, yellow, and gray. However, one can override these colors with group.cols.
#'
#' @param f1 A formula for the linear model (e.g., y~x1+x2)
#' @param estimate.traj A logical value indicating whether trajectories are estimated from original data; 
#'   described below
#' @param iter Number of iterations for significance testing
#' @param traj.pts An optional value specifying the number of points in each trajectory (if estimate.traj=FALSE)
#' @param verbose A logical indicator for verbose (random) output (observed cases always first)
#' @param group.cols A vector of colors to use for trajectories, in the order of the levels of the grouping variable. E.g., c("red", "blue, "orange",...)
#' @param pca A logical indicator if a principal component analysis should be performed on data
#' @export
#' @keywords analysis
#' @author Dean Adams and Michael Collyer
#' @return If "estimate.traj=TRUE", the function returns a list with the following components: 
#'   \item{anova.table}{Procrustes ANOVA table}
#'   \item{Size$Obs.dif}{A matrix of pairwise differences in trajectory size}
#'   \item{Size$Z}{A matrix of pairwise effect sizes for differences in trajectory size}
#'   \item{Size$P}{A matrix of pairwise significance levels for differences in trajectory size}
#'   \item{Direction$Obs.dif}{A matrix of pairwise differences in trajectory orientation}
#'   \item{Direction$Z}{A matrix of effect sizes for differences in trajectory orientation}
#'   \item{Direction$P}{A matrix of pairwise significance levels for differences in trajectory orientation}
#'   \item{Shape$Obs.dif}{A matrix of pairwise differences in trajectory shape (if applicable)}
#'   \item{Shape$Z}{A matrix of pairwise effect sizes for differences in trajectory shape (if applicable)}
#'   \item{Random.values}{All random values for all RRPP permutations (when {verbose=TRUE})}
#' @return If "estimate.traj=FALSE", the function returns a list with the following components: 
#'   \item{MANOVA.location.covariation}{Procrustes ANOVA table}
#'   \item{ANOVA.Size}{Results of permutational-ANOVA assessing variation in trajectory size}
#'   \item{ANOVA.Dir}{Results of permutational-ANOVA assessing variation in trajectory orientation}
#'   \item{ANOVA.Shape}{Results of permutational-ANOVA assessing variation in trajectory shape (if applicable)}
#'   \item{random.SS.location}{Random SS from RRPP permutations (when {verbose=TRUE})}
#'   \item{random.SS.size}{Random SS from RRPP permutations (when {verbose=TRUE})}
#'   \item{random.SS.dir}{Random SS from RRPP permutations (when {verbose=TRUE})}
#'   \item{random.SS.shape}{Random SS from RRPP permutations (when {verbose=TRUE})}
#'  
#' @references Collyer, M.L., and D.C. Adams. 2013. Phenotypic trajectory analysis: Comparison of 
#'  shape change patterns in evolution and ecology. Hystrix. 24:75-83.
#' @references Adams, D. C. 2010. Parallel evolution of character displacement driven by competitive 
#'   selection in terrestrial salamanders. BMC Evol. Biol. 10:1-10.
#' @references Adams, D. C., and M. M. Cerney. 2007. Quantifying biomechanical motion using Procrustes 
#'   motion analysis. J. Biomech. 40:437-444.
#' @references Adams, D. C., and M. L. Collyer. 2007. The analysis of character divergence along environmental 
#'   gradients and other covariates. Evolution 61:510-515.
#' @references Adams, D. C., and M. L. Collyer. 2009. A general framework for the analysis of phenotypic 
#'   trajectories in evolutionary studies. Evolution 63:1143-1154.
#' @references Collyer, M. L., and D. C. Adams. 2007. Analysis of two-state multivariate phenotypic change 
#'   in ecological studies. Ecology 88:683-692.
#' @examples
#' #1: Estimate trajectories from LS means in 2-factor model
#' data(plethodon) 
#' Y.gpa<-two.d.array(gpagen(plethodon$land)$coords)    
#'
#' trajectory.analysis(Y.gpa~plethodon$species*plethodon$site,iter=15)
#' 
#' # Retaining random values (first sets are always observed)
#' tra <- trajectory.analysis(Y.gpa~plethodon$species*plethodon$site,iter=15, verbose = TRUE)
#' tra$anova.table
#' tra$Random.values
#'
#' #2: Compare motion trajectories
#' data(motionpaths) 
#'
#' #Motion paths represented by 5 time points per motion 
#'
#' trajectory.analysis(motionpaths$trajectories~motionpaths$groups,
#' estimate.traj=FALSE, traj.pts=5,iter=15)
#' 
#' trajectory.analysis(motionpaths$trajectories~motionpaths$groups,
#' estimate.traj=FALSE, traj.pts=5,iter=15, verbose=TRUE)
trajectory.analysis<-function(f1,estimate.traj=TRUE,traj.pts=NULL,iter=999, pca=TRUE,verbose=FALSE, group.cols=NULL){
  pf= procD.fit(f1)
  Y <- pf$Y
  dat <- pf$mf
  if(estimate.traj==TRUE){
    anova.parts.obs <- anova.parts(pf, keep.order=ko)
    anova.tab <-anova.parts.obs$table  
    Xs <- pf$Xs
    k <- length(pf$Terms)
    if(pca==TRUE) y<-prcomp(Y)$x else y <- Y
    if(k < 3) stop("Model does not appear to be factorial.  Check model formula (see help file).") 
    int.term<-grep(":", pf$Terms[k])
    if(int.term!=1) stop("Last col of X-matrix does not contain interaction between main effects (see help file).")        
    p<-ncol(y) 
    n1<-length(levels(dat[,k-1]))
    k1<-length(levels(dat[,k]))
    Plm <-array(, c(k, 1, iter+1))
    SS.obs <-anova.parts.obs$SS[1:k]
    Plm[,,1] <- SS.obs    
    fac12<-single.factor(formula(f1))
    covs <- cov.extract(formula(f1))
    if(ncol(covs) == 0) covs <- NULL else covs <- model.matrix(~covs)
    lsmeans.obs <- ls.means(fac12, cov.mf=covs, y)
    traj.specs.obs<- aperm(array(t(lsmeans.obs), c(p,k1,n1)), c(2,1,3)) 
    trajsize.obs<-trajsize(traj.specs.obs,n1,k1) 
    trajdir.obs<-trajorient(traj.specs.obs,n1,p); diag(trajdir.obs)<-0 
    trajshape.obs<-trajshape(traj.specs.obs) 
    PSize<-POrient<-PShape<-array(,dim=c(n1,n1,iter+1))
    PSize[,,1] <- trajsize.obs
    POrient[,,1] <- trajdir.obs
    PShape[,,1] <- trajshape.obs
    for(i in 1:iter){
      SS.r <- SS.traj.random(pf, Yalt = "RRPP", iter=1)
      Plm[,,i+1] <- SS.r$SS[,,2]
      lsmeans.r <- ls.means(fac12, cov.mf=covs, SS.r$Y)
      traj.specs.r<- aperm(array(t(lsmeans.r), c(p,k1,n1)), c(2,1,3)) 
      trajsize.r<-trajsize(traj.specs.r,n1,k1) 
      trajdir.r<-trajorient(traj.specs.r,n1,p); diag(trajdir.r)<-0 
      trajshape.r<-trajshape(traj.specs.r) 
      PSize[,,i+1] <- trajsize.r
      POrient[,,i+1] <- trajdir.r
      PShape[,,i+1] <- trajshape.r
    }  
    P.val.lm <- Pval.matrix(Plm)
    P.val.size <- Pval.matrix(PSize)
    P.val.dir <- Pval.matrix(POrient)
    P.val.shape <- Pval.matrix(PShape)
    Z.lm <- Effect.size.matrix(Plm)
    Z.size <- Effect.size.matrix(PSize); diag(Z.size) <- 0
    Z.dir <- Effect.size.matrix(POrient); diag(Z.dir) <- 0
    Z.shape <- Effect.size.matrix(PShape); diag(Z.shape) <- 0 
    rownames(P.val.size) <- colnames(P.val.size) <- rownames(P.val.dir) <- colnames(P.val.dir) <- 
      rownames(P.val.shape) <- colnames(P.val.shape) <- rownames(Z.size) <- colnames(Z.size) <-
      rownames(Z.dir) <- colnames(Z.dir) <- rownames(Z.shape) <- colnames(Z.shape) <-
      rownames(trajsize.obs) <- colnames(trajsize.obs) <- rownames(trajdir.obs) <- colnames(trajdir.obs) <-
      rownames(trajshape.obs) <- colnames(trajshape.obs) <- levels(dat[,k-1])
    dimnames(Plm)[[1]] <- dimnames(anova.tab)[[1]][1:k]
    dimnames(PSize)[[1]] <- dimnames(PSize)[[2]] <- 
      dimnames(POrient)[[1]] <- dimnames(POrient)[[2]] <- 
      dimnames(PShape)[[1]] <- dimnames(PShape)[[2]] <- levels(dat[,k-1])
    
    anova.tab <- data.frame(anova.tab, Z = c(Z.lm, NA, NA), P.value = c(P.val.lm, NA, NA))
    anova.title = "\nRandomized Residual Permutation Procedure used\n"
    attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
    class(anova.tab) <- c("anova", class(anova.tab))

    if (k1 > 2) Random.values=list(random.SS=Plm,random.size.dif = PSize, random.dir.dif = POrient, random.shape.dif = PShape) else 
      Random.values=list(random.SS=Plm,random.size.dif = PSize, random.dir.dif = POrient)
    if(verbose==FALSE) results <- list(anova.table = anova.tab, 
         Size=list(Obs.dif=trajsize.obs,Z=Z.size,P = P.val.size),
         Direction=list(Obs.dif=trajdir.obs,Z=Z.dir,P = P.val.dir),
         Shape=list(Obs.dif=trajshape.obs,Z=Z.shape,P = P.val.shape))
    else results <- list(anova.table = anova.tab, 
                         Size=list(Obs.dif=trajsize.obs,Z=Z.size,P = P.val.size),
                         Direction=list(Obs.dif=trajdir.obs,Z=Z.dir,P = P.val.dir),
                         Shape=list(Obs.dif=trajshape.obs,Z=Z.shape,P = P.val.shape),
                         Random.values=Random.values)
    
    trajplot(y,traj.specs.obs, groups = dat[,k-1], group.cols = group.cols)
    if(k1 == 2) return(results[-4]) else 
        return(results)
  }
  
  if(estimate.traj==FALSE){
    if(is.null(traj.pts)) stop("Number of points in the trajectory not specified.") 
    if(length(pf$Terms) > 1) stop("If data are already trajectories, only a single-factor model is currently supported. (See Help file)")
    form.in = formula(f1)
    X <- model.matrix(form.in)
    Xs <- pf$Xs
    k <- length(pf$Terms) 
    k1<-traj.pts
    y <- Y
    n1<-nrow(y)
    p1<-ncol(y)/k1  
    if (p1>2){
      y.2d<-matrix(t(y),ncol=k1,byrow=TRUE)
      y.2d<-prcomp(y.2d)$x
      y<-two.d.array(arrayspecs(y.2d,k1,p1))
    }        
    form.new <- as.formula(paste(c("~", form.in[[3]])))
    traj.specs.obs<-arrayspecs(y,k1,p1) 
    size.obs<-trajsize(traj.specs.obs,n1,k1) 
    dir.obs<-trajorient(traj.specs.obs,n1,p1) 
    diag(dir.obs)<-0
    shape.obs<-trajshape(traj.specs.obs) 
    size.tab<-Hat.anova.tab(size.obs, form.new)
    shape.tab<-Hat.anova.tab(shape.obs, form.new)
    dir.tab<-Hat.anova.tab(dir.obs, form.new)
    size.SS.obs <- Hat.SS.model(Gower.center(size.obs), X)
    shape.SS.obs <- Hat.SS.model(Gower.center(shape.obs), X)
    dir.SS.obs <- Hat.SS.model(Gower.center(dir.obs), X)
    anova.parts.obs <- anova.parts(pf,Yalt="observed")
    SS.obs <-anova.parts.obs$SS[1:k]
    anova.tab <-anova.parts.obs$table  
    Plm <- array(,c(k,1,iter+1))
    Plm[,,1] <- SS.obs
    PSize <- PShape <- POrient <- array(,c(k,k,iter+1))
    PSize[,,1] <- size.SS.obs[1:k]
    PShape[,,1] <- shape.SS.obs[1:k]
    POrient[,,1] <- dir.SS.obs[1:k]
      perm.index <- 1:nrow(y)
      SSY <- SSE(list(lm(y~1)$resdiuals))
      Dy <- as.matrix(dist(y))
      for(i in 1:iter){
        pr <- sample(perm.index)
        Plm[,,i+1] <- SSY - Hat.SSE(Gower.center(Dy[pr,pr]), X)
        size.r<-size.obs[pr,pr]
        dir.r <-dir.obs[pr,pr]
        shape.r<-shape.obs[pr,pr]
        PSize[,,i+1] <- Hat.SS.model(Gower.center(size.r), X)
        PShape[,,i+1] <- Hat.SS.model(Gower.center(shape.r), X)
        POrient[,,i+1] <- Hat.SS.model(Gower.center(dir.r), X)
      }
    
    size.tab <- data.frame(cbind(size.tab, Z=c(Effect.size.matrix(PSize),NA,NA),P.value=c(Pval.matrix(PSize),NA,NA)))
    shape.tab <- data.frame(cbind(shape.tab, Z=c(Effect.size.matrix(PShape),NA,NA),P.value=c(Pval.matrix(PShape),NA,NA)))
    dir.tab <- data.frame(cbind(dir.tab, Z=c(Effect.size.matrix(POrient),NA,NA),P.value=c(Pval.matrix(POrient),NA,NA)))
    
    anova.tab <- data.frame(anova.tab, Z = c(Effect.size.matrix(Plm), NA, NA), P.value = c(Pval.matrix(Plm), NA, NA))
    anova.title = "\nRandomized Residual Permutation Procedure used\n"
    attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n",anova.title)
    class(anova.tab) <- c("anova", class(anova.tab))
    class(size.tab) <- c("anova", class(size.tab))
    class(dir.tab) <- c("anova", class(dir.tab))
    class(shape.tab) <- c("anova", class(shape.tab))
    if(verbose==TRUE) {
      results <- list(MANOVA.location.covariation = anova.tab, 
                      ANOVA.Size=size.tab,
                      ANOVA.Dir = dir.tab,
                      ANOVA.shape=shape.tab,
                      random.SS.location=as.vector(Plm), 
                      random.SS.size=as.vector(PSize), 
                      random.SS.dir=as.vector(POrient), 
                      random.SS.shape=as.vector(PShape))
    } else
    {
      results <- list(MANOVA.location.covariation = anova.tab, 
                      ANOVA.Size=size.tab,
                      ANOVA.Dir = dir.tab,
                      ANOVA.shape=shape.tab)
    }
    
    y.plot<-matrix(t(two.d.array(traj.specs.obs)), ncol=p1,byrow=TRUE)
    if(pca==T) {
      y.plot <-prcomp(y.plot)$x
      traj.specs.obs <- arrayspecs(y.plot, traj.pts, p1)
    }
      
    trajplot(y.plot,traj.specs.obs, groups = as.factor(dat[[2]]))
    if(k1 == 2) return(results[-4]) else return(results)
  }
  results
}
