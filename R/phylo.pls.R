#' Quantify phylogenetic morphological integration between two sets of variables
#' 
#' Function quantifies the degree of phylogenetic morphological covariation between two sets of
#' Procrustes-aligned coordinates using partial least squares. 
#' 
#' The function quantifies the degree of phylogenetic morphological integration between two sets of shape data as 
#'   defined by landmark coordinates. It is assumed that the landmarks have previously been aligned using 
#'   Generalized Procrustes Analysis (GPA) [e.g., with \code{\link{gpagen}}].   
#' 
#' The function estimates the degree of  morphological covariation between two sets of variables 
#' while accounting for phylogeny using partial least squares (Adams and Felice 2014). The observed value is statistically assessed 
#' using permutation, where data for one block are permuted across the tips of the phylogeny, 
#' an estimate of the covariation between sets of variables, and compared to the observed value. 
#' 
#'   A plot of PLS scores from Block1 versus Block2 is provided for the first set of PLS axes. Thin-plate spline 
#'   deformation grids along these axes are also shown (if data were input as a 3D array).
#' 
#' @param A1 A 2D array (n x [p1 x k1]) or 3D array (p1 x k1 x n) containing landmark coordinates for the first block
#' @param A2 A 2D array (n x [p2 x k2]) or 3D array (p2 x k2 x n) containing landmark coordinates for the second block 
#' @param phy A phylogenetic tree of {class phylo} - see \code{\link[ape]{read.tree}} in library ape
#' @param warpgrids A logical value indicating whether deformation grids for shapes along PC1 should be displayed
#'  (only relevant if data for A1 or A2 [or both] were input as 3D array)
#' @param iter Number of iterations for significance testing
#' @param label An optional vector indicating labels for each specimen that are to be displayed
#' @param verbose A logical value indicating whether the output is basic or verbose (see Value below)
#' @export
#' @keywords analysis
#' @author Dean Adams
#' @return Function returns a list with the following components: 
#'   \item{PLS Correlation}{The estimate of phylogenetic morphological covariation}
#'   \item{pvalue}{The significance level of the observed signal}
#'   \item{Block 1 PLS Scores}{PLS scores for the first block of landmarks (when {verbose=TRUE})}
#'   \item{Block 2 PLS Scores}{PLS scores for the second block of landmarks (when {verbose=TRUE})}
#' @references  Adams, D.C. and R. Felice. 2014. Assessing phylogenetic morphological 
#' integration and trait covariation in morphometric data using evolutionary covariance 
#' matrices. PLOS ONE. 9(4):e94335.
#' @examples
#' 
#' data(plethspecies)
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' 
#' phylo.pls(Y.gpa$coords[1:5,,],Y.gpa$coords[6:11,,],plethspecies$phy,iter=99)

phylo.pls <-function(A1, A2, phy, warpgrids=TRUE,iter=999, label=NULL,verbose=FALSE){ 
  if(any(is.na(A1))==T){
    stop("Data matrix 1 contains missing values. Estimate these first(see 'estimate.missing').")  } 
  if(any(is.na(A2))==T){
    stop("Data matrix 2 contains missing values. Estimate these first(see 'estimate.missing').") }
  if (is.numeric(A1)== FALSE){stop("A1 is not numeric. see ?numeric")}
  if (is.numeric(A2)== FALSE){stop("A2 is not numeric. see ?numeric")}
  if (class(phy) != "phylo")
   stop("phy must be of class 'phylo.'") 
  if (length(dim(A1)) == 3){x<-two.d.array(A1)}
  if (length(dim(A1)) == 2){x<-A1}
  if (length(dim(A2)) == 3){y<-two.d.array(A2)}
  if (length(dim(A2)) == 2){y<-A2}
  num.taxa.X<-nrow(x)
  namesX<-rownames(x)
  num.taxa.Y<-nrow(y)
  namesY<-rownames(y)
  if (is.null(namesX)){
    stop("No specimen names in data matrix 1. Please assign specimen names.")  } 
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && length(phy$tip.label) < num.taxa.X)
    stop("Tree is missing some taxa present in the data matrix") 
  if (length(match(phy$tip.label, namesX)) != num.taxa.X && num.taxa.X < length(phy$tip.label)) 
    stop("Tree contains some taxa not present in present in the data matrix")  
  if (is.null(namesY)){
    stop("No specimen names in data matrix 2. Please assign specimen names")  } 
  if (is.null(namesX) == FALSE && is.null(namesY) == FALSE) {
    mtch.A <- namesX[is.na(match(namesX, namesY))]
    if (length(mtch.A) > 0) {
      stop("Specimen names in data sets are not the same.")   } 
  }
  mtch.B <- namesX[is.na(match(namesX, phy$tip.label))]
  if (length(mtch.B) > 0) {
    stop("Taxa labels on tree and taxa matrix are not the same.")} 
  data.all<-cbind(x,y) 
  Nspec<-nrow(x) 
  C<-vcv.phylo(phy,anc.nodes=FALSE) 
  C<-C[rownames(y),rownames(y)] 
  x<-x[rownames(y),]  
  det.C<-det(C)
  if(det.C>0){invC<-solve(C)}
  if(det.C==0){svd.C<-svd(C)
               Positive <- svd.C$d > max(1e-08 * svd.C$d[1L], 0)
               invC<- svd.C$v[, Positive, drop = FALSE] %*% ((1/svd.C$d[Positive]) *t(svd.C$u[, Positive, drop = FALSE]))}
  one<-matrix(1,Nspec,1)  
  a<-t(t(one)%*%invC%*% data.all)*sum(sum(invC))^-1  
  R<- t(data.all-one%*%t(a))%*%invC%*%(data.all-one%*%t(a))*(Nspec-1)^-1 
  R12<- R[1:dim(x)[2], (dim(x)[2] + 1):(dim(x)[2] +  dim(y)[2])] 
  pls <- svd(R12) 
  U <- pls$u 
  V <- pls$v 
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if(any(lambda == 0)){
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[,1:(length(lambda))]
  D.mat <- solve(eigC.vect%*% diag(sqrt(lambda)) %*% t(eigC.vect))
  Phy.X<-D.mat%*%( data.all-one%*%t(a)) 
  x.phy <- Phy.X[, c(1:dim(x)[2])] 
  y.phy <- Phy.X[, c((dim(x)[2] + 1):(dim(x)[2] +  dim(y)[2]))] 
  XScores <- x.phy %*% U 
  YScores <- y.phy %*% V
  pls.obs <- cor(XScores[, 1], YScores[, 1]) 
  P.val <- 1
  pls.val <- rep(0, iter)
  for (ii in 1:iter) {
    y.r <- y[sample(nrow(y)), ]
    data.all.r<-cbind(x,y.r)
    a.r<-t(t(one)%*%invC%*%data.all.r)*sum(sum(invC))^-1  
    R.r<- t(data.all.r-one%*%t(a.r))%*%invC%*%(data.all.r-one%*%t(a.r))*(Nspec-1)^-1  
    R12.r <- R.r[1:dim(x)[2], (dim(x)[2] + 1):(dim(x)[2] +  dim(y.r)[2])]
    pls.r <- svd(R12.r)
    U.r <- pls.r$u
    V.r <- pls.r$v
    Phy.X.r<-D.mat%*%( data.all.r-one%*%t(a.r))
    x.phy.r <- Phy.X.r[, c(1:dim(x)[2])]
    y.phy.r <- Phy.X.r[, c((dim(x)[2] + 1):(dim(x)[2] +  dim(y.r)[2]))]
    XScores.r <- x.phy.r %*% U.r[, 1]
    YScores.r <- y.phy.r %*% V.r[, 1]
    pls.r <- cor(XScores.r, YScores.r)
    pls.val[ii] <- pls.r
    P.val <- ifelse(pls.r >= pls.obs, P.val + 1, P.val)
  }
  pls.val[iter + 1] = pls.obs
  P.val <- P.val/(iter + 1) 
  if (length(dim(A1))==2 && length(dim(A2))==2){
    plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS Plot",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
    if(length(label!=0)){text(XScores[,1],YScores[,1],label,adj=c(-.7,-.7))}
  }
  if (length(dim(A1))==3){A1.ref<-mshape(A1);
                          pls1.min<-A1[,,which.min(XScores[,1])];pls1.max<-A1[,,which.max(XScores[,1])]}
  if (length(dim(A2))==3){A2.ref<-mshape(A2);
                          pls2.min<-A2[,,which.min(XScores[,1])];pls2.max<-A2[,,which.max(XScores[,1])]}
  if (dim(A1)[2] == 2 ||dim(A2)[2] == 2) {par(mar=c(1,1,1,1)+0.1)
    split.screen(matrix(c(0.22,1,0.22,1,.19,.39,0,.19,.8,1,0,.19,0,.19,.19,.39,0,.19,.8,1), byrow=T, ncol=4))
    screen(1)
      plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS1 Plot: Block 1 (X) vs. Block 2 (Y) ",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
      if(length(label!=0)){text(XScores[,1],YScores[,1],label,adj=c(-.7,-.7))}
    if(warpgrids==TRUE){
      if (length(dim(A1))==3  && dim(A1)[2]==2){
        screen(2);       tps(A1.ref, pls1.min, 20,sz=.7)
        screen(3);       tps(A1.ref, pls1.max, 20,sz=.7)     
      }
      if (length(dim(A2))==3  && dim(A2)[2]==2){
        screen(4);       tps(A2.ref, pls2.min, 20,sz=.7)
        screen(5);       tps(A2.ref, pls2.max, 20,sz=.7)
      }
    }
    close.screen(all.screens=TRUE)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  if (length(dim(A1))==3  && dim(A1)[2] == 3) {
    plot(XScores[,1],YScores[,1],pch=21,bg="black",main="PLS Plot",xlab = "PLS1 Block 1",ylab = "PLS1 Block 2")
    if(length(label!=0)){text(XScores[,1],YScores[,1],label,adj=c(-.7,-.7))}
    open3d()
    plot3d(pls1.min, type = "s", col = "gray", main = paste("PLS Block1 negative"),size = 1.25, aspect = FALSE)
    open3d()
    plot3d(pls1.max, type = "s", col = "gray", main = paste("PLS Block1 positive"),size = 1.25, aspect = FALSE)
  }  
  if (length(dim(A2))==3  && dim(A2)[2] == 3){
    open3d()
    plot3d(pls2.min, type = "s", col = "gray", main = paste("PLS Block2 negative"),size = 1.25, aspect = FALSE)
    open3d()
    plot3d(pls2.max, type = "s", col = "gray", main = paste("PLS Block2 positive"),size = 1.25, aspect = FALSE)
  }
  if(verbose==TRUE){return(list("PLS Correlation" = pls.obs, pvalue = P.val, 
                                "Block 1 PLS Scores" = XScores[,1], "Block 2 PLS Scores" = YScores[,1]))}
  if(verbose==FALSE){return(list("PLS Correlation" = pls.obs, pvalue = P.val)) }
}

