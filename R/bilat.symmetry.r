#' Analysis of bilateral symmetry
#'
#' Function performs an analysis of directional and fluctuating asymmetry for bilaterally symmetric objects 
#'
#' The function quantifies components of shape variation for a set of specimens as described by their patterns of symmetry
#'  and asymmetry. Here, shape variation is decomposed into variation among individuals, variation among sides (directional 
#'  asymmetry), and variation due to an individual x side interaction (fluctuating symmetry). These components are then 
#'  statistically evaluated using Procrustes ANOVA. Statistical assessment of model effects for shape variation is accomplished using permutation procedures. 
#'  Methods for both matching symmetry and object symmetry can be implemented. Matching symmetry is when each object contains mirrored 
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
#'  Compared to older versions of geomorph, some results can be expected to be slightly different.  Starting with geomorph 3.0,
#'  results use only type I sums of squares (SS) with either full ranomization of raw shape values or RRPP (prefered with nested terms)
#'  for analysis of variance (ANOVA).  Older versions used a combination of parametric and non-parametric results, as well as a combinaton
#'  of type I and type III SS.  While analytical conclusions should be consistent (i.e., "significance" of effects is the same), these
#'  updates maintain consistency in analytical philosophy.
#'  
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{bilat.symmetry}}.
#'
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens [for "object.sym=FALSE, A is of dimension (n x k x 2n)]
#' @param ind A vector containing labels for each individual. For matching symmetry, the matched pairs receive the same 
#' label (replicates also receive the same label).
#' @param side An optional vector (for matching symmetry) designating which object belongs to which 'side-group'
#' @param replicate An optional vector designating which objects belong to which group of replicates.
#' @param object.sym A logical value specifying whether the analysis should proceed based on object symmetry {=TRUE} or matching symmetry {=FALSE}
#' @param land.pairs An optional matrix (for object symmetry) containing numbers for matched pairs of landmarks across the line of symmetry 
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}}. It is imperative
#' that the variables "ind", "side", and "replicate" in the data frame match these names exactly (as shown in examples below).  
#' @param iter Number of iterations for significance testing.
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing.
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return An object of class "bilat.symmetry" returns a list of the following
#' \item{shape.anova}{An analysis of variance table for the shape data.}
#' \item{size.anova}{An analysis of variance table for the shape data (when object.sym=FALSE).}
#' \item{symm.component}{The symmetric component of shape variation.}
#' \item{asym.component}{The asymmetric component of shape variation.}
#' \item{DA.component}{The directional asymmetry component, found as the mean shape for each side.}
#' \item{FA.compnent}{The fluctuating asymmetry component for each specimen, found as the specimen-specific side deviation adjusted for the mean
#'  directional asymmetry in the dataset.}
#' \item{data.type}{A value indicating whether the analysis was performed as Object or Matching symmetry.}
#' \item{permutations}{The number of random permutations used.}
#' \item{random.shape.SS}{A matrix of random outcomes for SS from the Shape analysis.}
#' \item{random.size.SS}{A matrix of random outcomes for SS from the Cntroid Size analysis.}
#' \item{perm.method}{A value indicating whether "Raw" values were shuffled or "RRPP" performed.}
#' \item{call}{The matched call.}
#' 
#' @references Klingenberg, C.P. and G.S. McIntyre. 1998. Quantitative genetics of geometric shape in the mouse mandible. Evolution. 55:2342-2352.
#' @references Mardia, K.V., F.L. Bookstein, and I.J. Moreton. 2000. Statistical assessment of bilateral symmetry of shapes. Biometrika. 87:285-300.
#' @references Klingenberg, C.P., M. Barluenga, and A. Meyer. 2002. Shape analysis of symmetric structures: quantifying variation among
#' individuals and asymmetry. Evolution. 56:1909-1920.
#' @examples
#' #Example of matching symmetry
#'
#' data(mosquito)
#' gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, ind=mosquito$ind, side=mosquito$side,
#' replicate=mosquito$replicate)
#' mosquito.sym <- bilat.symmetry(A = wingshape, ind = ind, side = side,
#' replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 499, data = gdf)
#' summary(mosquito.sym)
#' plot(mosquito.sym, warpgrids = TRUE)
#' mosquito.sym$shape.anova # extract just the anova table on shape
#'
#' #Example of object symmetry
#'
#' data(scallops)
#' gdf <- geomorph.data.frame(shape = scallops$coorddata, ind=scallops$ind)
#' scallop.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
#' land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 499)
#' summary(scallop.sym)
#' # NOTE one can also: plot(scallop.sym, warpgrids = TRUE, mesh = NULL)
#' # NOTE one can also: scallop.sym$data.type # recall the symmetry type

bilat.symmetry<-function(A,ind=NULL,side=NULL,replicate=NULL,object.sym=FALSE,land.pairs=NULL,
      data = NULL, iter=999,RRPP=TRUE){
  if(!is.null(data)) {
    A.name <- deparse(substitute(A))
    A.name.match <- match(names(data), A.name)
    if(all(is.na(A.name.match))) stop("Coordinates are not part of the data frame provided")
    A.name.match <- which(!is.na(A.name.match))
    A <- data[[A.name.match]]
  }
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  if(is.null(data)){
    if(is.null(ind)) stop("Individuals not specified.") else ind <- factor(ind)
    if(!is.null(side)) side <- factor(side)
    if(!is.null(replicate)) replicate <- factor(replicate)
  } else {
    ind.match <- match(names(data), "ind")
    if(all(is.na(ind.match))) stop("Individuals not specified in geomorph data frame")
    ind <- factor(data[[which(!is.na(ind.match))]])
    side.match <- match(names(data), "side")
    if(all(is.na(side.match))) side <- NULL else 
      side <- factor(data[[which(!is.na(side.match))]])
    replicate.match <- match(names(data), "replicate")
    if(all(is.na(replicate.match))) replicate <- NULL else 
      replicate <- factor(data[[which(!is.na(replicate.match))]])
  }
  n<-dim(A)[3]; k<-dim(A)[2]; p<-dim(A)[1]; nind<-nlevels(ind); spec.names<-dimnames(A)[[3]]
  if(object.sym == FALSE && is.null(side)) stop("Sides not specified.") 
  if(object.sym==TRUE){
    if(is.null(land.pairs)){stop("Landmark pairs not specified.")} 
    npairs <- nrow(land.pairs); nl<-p-2*npairs
    A2 <- A
    A2[land.pairs[,1],,] <- A[land.pairs[,2],,]
    A2[land.pairs[,2],,] <- A[land.pairs[,1],,]
    A2[,1,]<-A2[,1,]*-1
    A <- array(c(A,A2), c(p,k, 2*n))
    ind <- factor(rep(ind,2)); side <- gl(2,n); if(!is.null(replicate)) replicate <- rep(replicate,2)
  }
  gpa.res <- gpagen(A)
  Y <- two.d.array(gpa.res$coords)
  if(!is.null(replicate)) {
    form.shape <- Y ~ ind*side + ind:side:replicate 
    dat.shape <- geomorph.data.frame(Y = Y, ind = ind, side = side, replicate = replicate)
  } else {
    form.shape <- Y ~ ind*side
    dat.shape <- geomorph.data.frame(Y = Y, ind = ind, side = side)
  }
  pfitSh <- procD.fit(form.shape, data = dat.shape, keep.order = TRUE)
  kSh <- length(pfitSh$term.labels)
  if(RRPP == TRUE) PSh <- SS.iter(pfitSh,Yalt="RRPP", iter=iter) else PSh <- SS.iter(pfitSh, Yalt="resample", iter=iter)
  anova.parts.Sh <- anova.parts.symmetry(pfitSh, PSh, object.sym)
  anova.tab <-anova.parts.Sh$anova.table 
     if(is.matrix(PSh)){
       P.val <- apply(PSh,1,pval)
       rownames(PSh) <- pfitSh$term.labels
       colnames(PSh) <- c("obs", paste("iter", 1:iter, sep=":"))
     } else {
       P.val <- pval(PSh)
       names(PSh) <-c("obs", paste("iter", 1:iter, sep=":"))
     }
  P.val[kSh] = NA; if(kSh==3) P.val[3] = NA  #b/c nested
  anovaSh <- data.frame(anova.tab, Pr = c(P.val))
  
  if(object.sym==FALSE){  
    size <- gpa.res$Csize
    if(!is.null(replicate)) {
      form.size <- size~ind*side+ind:side:replicate 
      dat.size <- geomorph.data.frame(size = size, ind = ind, side = side, replicate = replicate)
    } else {
      form.size <- size~ind*side
      dat.size <- geomorph.data.frame(size = size, ind = ind, side = side)
    }
    pfitSz=procD.fit(form.size, data=dat.size, keep.order=TRUE)
    if(RRPP == TRUE) PSz <- SS.iter(pfitSz,Yalt="RRPP", iter=iter) else PSz <- SS.iter(pfitSz, Yalt="resample", iter=iter)
    anova.parts.Sz <- anova.parts.symmetry(pfitSz, PSz,object.sym)
    anova.tab.Sz <-anova.parts.Sz$anova.table 
    if(is.matrix(PSz)){
      P.val.Sz <- apply(PSz,1,pval)
      rownames(PSz) <- pfitSz$term.labels
      colnames(PSz) <- c("obs", paste("iter", 1:iter, sep=":"))
    } else {
      P.val.Sz <- pval(PSz)
      names(PSz) <-c("obs", paste("iter", 1:iter, sep=":"))
    }
    P.val.Sz[kSh]=NA; if(kSh==3) P.val.Sz[3]=NA     #b/c nested
    anovaSz <- data.frame(anova.tab.Sz, Pr = c(P.val.Sz))
  }
# build shape components for output
  if(object.sym==FALSE){
    X.ind <- model.matrix(~ind + 0, data = dat.shape)
    symm.component <- arrayspecs(coef(.lm.fit(X.ind, Y)),p,k)
    X.side <- model.matrix(~(side:ind) + 0, data = dat.shape)
    avg.side.symm <- coef(.lm.fit(X.side, Y))
    n.ind <- nlevels(ind)
    n.side <- nlevels(side)
    indsq <- seq(n.side, (n.ind*n.side), n.side)
    asymm.component <- avg.side.symm[indsq,] - avg.side.symm[-indsq,]
    mn.shape <- gpa.res$consensus
    asymm.component<-simplify2array(lapply(1:n.ind, function(j) 
      {t(matrix(asymm.component[j,],k,p)) + mn.shape}))
#    dimnames(asymm.component)[[3]] <- dimnames(symm.component)[[3]] <- spec.names

  }
  if(object.sym==TRUE){
    X.ind <- model.matrix(~ind + 0, data = dat.shape)
    symm.component <- arrayspecs(coef(.lm.fit(X.ind, Y)),p,k)
    mn.shape<-gpa.res$consensus
    n.ind <- nlevels(ind)
    asymm.component <-simplify2array(lapply(1:n.ind, function(j) 
    {mn.shape + gpa.res$coords[,,j] - symm.component[,,j]}))
#    dimnames(asymm.component)[[3]] <- dimnames(symm.component)[[3]] <- spec.names
  }
  X.side <- model.matrix(~side + 0, data = dat.shape)
  DA.mns <- arrayspecs(coef(.lm.fit(X.side, Y)),p,k)
  X.ind.side <- model.matrix(~(side:ind) + 0, data = dat.shape)
  ind.mns <- coef(.lm.fit(X.ind.side, Y))
  mn.DA <- DA.mns[,,1] - DA.mns[,,2]
  n.ind <- nlevels(ind)
  n.side <- nlevels(side)
  indsq <- seq(n.side, (n.ind*n.side), n.side)
  FA.component <- ind.mns[indsq,] - ind.mns[-indsq,]
  mn.shape<-gpa.res$consensus
  FA.component<-simplify2array(lapply(1:n.ind, function(j) 
  {t(matrix(FA.component[j,],k,p)) + mn.shape - mn.DA}))
#  dimnames(FA.component)[[3]] <- spec.names
  colnames(anovaSh)[1] <- "Df"
  colnames(anovaSh)[ncol(anovaSh)] <- "Pr(>F)"
  class(anovaSh) <- c("anova", class(anovaSh))
  if(object.sym==FALSE){
    colnames(anovaSz)[1] <- "Df"
    colnames(anovaSz)[ncol(anovaSz)] <- "Pr(>F)"
    class(anovaSz) <- c("anova", class(anovaSz))
    out<-list(size.anova = anovaSz, shape.anova = anovaSh, symm.component = symm.component,
              asymm.component = asymm.component, DA.component = DA.mns, FA.compnent = FA.component,
              data.type = ifelse(object.sym==TRUE,"Object", "Matching"),
              FA.mns = FA.component, DA.mns = DA.mns,
              permutations = iter+1,
              random.shape.SS = PSh, random.size.SS = PSz,
              perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"),
              call=match.call()) }
  if(object.sym==TRUE){
    out<-list(size.anova = NULL, shape.anova = anovaSh, symm.component = symm.component,
              asymm.component = asymm.component, DA.component = DA.mns, FA.compnent = FA.component,
              data.type = ifelse(object.sym==TRUE,"Object", "Matching"),
              FA.mns = FA.component, DA.mns = DA.mns,
              permutations = iter+1,
              random.shape.SS = PSh, random.size.SS = NULL,
              perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"),
              call=match.call()) }
  class(out) <- "bilat.symmetry"
  out  
}
