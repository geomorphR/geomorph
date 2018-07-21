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
#'  about a midline (e.g., right and left sides of human faces). Details on general approaches for the study of symmetry in geometric 
#'  morphometrics may be found in: Mardia et al. 2000; Klingenberg et al. 2002. 
#'  
#'  As input, the function receives either A 3D array (p x k x n) containing raw landmarks (requiring 
#'  GPA to be performed) or a gpagen object (if GPA has been previously performed) or a geomorphShapes object.
#'  If one wishes to incorporate semilandmarks, GPA can either be performed first using gpagen, or within bilat.symmetry
#'  by passing adequate GPA arguments (i.e. curves, surfaces, ProcD etc, see \code{\link{gpagen}}. 
#'  For "object.sym = FALSE, landmarks should be of dimension (p x k x 2n), as each specimen is 
#'  represented by both left and right configurations.
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
#'   \subsection{Notes for geomorph 3.0}{ 
#'  Compared to older versions of geomorph, some results can be expected to be slightly different.  Starting with geomorph 3.0,
#'  results use only type I sums of squares (SS) with either full randomization of raw shape values or RRPP (preferred with nested terms)
#'  for analysis of variance (ANOVA).  Older versions used a combination of parametric and non-parametric results, as well as a combination
#'  of type I and type III SS.  While analytical conclusions should be consistent (i.e., "significance" of effects is the same), these
#'  updates maintain consistency in analytical philosophy.  This change will require longer computation time for large datasets, but the trade-off
#'  allows users to have more flexibility and eliminates combining disparate analytical philosophies. 
#'  
#'  Note also that significance of terms in the 
#'  model are found by comparing F-values for each term to those obtained via permutation.  F-ratios and df are not strictly necessary (a ratio of SS would suffice), 
#'  but they are reported as is standard for anova tables. Additionally, users will notice that the df reported are based on the number of observations rather than 
#'  a combination of objects * coordinates * dimensions, as is sometimes found in morphometric studies of symmetry. However, this change has no effect 
#'  on hypothesis testing, as only SS vary among permutations (df, coordinates, and dimensions are constants). 
#'  }
#'  
#'  The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work with \code{\link{bilat.symmetry}}.
#'
#' @param A Either A 3D array (p x k x n) containing raw landmarks (requiring GPA to be performed) or a gpagen object (if GPA has been previously performed).  If one 
#' wishes to incorporate semilandmarks, GPA should be performed first using \code{\link{gpagen}}.  Otherwise, bilat.symmetry can perform the initial GPA, assuming all landmarks
#' are fixed.  For "object.sym = FALSE, landmarks should be of dimension (p x k x 2n), as each specimen is represented by both left and right configurations.
#' @param ind A vector containing labels for each individual. For matching symmetry, the matched pairs receive the same 
#' label (replicates also receive the same label).
#' @param side An optional vector (for matching symmetry) designating which object belongs to which 'side-group'
#' @param replicate An optional vector designating which objects belong to which group of replicates.
#' @param object.sym A logical value specifying whether the analysis should proceed based on object symmetry {=TRUE} or matching symmetry {=FALSE}
#' @param land.pairs An optional matrix (for object symmetry) containing numbers for matched pairs of landmarks across the line of symmetry 
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}}. It is imperative
#' that the variables "ind", "side", and "replicate" in the data frame match these names exactly (as shown in examples below).  
#' @param iter Number of iterations for significance testing.
#' @param seed An optional argument for setting the seed for random permutations of the resampling procedure.  
#' If left NULL (the default), the exact same P-values will be found for repeated runs of the analysis (with the same number of iterations).
#'   If seed = "random", a random seed will be used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param RRPP A logical value indicating whether residual randomization should be used for significance testing.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to the screen.  
#' This is helpful for long-running analyses.
#' @param ... Arguments to pass onto gpagen
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return An object of class "bilat.symmetry" returns a list of the following
#' \item{shape.anova}{An analysis of variance table for the shape data.}
#' \item{size.anova}{An analysis of variance table for the shape data (when object.sym=FALSE).}
#' \item{symm.shape}{The symmetric component of shape variation.}
#' \item{asymm.shape}{The asymmetric component of shape variation.}
#' \item{DA.component}{The directional asymmetry component, found as the mean shape for each side.}
#' \item{FA.component}{The fluctuating asymmetry component for each specimen, found as the specimen-specific side deviation adjusted for the mean
#'  directional asymmetry in the dataset.}
#' \item{data.type}{A value indicating whether the analysis was performed as Object or Matching symmetry.}
#' \item{permutations}{The number of random permutations used.}
#' \item{random.shape.F}{A matrix of random F-values from the Shape analysis.}
#' \item{random.size.F}{A matrix of random F-values from the Centroid Size analysis.}
#' \item{perm.method}{A value indicating whether "Raw" values were shuffled or "RRPP" performed.}
#' \item{procD.lm.shape}{A list of typical output from an object of class procD.lm, for shape}
#' \item{procD.lm.size}{If applicable, a list of typical output from an object of class procD.lm, for size.}
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
#' # Previous example, performing GPA first
#' Y.gpa <- gpagen(mosquito$wingshape)
#' mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = ind, side = side,
#' replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 499, data = gdf)
#' summary(mosquito.sym2)
#' summary(mosquito.sym) # same results
#'
#' #Example of object symmetry
#'
#' data(scallops)
#' gdf <- geomorph.data.frame(shape = scallops$coorddata, ind=scallops$ind)
#' scallop.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
#' land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 499)
#' summary(scallop.sym)
#' 
#' # Previous example, incorporating semilandmarks
#' 
#' scallop.sym <- bilat.symmetry(A = shape, ind = ind, object.sym = TRUE, 
#' curves= scallops$curvslide, surfaces = scallops$surfslide,
#' land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 499)
#' summary(scallop.sym)
#' # NOTE one can also: plot(scallop.sym, warpgrids = TRUE, mesh = NULL)
#' # NOTE one can also: scallop.sym$data.type # recall the symmetry type

bilat.symmetry<-function(A, ind=NULL, side=NULL, replicate=NULL, object.sym=FALSE, land.pairs=NULL,
                         data = NULL, iter = 999, seed = NULL, RRPP = TRUE, print.progress = TRUE, ...){

  if(!is.null(data)){
    data <- droplevels(data)
    A.name <- deparse(substitute(A))
    A.name.match <- match(A.name, names(data))[1]
    if(is.na(A.name.match)) A.name.match <- NULL
    } else A.name.match <- NULL
  
    if(is.null(A.name.match) && is.gpagen(A)) {
    size <- A$Csize
    A <- A$coords
    } else if(is.null(A.name.match) && inherits(A, "geomorphShapes")) {
      A <- gpagen(A, ...)
      size <- A$Csize
      A <- A$coords
    } else {
      if(!is.null(data)) {
        if(is.null(A.name.match)) stop("Coordinates are not part of the data frame provided")
          A <- data[[A.name.match]]
          if(any(is.na(A))) stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") 
          if (length(dim(A))!=3) stop("Data matrix not a 3D array (see 'arrayspecs').") 
          if(print.progress) cat("\nInitial GPA\n")
          A <- gpagen(A, print.progress = print.progress, ...)
          size <- A$Csize
          A <- A$coords
        } else {
          if(any(is.na(A))) stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") 
          if (length(dim(A))!=3) stop("Data matrix not a 3D array (see 'arrayspecs').") 
          if(print.progress) cat("\nInitial GPA\n")
          A <- gpagen(A, print.progress = print.progress, ...)
          size <- A$Csize
          A <- A$coords
        }
      }
  
  if(is.null(data)){
    if(is.null(ind)) stop("Individuals not specified.") else ind <- factor(ind, levels = unique(ind))
    if(!is.null(side)) side <- factor(side)
    if(!is.null(replicate)) replicate <- factor(replicate)
  } else {
    ind.match <- match(names(data), "ind")
    if(all(is.na(ind.match))) stop("Individuals not specified in geomorph data frame")
    ind <- factor(data[[which(!is.na(ind.match))]])
    ind <- factor(ind, levels = unique(ind))  
    side.match <- match(names(data), "side")
    if(all(is.na(side.match))) side <- NULL else 
      side <- factor(data[[which(!is.na(side.match))]])
    replicate.match <- match(names(data), "replicate")
    if(all(is.na(replicate.match))) replicate <- NULL else 
      replicate <- factor(data[[which(!is.na(replicate.match))]])
  }
  
  n<-dim(A)[3]; k<-dim(A)[2]; p<-dim(A)[1]; nind<-nlevels(ind); spec.names<-dimnames(A)[[3]]
  if(object.sym == FALSE && is.null(side)) stop("Sides not specified.") 
  if(object.sym == TRUE){
    if(is.null(land.pairs)) {stop("Landmark pairs not specified.")} 
    npairs <- nrow(land.pairs); nl <- p-2*npairs
    A2 <- A
    A2[land.pairs[,1],,] <- A[land.pairs[,2],,]
    A2[land.pairs[,2],,] <- A[land.pairs[,1],,]
    A2[,1,] <- A2[,1,]*(-1)
    A <- array(c(A,A2), c(p,k, 2*n))
    ind <- factor(rep(ind,2)); side <- gl(2,n); if(!is.null(replicate)) replicate <- rep(replicate,2)
    if(print.progress) cat("\nObject Symmetry GPA\n")
    gpa.res <- gpagen(A, print.progress = print.progress)
    A <- gpa.res$coords
  }
  Y <- two.d.array(A)
  if(!is.null(replicate)) {
    form.shape <- Y ~ ind + side + ind/side 
    form.names <- c("ind", "side", "ind:side", "ind:side:replicate", "Total")
    dat.shape <- geomorph.data.frame(Y = Y, ind = ind, side = side, replicate = replicate)
  } else {
    form.shape <- Y ~ ind + side 
    form.names <- c("ind", "side", "ind:side", "Total")
    dat.shape <- geomorph.data.frame(Y = Y, ind = ind, side = side)
  }
  if(print.progress) cat("\nShape Analysis")
  PSh <- procD.lm(form.shape, data = dat.shape, RRPP = RRPP, 
            seed = seed, iter= iter, print.progress = print.progress, 
            effect.type = "F")
  random.shape.F <- PSh$random.F
  if(length(form.names) > 4) {
    SS <- PSh$random.SS
    df.mat <- matrix(PSh$df, NROW(SS), NCOL(SS))
    MS <- SS/df.mat
    random.shape.F[1,] <- MS[1,]/MS[3,]
    random.shape.F[2,] <- MS[2,]/MS[3,]
    random.shape.F[3,] <- MS[3,]/MS[4,]
    PSh$random.F <- random.shape.F
    newZ <- apply(log(random.shape.F + 0.000001), 1, effect.size)
    PSh$aov.table$F[1:3] <- random.shape.F[1:3, 1]
    PSh$aov.table$Z[1:3] <- newZ[1:3]
    rownames(PSh$aov.table) <- form.names
  }
  shape.anova <- PSh$aov.table
  if(object.sym==FALSE){  
    if(!is.null(replicate)) {
      form.size <- size ~ ind + side + ind/side 
      dat.size <- geomorph.data.frame(size = size, ind = ind, side = side, replicate = replicate)
    } else {
      form.size <- size ~ ind + side 
      dat.size <- geomorph.data.frame(size = size, ind = ind, side = side)
    }
    if(print.progress) cat("\nSize Analysis")
    PSz <- procD.lm(form.size, data = dat.size, RRPP = RRPP, 
                    seed = seed, iter= iter, print.progress = print.progress, 
                    effect.type = "F")
    random.size.F <- PSz$random.F
    if(length(form.names) > 4) {
      SS <- PSz$random.SS
      df.mat <- matrix(PSz$df, NROW(SS), NCOL(SS))
      MS <- SS/df.mat
      random.size.F[1,] <- MS[1,]/MS[3,]
      random.size.F[2,] <- MS[2,]/MS[3,]
      random.size.F[3,] <- MS[3,]/MS[4,]
      PSz$random.F <- random.size.F
      newZ <- apply(log(random.size.F + 0.000001), 1, effect.size)
      PSz$aov.table$F[1:3] <- random.size.F[1:3, 1]
      PSz$aov.table$Z[1:3] <- newZ[1:3]
      rownames(PSz$aov.table) <- form.names
    }
    size.anova <- PSz$aov.table
  }
  # build shape components for output
  if(object.sym==FALSE){
    X.ind <- model.matrix(~ind + 0, data = as.data.frame(dat.shape[-1]))
    symm.component <- arrayspecs(coef(lm.fit(X.ind, Y)),p,k)
    dimnames(symm.component)[[3]] <- substr(dimnames(symm.component)[[3]], start=4,
               stop=nchar(dimnames(symm.component)[[3]]))
    X.side <- model.matrix(~(side:ind) + 0, data = as.data.frame(dat.shape[-1]))
    avg.side.symm <- coef(lm.fit(X.side, Y))
    n.ind <- nlevels(ind)
    n.side <- nlevels(side)
    indsq <- seq(n.side, (n.ind*n.side), n.side)
    asymm.component <- avg.side.symm[indsq,] - avg.side.symm[-indsq,]
    mn.shape <- mshape(A)
    asymm.component<-simplify2array(lapply(1:n.ind, function(j) 
    {t(matrix(asymm.component[j,],k,p)) + mn.shape}))
    dimnames(asymm.component)[[3]] <- dimnames(symm.component)[[3]]
    
  }
  if(object.sym==TRUE){
    X.ind <- model.matrix(~ind + 0, data = as.data.frame(dat.shape[-1]))
    symm.component <- arrayspecs(coef(lm.fit(X.ind, Y)),p,k)
    dimnames(symm.component)[[3]] <- substr(dimnames(symm.component)[[3]], start=4,
                stop=nchar(dimnames(symm.component)[[3]]))
    mn.shape<-mshape(A)
    n.ind <- nlevels(ind)
    asymm.component <-simplify2array(lapply(1:n.ind, function(j) 
    {mn.shape + A[,,j] - symm.component[,,j]}))
    dimnames(asymm.component)[[3]] <- dimnames(symm.component)[[3]] 
  }
  X.side <- model.matrix(~side + 0, data = as.data.frame(dat.shape[-1]))
  DA.mns <- arrayspecs(coef(.lm.fit(X.side, Y)),p,k)
  X.ind.side <- model.matrix(~(side:ind) + 0, data = as.data.frame(dat.shape[-1]))
  ind.mns <- coef(.lm.fit(X.ind.side, Y))
  mn.DA <- DA.mns[,,1] - DA.mns[,,2]
  n.ind <- nlevels(ind)
  n.side <- nlevels(side)
  indsq <- seq(n.side, (n.ind*n.side), n.side)
  FA.component <- ind.mns[indsq,] - ind.mns[-indsq,]
  mn.shape<-mshape(A)
  FA.component<-simplify2array(lapply(1:n.ind, function(j) 
  {t(matrix(FA.component[j,],k,p)) + mn.shape - mn.DA}))
  dimnames(FA.component)[[3]] <- dimnames(symm.component)[[3]] 

  if(object.sym==FALSE){
    out<-list(size.anova = size.anova, shape.anova = shape.anova, symm.shape = symm.component,
              asymm.shape = asymm.component, DA.component = DA.mns, FA.component = FA.component,
              data.type = ifelse(object.sym==TRUE,"Object", "Matching"),
              FA.mns = FA.component, DA.mns = DA.mns,
              permutations = iter+1,
              random.shape.F = random.shape.F, random.size.F = random.size.F,
              perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"),
              procD.lm.shape = PSh, procD.lm.size = PSz,
              call=match.call()) }
  if(object.sym==TRUE){
    out<-list(size.anova = NULL, shape.anova = shape.anova, symm.shape = symm.component,
              asymm.shape = asymm.component, DA.component = DA.mns, FA.component = FA.component,
              data.type = ifelse(object.sym==TRUE,"Object", "Matching"),
              FA.mns = FA.component, DA.mns = DA.mns,
              permutations = iter+1,
              random.shape.F = random.shape.F, random.size.F = NULL,
              perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"),
              procD.lm.shape = PSh,
              call=match.call()) }
  class(out) <- "bilat.symmetry"
  out  
}
