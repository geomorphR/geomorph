#' Analysis of bilateral symmetry
#'
#' Function performs an analysis of directional and fluctuating asymmetry for bilaterally symmetric objects 
#'
#' The function quantifies components of shape variation for a set of specimens as described by 
#' their patterns of symmetry and asymmetry. Here, shape variation is decomposed into variation 
#' among individuals, variation among sides (directional asymmetry), and variation due to an 
#' individual x side interaction (fluctuating symmetry). These components are then statistically 
#' evaluated using Procrustes ANOVA. Statistical assessment of model effects for shape variation 
#' is accomplished using permutation procedures. Methods for both matching symmetry and object 
#' symmetry can be implemented. Matching symmetry is when each object contains mirrored pairs of 
#' structures (e.g., right and left hands) while object symmetry is when a single object is 
#' symmetric about a midline (e.g., right and left sides of human faces). Details on general 
#' approaches for the study of symmetry in geometric morphometrics may be found in: Mardia et 
#' al. 2000; Klingenberg et al. 2002. 
#'  
#' As input, the function receives either A 3D array (p x k x n) containing raw landmarks (requiring 
#' GPA to be performed) or a gpagen object (if GPA has been previously performed) or a geomorphShapes 
#' object. If one wishes to incorporate semilandmarks, GPA can either be performed first using gpagen,
#' or within bilat.symmetry by passing adequate GPA arguments (i.e. curves, surfaces, ProcD etc, 
#' see \code{\link{gpagen}}. If a geomorphShapes object is provided, semilandmarks are automatically 
#' identified and slid during GPA. For "object.sym = FALSE, landmarks should be of dimension (p x k 
#' x 2n), as each specimen is represented by both left and right configurations.
#'    
#' Analyses of symmetry for matched pairs of objects is implemented when {object.sym=FALSE}. Here, 
#' a 3D array [p x k x 2n] contains the landmark coordinates for all pairs of structures (2 
#' structures for each of n specimens). Because the two sets of structures are on opposite sides,
#' they represent mirror images, and one set must be reflected prior to the analysis to allow 
#' landmark correspondence. IT IS ASSUMED THAT THE USER HAS DONE THIS PRIOR TO PERFORMING THE
#' SYMMETRY ANALYSIS. Reflecting a set of specimens may be accomplished by multiplying one coordinate 
#' dimension by '-1' for these structures (either the x-, the y-, or the z-dimension). A vector 
#' containing information on individuals and sides must also be supplied. Replicates of each 
#' specimen may also be included in the dataset, and when specified will be used as measurement 
#' error (see Klingenberg and McIntyre 1998). 
#' 
#' Analyses of object symmetry is implemented when {object.sym=TRUE}. Here, a 3D array [p x k x n] 
#' contains the landmark coordinates for all n specimens. To obtain information about asymmetry, 
#' the function generates a second set of objects by reflecting them about one of their coordinate 
#' axes. The landmarks across the line of symmetry are then relabeled to obtain landmark 
#' correspondence. The user must supply a list of landmark pairs. A vector containing information 
#' on individuals must also be supplied. Replicates of each specimen may also be included in the 
#' dataset, and when specified will be used as measurement error. 
#' 
#' The function also provides individual measures of signed and unsigned asymmetry, calculated as the
#' Procrustes distance between the right and left element (for paired structures, as detailed in 
#' Klingenberg and McIntyre 1998) or side of the structure (for object symmetry, following Lazić 
#' et al 2015). The computational difference between the two approaches consists in that, for object
#' symmetry, only paired landmarks are considered, excluding the landmarks of the midline.
#'  
#' \subsection{Notes for geomorph 3.0}{ 
#' Compared to older versions of geomorph, some results can be expected to be slightly different. 
#' Starting with geomorph 3.0, results use only type I sums of squares (SS) with either full 
#' randomization of raw shape values or RRPP (preferred with nested terms) for analysis of variance
#' (ANOVA).  Older versions used a combination of parametric and non-parametric results, as well as 
#' a combination of type I and type III SS.  While analytical conclusions should be consistent 
#' (i.e., "significance" of effects is the same), these updates maintain consistency in analytical 
#' philosophy.  This change will require longer computation time for large datasets, but the 
#' trade-off allows users to have more flexibility and eliminates combining disparate analytical 
#' philosophies. 
#'  
#' Note also that significance of terms in the model are found by comparing F-values for each term 
#' to those obtained via permutation.  F-ratios and df are not strictly necessary (a ratio of SS 
#' would suffice), but they are reported as is standard for anova tables. Additionally, users will 
#' notice that the df reported are based on the number of observations rather than a combination 
#' of objects * coordinates * dimensions, as is sometimes found in morphometric studies of symmetry. 
#' However, this change has no effect on hypothesis testing, as only SS vary among permutations (df, 
#' coordinates, and dimensions are constants). 
#' }
#'  
#' The generic functions, \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} all work 
#' with \code{\link{bilat.symmetry}}.
#'
#' @param A One of either A 3D array (p x k x n) containing raw landmarks (requiring GPA to be 
#' performed) or a gpagen object (if GPA has been previously performed) or a geomorphShapes object 
#' (requiring GPA to be performed).  Any gpagen argument should work within bilat.symmetry.
#' @param ind A vector containing labels for each individual. For matching symmetry, the matched 
#' pairs receive the same label (replicates also receive the same label).
#' @param side An optional vector (for matching symmetry) designating which object belongs to which
#' 'side-group'
#' @param replicate An optional vector designating which objects belong to which group of replicates.
#' Alternatively, this can be a character value to indicate the name of the variable in the data frame to use.
#' @param object.sym A logical value specifying whether the analysis should proceed based on object 
#' symmetry {=TRUE} or matching symmetry {=FALSE}
#' @param land.pairs An optional matrix (for object symmetry) containing numbers for matched pairs 
#' of landmarks across the line of symmetry 
#' @param data A data frame for the function environment, see \code{\link{geomorph.data.frame}}. It 
#' is imperative that the variables "ind", "side", and "replicate" in the data frame match these 
#' names exactly (as shown in examples below).  
#' @param iter Number of iterations for significance testing.
#' @param seed An optional argument for setting the seed for random permutations of the resampling 
#' procedure. If left NULL (the default), the exact same P-values will be found for repeated runs 
#' of the analysis (with the same number of iterations). If seed = "random", a random seed will be 
#' used, and P-values will vary.  One can also specify an integer for specific seed values,
#' which might be of interest for advanced users.
#' @param RRPP A logical value indicating whether residual randomization should be used for 
#' significance testing.
#' @param SS.type A choice between type I (sequential), type II (hierarchical), or type III (marginal).
#' @param turbo A logical value that if TRUE, suppresses coefficient estimation in every random permutation, 
#' in order to speed up computation time.
#' @param Parallel Either a logical value to indicate whether parallel processing 
#' should be used or a numeric value to indicate the number of cores to use in 
#' parallel processing via the \code{parallel} library. 
#' If TRUE, this argument invokes forking of all processor cores, except one.  If
#' FALSE, only one core is used. A numeric value directs the number of cores to use,
#' but one core will always be spared.
#' @param print.progress A logical value to indicate whether a progress bar should be printed to 
#' the screen.  
#' This is helpful for long-running analyses.
#' @param ... Arguments to pass onto gpagen
#' @keywords analysis
#' @export
#' @author Dean Adams, Michael Collyer and Antigoni Kaliontzopoulou
#' @return An object of class "bilat.symmetry" returns a list of the following
#' \item{shape.anova}{An analysis of variance table for the shape data.}
#' \item{size.anova}{An analysis of variance table for the shape data (when object.sym = FALSE).}
#' \item{symm.shape}{The symmetric component of shape variation.}
#' \item{asymm.shape}{The asymmetric component of shape variation.}
#' \item{DA.component}{The directional asymmetry component, found as the mean shape for each side.}
#' \item{FA.component}{The fluctuating asymmetry component for each specimen, 
#' found as the specimen specific side deviation adjusted for the mean 
#' directional asymmetry in the dataset.}
#' \item{signed.AI}{Individual signed asymmetry index, as per Klingenberg and McIntyre, 1998;
#' Lazić et al 2015.}
#' #' \item{unsigned.AI}{Individual unsigned asymmetry index, as per Klingenberg and McIntyre, 1998;
#' Lazić et al 2015.}
#' \item{data.type}{A value indicating whether the analysis was performed as Object or Matching 
#' symmetry.}
#' \item{permutations}{The number of random permutations used.}
#' \item{random.shape.F}{A matrix of random F-values from the Shape analysis.}
#' \item{random.size.F}{A matrix of random F-values from the Centroid Size analysis (when 
#' object.sym = FALSE).}
#' \item{perm.method}{A value indicating whether "Raw" values were shuffled or "RRPP" performed.}
#' \item{procD.lm.shape}{A list of typical output from an object of class procD.lm, for shape}
#' \item{procD.lm.size}{If applicable, a list of typical output from an object of class procD.lm, 
#' for size (when object.sym = FALSE).}
#' \item{call}{The matched call.}
#' 
#' @references Klingenberg, C.P. and G.S. McIntyre. 1998. Quantitative genetics of geometric shape 
#' in the mouse mandible. Evolution. 55:2342-2352.
#' @references Mardia, K.V., F.L. Bookstein, and I.J. Moreton. 2000. Statistical assessment of 
#' bilateral symmetry of shapes. Biometrika. 87:285-300.
#' @references Klingenberg, C.P., M. Barluenga, and A. Meyer. 2002. Shape analysis of symmetric 
#' structures: quantifying variation among individuals and asymmetry. Evolution. 56:1909-1920.
#' @references Lazić, M. M., M. A. Carretero, J. Crnobrnja-Isailović, and A. Kaliontzopoulou. 2015.
#' Effects of environmental disturbance on phenotypic variation: an integrated assessment of 
#' canalization, developmental stability, modularity, and allometry in lizard head shape. The American
#' Naturalist 185:44–58.
#' @examples
#' #Example of matching symmetry
#' # NOT RUN
#' # data(mosquito)
#' # gdf <- geomorph.data.frame(wingshape = mosquito$wingshape, 
#' # ind=mosquito$ind, 
#' # side=mosquito$side,
#' # replicate=mosquito$replicate)
#' # mosquito.sym <- bilat.symmetry(A = wingshape, ind = ind, side = side,
#' # replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 149, 
#' # data = gdf)
#' # summary(mosquito.sym)
#' # plot(mosquito.sym, warpgrids = TRUE)
#' # mosquito.sym$shape.anova # extract just the anova table on shape
#' 
#' # Previous example, performing GPA first
#' # Y.gpa <- gpagen(mosquito$wingshape)
#' # mosquito.sym2 <- bilat.symmetry(A = Y.gpa, ind = ind, side = side,
#' # replicate = replicate, object.sym = FALSE, RRPP = TRUE, iter = 149, 
#' # data = gdf)
#' # summary(mosquito.sym2)
#' # summary(mosquito.sym) # same results
#'
#' #Example of object symmetry
#'
#' # data(lizards)
#' # gdf <- geomorph.data.frame(shape = lizards$coords, 
#' # ind = lizards$ind, 
#' # replicate = lizards$rep)
#' # liz.sym <- bilat.symmetry(A = shape, ind = ind, rep = rep, 
#' # object.sym = TRUE, 
#' # land.pairs = lizards$lm.pairs, data = gdf, RRPP = TRUE, iter = 149)
#' # summary(liz.sym)
#' 
#' # Example of object symmetry in 3D and including semilandmarks
#' 
#' # data(scallops)
#' # gdf <- geomorph.data.frame(shape = scallops$coorddata, 
#' # ind = scallops$ind)
#' # scallop.sym <- bilat.symmetry(A = shape, ind = ind, 
#' # object.sym = TRUE, 
#' # curves= scallops$curvslide, surfaces = scallops$surfslide,
#' # land.pairs=scallops$land.pairs, data = gdf, RRPP = TRUE, iter = 149)
#' # summary(scallop.sym)
#' # NOTE one can also: plot(scallop.sym, warpgrids = TRUE, mesh = NULL)
#' # NOTE one can also: scallop.sym$data.type # recall the symmetry type

bilat.symmetry <- function(A, ind = NULL, side = NULL, replicate = NULL, object.sym = FALSE, land.pairs = NULL,
                           data = NULL, iter = 999, seed = NULL, RRPP = TRUE, SS.type = c("I", "II", "III"),
                           turbo = TRUE, Parallel = FALSE, print.progress = TRUE, ...){
  
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
  
  n <- dim(A)[3]; k <- dim(A)[2]; p <- dim(A)[1]; nind <- nlevels(ind); spec.names <- dimnames(A)[[3]]
  if(!object.sym && is.null(side)) stop("Sides not specified.") 
  
  if(object.sym){
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
  PSh <- procD.lm(form.shape, data = dat.shape, RRPP = RRPP, SS.type = SS.type,
                  Parallel = Parallel, turbo = turbo,
                  seed = seed, iter = iter, print.progress = print.progress, 
                  effect.type = "F")
  shape.anova <- anova(PSh, print.progress = FALSE, effect.type = "F")$table
  shape.anova$Z[is.nan(shape.anova$Z)] <- NA
  
  MS <- PSh$ANOVA$MS
  RSS <- PSh$ANOVA$RSS
  MSE <- RSS / matrix(PSh$ANOVA$df[nrow(RSS) + 1], nrow(RSS), ncol(RSS))
  random.shape.F <- MS/MSE
  
  if(length(form.names) > 4) {
    MS <- random.shape.F <- PSh$ANOVA$MS
    MS.mod <- PSh$ANOVA$RSS.model[3,]/PSh$ANOVA$df[4]
    random.shape.F[1,] <- MS[1,]/MS[3,]
    random.shape.F[2,] <- MS[2,]/MS[3,]
    random.shape.F[3,] <- PSh$ANOVA$Fs[3,]
    
    newZ <- apply(random.shape.F, 1, effect.size)
    newP <- apply(random.shape.F, 1, pval)
    shape.anova$F[1:3] <- random.shape.F[1:3, 1]
    shape.anova$Z[1:3] <- newZ
    shape.anova[[ncol(shape.anova)]][1:3] <- newP
  }
  rownames(shape.anova) <- form.names
  
  
  if(!object.sym){  
    if(!is.null(replicate)) {
      form.size <- size ~ ind + side + ind/side 
      dat.size <- geomorph.data.frame(size = size, ind = ind, side = side, replicate = replicate)
    } else {
      form.size <- size ~ ind + side 
      dat.size <- geomorph.data.frame(size = size, ind = ind, side = side)
    }
    if(print.progress) cat("\nSize Analysis")
    PSz <- procD.lm(form.size, data = dat.size, RRPP = RRPP, 
                    SS.type = SS.type,
                    Parallel = Parallel, turbo = turbo,
                    seed = seed, iter= iter, print.progress = print.progress, 
                    effect.type = "F")
    size.anova <- anova(PSz, print.progress = FALSE, effect.type = "F")$table
    size.anova$Z[is.nan(size.anova$Z)] <- NA
    
    MS <- PSz$ANOVA$MS
    RSS <- PSz$ANOVA$RSS
    MSE <- RSS / matrix(PSz$ANOVA$df[nrow(RSS) + 1], nrow(RSS), ncol(RSS))
    random.size.F <- MS/MSE
    
    if(length(form.names) > 4) {
      MS <- random.size.F <- PSz$ANOVA$MS
      MS.mod <- PSz$ANOVA$RSS.model[3,]/PSz$ANOVA$df[4]
      random.size.F[1,] <- MS[1,]/MS[3,]
      random.size.F[2,] <- MS[2,]/MS[3,]
      random.size.F[3,] <- PSz$ANOVA$Fs[3,]
      newZ <- apply(random.size.F, 1, effect.size)
      newP <- apply(random.size.F , 1, pval)
      size.anova$F[1:3] <- random.size.F[1:3, 1]
      size.anova$Z[1:3] <- newZ
      size.anova[[ncol(size.anova)]][1:3] <- newP
    }
    rownames(PSz$aov.table) <- form.names
    size.anova <- PSz$aov.table
  }
  
  # build shape components for output
  X.ind <- model.matrix(~ind + 0, data = as.data.frame(dat.shape[-1]))
  symm.component <- arrayspecs(coef(lm.fit(X.ind, Y)),p,k)
  ind.names <- substr(dimnames(symm.component)[[3]], start=4,
                      stop=nchar(dimnames(symm.component)[[3]]))
  dimnames(symm.component)[[3]] <- ind.names
  X.side <- model.matrix(~(side:ind) + 0, data = as.data.frame(dat.shape[-1]))
  avg.side.symm <- coef(lm.fit(X.side, Y))
  n.ind <- nlevels(ind)
  n.side <- nlevels(side)
  indsq <- seq(n.side, (n.ind*n.side), n.side)
  asymm.component <- avg.side.symm[indsq,] - avg.side.symm[-indsq,]
  mn.shape <- mshape(A)
  asymm.component <- simplify2array(lapply(1:n.ind, function(j) {
    t(matrix(asymm.component[j,],k,p)) + mn.shape
    }))
  dimnames(asymm.component)[[3]] <- ind.names
  
  DA.est <- coef(.lm.fit(X.side, Y))
  DA.mns <- arrayspecs(rbind(apply(DA.est[-indsq,], 2, mean), apply(DA.est[indsq,], 2, mean)), p, k)
  mn.DA <- matrix(apply((DA.est[-indsq,] - DA.est[indsq,]), 2, mean), byrow=T, nrow=p, ncol=k)
  
  X.ind.side <- model.matrix(~(side:ind) + 0, data = as.data.frame(dat.shape[-1]))
  ind.mns <- coef(.lm.fit(X.ind.side, Y))
  FA.component <- ind.mns[-indsq,] - ind.mns[indsq,]
  FA.component <- simplify2array(lapply(1:n.ind, function(j) 
  {t(matrix(FA.component[j,],k,p)) + mn.shape - mn.DA}))
  dimnames(FA.component)[[3]] <- ind.names 
  
  # Calculate individual asymmetry indices
  signed.asymm <- two.d.array(asymm.component)
  signed.AI <- sqrt(apply(signed.asymm^2, 1, sum))
  names(signed.AI) <- ind.names
  
  asymm.mean <- apply(signed.asymm, 2, mean)
  unsigned.asymm <- matrix(NA, nrow=nrow(signed.asymm), ncol=ncol(signed.asymm))
  for (i in 1:nrow(signed.asymm)){
    unsigned.asymm[i,] <- ifelse(signed.asymm[i,]%*%asymm.mean > 0, signed.asymm[i,], signed.asymm[i,]*(-1))
  } 
  unsigned.AI <- sqrt(apply(unsigned.asymm^2, 1, sum))
  names(unsigned.AI) <- ind.names

  out <- list(shape.anova = shape.anova, symm.shape = symm.component,
              asymm.shape = asymm.component, DA.component = DA.mns, FA.component = FA.component,
              signed.AI = signed.AI, unsigned.AI = unsigned.AI,
              data.type = ifelse(object.sym == TRUE, "Object", "Matching"),
              permutations = iter+1,
              random.shape.F = random.shape.F, 
              perm.method = ifelse(RRPP==TRUE,"RRPP", "Raw"),
              procD.lm.shape = PSh)
  if(!object.sym) {
    out$size.anova = size.anova
    out$procD.lm.size = PSz
    out$random.size.F = random.size.F
  }
  out$call <- match.call()
  class(out) <- "bilat.symmetry"
  out  
}
