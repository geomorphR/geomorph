#' Update procD.lm objects with nested effects
#' 
#' This function is used to update \code{\link{procD.lm}} objects for fixed effects with 
#' nested random effects (nested design)
#'
#' This functions serves as a helper function when linear models have nested (hierarchical) structure.  It is used on
#' \code{\link{procD.lm}} objects that were formerly evaluated with type I sums of squares (SS), as is typical
#' with models with only fixed effects.  Using a formula for nested effects, this function identifies the fixed and random
#' SS in the random outcomes used to generate the \code{\link{procD.lm}} object, and updates the F-values,
#' Z-scores, and P-values based on F values adjusted to be MS fixed/MS random (nested).  This is accomplished by 
#' generating random values for each iteration previously used in the \code{\link{procD.lm}} object.
#' 
#' This function can be used recursively for multiple updates, when multiple nested effects are used.  The function can
#' currently only handle single factors nested within other single factors.  
#' 
#' Function returns the same list as \code{\link{procD.lm}} but with new random F values and Cohen's f-squared values
#' substituted  The ANOVA table is updated in terms of F-values, Z-scores, and P-values.  Z-scores are re-calculated 
#' for all effects to be consistent with the type of distribution used.  If either Cohen's f-squared values or F values 
#' were originally chosen, the same distributions are used in the update; if SS values were originally chosen, the 
#' distribution is changed to Cohen's f-squared to calculate Z-scores. This change assures consistency in effect size 
#' estimation, as the effect that is updated cannot have an effect size based on SS.
#' 
#' It is important that the formula is input correctly. It can be input as one of the following four styles:
#' 
#'  ~ fixed/random
#' 
#'  ~ fixed + fixed/random
#' 
#' The two formulae above achieve the same model terms for the expanded model: ~ fixed + fixed:random
#' 
#'  ~ random + fixed/random
#' 
#'  ~ random + fixed + fixed/random
#' 
#' The two formulae above achieve the same model terms for the expanded model: ~ random + fixed + random:fixed
#' 
#' The \code{\link{procD.lm}} object will be updated in the same way with either of the approaches.  First, the F-value
#' for the fixed term will be adjusted as MS-fixed/MS-interaction for every random permutation.  Second, the P-value for the fixed
#' effect will be estimated from this new distribution of F-values.  Although the function will try to catch improper 
#' formulae and alert the user, it is possible the function will work with an improper formula.  Thus, adherence to one of the 
#' formulae above is recommended for best results.
#' 
#' Effect sizes (Z scores) are based on either the distribution of random F values or a distribution of
#' Cohen's f-squared values, calculated in every permutation.  An attempt will be made to preserve the effect size
#' type used in the previous \code{\link{procD.lm}} or \code{\link{procD.pgls}} analysis.  However, an analysis
#' performed in \code{\link{procD.lm}} using effect size calculated from random SS values will be updated
#' to use random Cohen's f-squared values for all effects, to avoid having effect sizes measured from different
#' distributions in the same analysis.
#' 
#' @param P A \code{\link{procD.lm}} object
#' @param f1 A right-hand or full formula for one factor nested within another; e.g., ~ A/B or ~B + A/B
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{procD.lm}}
#' @keywords utilities
#' @references Collyer, M.L., D.J. Sekora, and D.C. Adams. 2015. A method for analysis of phenotypic change for phenotypes described 
#' by high-dimensional data. Heredity. 115:357-365.
#' @examples
#' 
#'data(larvalMorph)
#'Y.gpa <- gpagen(larvalMorph$tailcoords)
#'gdf <- geomorph.data.frame(Y.gpa, Treatment = larvalMorph$treatment, Family = larvalMorph$family)
#'
#'# Model with fixed and nested effects
#'tailANOVA <- procD.lm(coords ~ Treatment/Family, iter = 99, RRPP=TRUE, data=gdf)
#'summary(tailANOVA)
#'
#'# Update for nested effects
#'tailANOVA <- nested.update(tailANOVA, ~ Treatment/Family)
#'summary(tailANOVA)
#'
#'# Model with random, fixed, and nested effects
#'tailANOVA <- procD.lm(coords ~ Family + Treatment/Family, iter=99, RRPP=TRUE, data=gdf)
#'summary(tailANOVA)
#'
#'# Update for nested effects
#'tailANOVA <- nested.update(tailANOVA, ~ Family + Treatment/Family)
#'summary(tailANOVA)
#'
#' # One needs to be careful using this function!
#' 
#'tailANOVA <- procD.lm(coords ~ Csize * Treatment/Family, iter=99, RRPP=TRUE, data=gdf)
#'
#' # This will not work: tailANOVA <- nested.update(tailANOVA, ~ Treatment/Family) 
#' # The updated terms must be included as part of the original terms
#' 
#'tailANOVA <- procD.lm(coords ~ Csize + Treatment/Family, iter=99, RRPP=TRUE, data=gdf)
#'summary(tailANOVA)
#'
#'# Now the format will allow an update
#'
#'tailANOVA <- nested.update(tailANOVA, ~ Treatment/Family) 
#'
#'summary(tailANOVA)
nested.update <- function(P, f1){
  Terms <- terms(f1)
  term.lbls <- attr(Terms, "term.labels")
  AT <- P$aov.table
  rFs <- P$random.F
  rSS <- P$random.SS
  rcohenf <- P$cohenf
  effect.type <- P$effect.type
  SS.type <- P$SS.type
  if(P$effect.type == "SS") P$effect.type = "cohenf"
  k <- length(term.lbls)
  n <- NROW(P$Y)
  if(k < 2  | k > 3) stop("Formula is not correct.  Should be of form ~ A/B or ~ B + A/B; where A is fixed and B is random")
  if(k == 2) {
    if(length(strsplit(term.lbls[2], ":")[[1]]) != 2) 
      stop("Formula is not correct.  Should be of form ~ A/B or ~ A + B + A/B; where A is fixed and B is random")
    fixedrow = which(rownames(AT) == term.lbls[1])
    randomrow = which(rownames(AT) == term.lbls[2])
  }
  if(k == 3) {
    if(length(strsplit(term.lbls[3], ":")[[1]]) != 2) 
      stop("Formula is not correct.  Should be of form ~ A/B or ~ A + B + A/B; where A is fixed and B is random")
    fixedrow = which(rownames(AT) == term.lbls[2])
    randomrow = which(rownames(AT) == term.lbls[3])
  } 
  
  if(length(fixedrow) == 0 | length(randomrow) == 0)
    stop("Formula provided is not consitent with formula used in procD.lm object")
  MSF <- AT[fixedrow,3]
  MSR <- AT[randomrow,3]
  AT[fixedrow, 5] <- MSF/MSR
  SSf <- P$random.SS[fixedrow,]
  SSr <- P$random.SS[randomrow,]
  dff <- AT[fixedrow, 1]
  dfr <- AT[randomrow,1]
  Fs <- (SSf/dff)/(SSr/dfr)
  rFs[fixedrow,] <- Fs
  if(!is.null(P$PGLS)){
    SSE <- P$random.SSE
    SSY <- colSums(rSS) + SSE 
    SSE.mat <- matrix(SSE, NROW(rSS), NCOL(rSS), byrow = TRUE)
  } else {
    SSY <- sum(center(P$Y)^2)
    SSE <- SSY - colSums(rSS)
    SSE.mat <- matrix(SSE, NROW(rSS), NCOL(rSS), byrow = TRUE)
  }
  if(SS.type == "III") {
    etas <- rSS/(rSS+SSE.mat)
    rcohenf <- etas/(1-etas)
  } else {
    etas <- rSS/SSY
    unexp <- 1 - apply(etas, 2, cumsum)
    rcohenf <- etas/unexp
  }
  newP <- pval(Fs)
  if(effect.type == "F") newZ <- apply(log(rFs), 1, effect.size) else
    newZ <- apply(log(rcohenf), 1, effect.size)
  AT[fixedrow,7] <- newP
  AT[,6] <- c(newZ, NA, NA)
  P$aov.table <- AT
  P$random.F <- rFs
  P$random.cohenf <- rcohenf
  P$nested.update <- TRUE
  return(P)
}