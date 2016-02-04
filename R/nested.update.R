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
#' Function returns the same list as \code{\link{procD.lm}} but with random F values appended.  The 
#' ANOVA table is updated in terms of F-values, Z-scores, and P-values.  It is important that the formula is input correctly.  
#' It can be input as one of the following six styles:
#' 
#' ~ fixed/random
#' ~ fixed + fixed/random
#' ~ fixed + random%in%fixed
#' 
#' The three formulae above achieve the same model terms for the expanded model: ~ fixed + fixed:random
#' 
#' ~ random + fixed/random
#' ~ random + fixed + fixed/random
#' ~ random + fixed + random%in%fixed
#' 
#' The three formulae above achieve the same model terms for the expanded model: ~ random + fixed + random:fixed
#' 
#' The \code{\link{procD.lm}} object will be updated in the same way with either of the approached.  First, the F-value
#' for the fixed term will be adjusted as MS-fixed/MS-interaction for every random permutation.  Second, the P-value for the fixed
#' effect will be estimated from this new distribution of F-values.  Although the function will try to catch improper 
#' formulae and alert the user, it is possible the function will work with an improper formula.  Thus, adherence to one of the 
#' formulae above is recommended for best results.
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
#'data(larvalTails)
#'Y.gpa <- gpagen(larvalTails$landmarks)
#'gdf <- geomorph.data.frame(Y.gpa, Treatment = larvalTails$Treatment, Family = larvalTails$Family)
#'
#'# Model with fixed and nested effects
#'tailANOVA <- procD.lm(coords ~ Treatment/Family, iter=99, RRPP=TRUE, data=gdf)
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
  k <- length(term.lbls)
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
  newP <- pval(Fs)
  newZ <- effect.size(Fs)
  AT[fixedrow,7] <- newP
  AT[fixedrow,6] <- newZ
  P$aov.table <- AT
if(is.null(P$randomF))
     {Fs <- matrix(Fs,1,length(Fs)); rownames(Fs) <- term.lbls[1]} else
     {
       newnames <- c(rownames(P$randomF), term.lbls[1])
       Fs <- rbind(P$randomF, Fs)
       rownames(Fs) <- newnames
     }
  P$random.F <- Fs
  return(P)
}