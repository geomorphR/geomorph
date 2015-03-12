#' Morphological disparity for one or more groups of specimens
#'
#' Function estimates morphological disparity and performs pairwise comparisons among groups. 
#'
#' The function estimates morphological disparity and performs pairwise comparisons to identify differences
#' between groups. The function takes as input GPA-aligned shape data [e.g., \code{\link{gpagen}}] and a grouping factor, and
#' estimates disparity as the Procrustes variance for each group, which is the sum of the diagonal elements 
#' of the group covariance matrix (e.g., Zelditch et al. 2012). The group Procrustes variances are used as 
#' test values, and these are then statistically evaluated through permutation, where the rows of the shape 
#' matrix are randomized relative to the grouping variable. The function can be used to obtain disparity for the whole
#' dataset by using a dummy group factor assigning all specimens to one group, in which case only Procrustes variance is returned.
#'
#' @param A A matrix (n x [p x k]) or 3D array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param groups A factor defining groups
#' @param iter Number of iterations for permutation test
#' @keywords analysis
#' @export
#' @author Emma Sherratt
#' @return Function returns a list with the following components: 
#'   \item{Disp.obs}{A matrix of Procrustes variance for each group}
#'   \item{Prob.Dist}{A matrix of pairwise significance levels based on permutation}
#' @references Zelditch, M. L., D. L. Swiderski, H. D. Sheets, and W. L. Fink. 2012. Geometric morphometrics 
#'   for biologists: a primer. 2nd edition. Elsevier/Academic Press, Amsterdam.
#' @examples
#' data(plethodon)
#' Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
#' morphol.disparity(Y.gpa$coords, groups=plethodon$site, iter = 99)
morphol.disparity <- function(A, groups, iter = 999){
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').") }
  if (length(dim(A))==3){ y<-two.d.array(A)}
  if (length(dim(A))==2){ y<-A }
  if (!is.factor(groups)) {
    groups <- factor(groups)
    cat("groups variable coerced into factor.\n")
  }
  m <- length(levels(groups))
  y = resid(lm(y ~ groups))
  procvar <- function(x){sum(dist(x)^2)/(nrow(x)^2)}
  if(m ==1){ d.obs <- procvar(y)
             return(d.obs)}
  if(m > 1){
    d.obs <- by(y, groups, procvar)
    diff.d.obs <- as.matrix(dist(d.obs))
    PDisp <- array(1, dim = c(m, m))
    for (i in 1:iter){
      y.r <- y[sample(nrow(y)),]
      d.rand <- by(y.r, groups, procvar)
      diff.d.rand <- as.matrix(dist(d.rand))
      PDisp <- ifelse(diff.d.rand >= diff.d.obs, PDisp + 1, PDisp)
    }
    PDisp <- PDisp/(iter + 1)
    d.obs <- as.matrix(d.obs)
    colnames(d.obs) <- "ProcVar"
    return(list(Disp.obs = d.obs, Prob.Disp = PDisp))
  }
}
