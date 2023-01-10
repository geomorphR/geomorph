#' Comparisons of Phylogenetic Signal Effect Sizes 
#'
#' Function performs an analysis to compare the effect sizes of two or more phylogenetic effect sizes
#'
#' The function statistically compares the effect sizes of two or more \code{\link{physignal.z}} analyses.  
#' This can be performed on different traits from the same tree, same or different traits from different
#' trees, or modules of landmark configurations.
#'  
#' To use this function, perform \code{\link{physignal.z}} on as many samples as desired.  
#' Any number of objects of class physignal.z can be input.  Note that some values of Z can be NaN,
#' if the scaling parameter, lambda, is optimized at 0.  In these cases, the standard error is also 0,
#' and pairwise comparisons might not make sense.
#'  
#' 
#' @param ... saved analyses of class physignal.z
#' @param two.tailed A logical value to indicate whether a two-tailed test (typical and default) should be performed.
#' @keywords analysis
#' @export
#' @author Michael Collyer
#' @return An object of class compare.physignal.z, returns a list of the following
#' \item{sample.z}{A vector of effect sizes for each sample.}
#' \item{sample.r.sd}{A vector of standard deviations for each sampling distribution (following Box-Cox transformation).}
#' \item{pairwise.z}{A matrix of pairwise, two-sample z scores between all pairs of effect sizes.}
#' \item{pairwise.p}{A matrix of corresponding P-values.}
#' @references Collyer,  M.L., E.K. Baken, & D.C. Adams.  2022. A standardized effect size for evaluating
#' and comparing the strength of phylogenetic signal. Methods in Ecology and Evolution. 13:367-382.
#' @examples
#'
#' # Example: Compare phylogenetic signal of head components in Plethodon
#' 
#' data(plethspecies) 
#' Y.gpa<-gpagen(plethspecies$land)    #GPA-alignment
#' 
#' ## landmarks of the jaw and cranium
#' jaw <- 1:5
#' cranium <- 6:11
#' 
#' PS.jaw <- physignal.z(A = Y.gpa$coords[jaw,,], phy = plethspecies$phy, 
#' lambda = "front", PAC.no = 7, iter=999)
#' 
#' PS.cranium <- physignal.z(A = Y.gpa$coords[cranium,,], phy = plethspecies$phy, 
#' lambda = "front", PAC.no = 7, iter=999)
#' 
#' PS.list <-list(PS.jaw, PS.cranium)
#' names(PS.list) <- c("jaw", "cranium")
#' 
#' PS.Z <- compare.physignal.z(PS.list)
#' summary(PS.Z)
#' 
compare.physignal.z <- function(..., two.tailed = TRUE){
   dots <- list(...)
   tails <- if(two.tailed) 2 else 1
   if(length(dots) == 1) n <- length(dots[[1]]) else n <- length(dots)
   if(n == 1) stop("At least two objects of class physignal.z are needed")
   if(length(dots) == 1) {
     list.names <- names(dots[[1]]) 
     dots <- lapply(1:n, function(j) dots[[1]][[j]])
     names(dots) <- list.names
     } else list.names <- names(dots)
   if(length(dots) < 2) stop("At least two objects of class pls are needed")
   
   is.psz <- function(x) inherits(x, "physignal.z")
   
   sdn <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
   
   list.check <- sapply(1:length(dots), function(j) any(is.psz(dots[[j]])))
   
   if(any(!list.check)) stop("Not all objects are class pls")
   
   k <- length(list.check)
   if(is.null(list.names)) list.names <- as.list(substitute(list(...)))[-1L]
   k.combn <- combn(k,2)
   
   defW <- getOption("warn") 
   options(warn = -1) 
   bct <- lapply(dots, function(x) {
     s <- sd(x$rand.logL)
     res <- if(s > 0) box.cox(x$rand.logL)$transformed else 
       rep(0, length(x$rand.logL))
     })
     
   list.drs <- sapply(1:k, function(j) bct[[j]][1] - mean(bct[[j]])) 
   list.sds <- sapply(1:k, function(j) sdn(bct[[j]]))
   list.zs <- sapply(1:k, function(j) effect.size(dots[[j]]$rand.logL, center=TRUE))
   options(warn = defW)

   z12 <- sapply(1:ncol(k.combn), function(j){
     a <- k.combn[1,j]; b <- k.combn[2,j]
     r1 <- list.drs[a]; r2 <- list.drs[b] 
     sd1 <- list.sds[a]; sd2 <- list.sds[b]
     denom <- sqrt(sd1^2+sd2^2)
     res <- if (denom > 0) (r1-r2)/denom else 0
   })
   
   z12.p <- sapply(1:length(z12), function(j) 
     pnorm(abs(z12[[j]]), lower.tail = FALSE) * tails)
   
   d <- rep(0,k); names(d) <- list.names
   D <-dist(d)
   z12.pw <- p12.pw <- D
   for(i in 1:length(z12)) z12.pw[i] <-z12[i]
   for(i in 1:length(z12)) p12.pw[i] <-z12.p[i]
   names(list.zs) <- names(list.sds) <-list.names
   pairwise.z <- as.matrix(z12.pw)
   pairwise.P <- as.matrix(p12.pw)
   diag(pairwise.P) <- 1
   
   out <- list(sample.z = list.zs,
               sample.logL.sd = list.sds,
               pairwise.z = abs(pairwise.z),
               pairwise.P = pairwise.P)
   
   class(out) <- "compare.physignal.z"
   out
 }
 