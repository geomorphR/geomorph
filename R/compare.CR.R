#' Comparisons of Effect Sizes from Modularity Analyses
#'
#' Function performs an analysis to compare the effect sizes of two or more CR effects
#'
#' The function statistically compares the effect sizes of two or more CR analyses.  Typically, this
#' function might be used to compare levels of modularity between two or more samples, each measuring the degree 
#' of morphological modularity in each.  Alternatively, the approach can compare the degree of modular signal as 
#' expressed by alternative modular hypotheses for the same dataset. 
#' 
#' The analysis calculates effect sizes as standard deviates, z, and performs two-sample z-tests, using the pooled 
#' standard error from the sampling distributions of the CR analyses. The method follows that of Adams and Collyer (2019) used 
#' to compare patterns of modularity across datasets.
#'  
#' To use this function, simply perform \code{\link{modularity.test}}, or \code{\link{phylo.modularity}} on as many samples or 
#' alternative modular hypotheses as desired.  Any number of objects of class CR can be input. For the case of the latter, one may wish to
#' include the null hypothesis of no modularity (i.e., that all variables belong to a single module). For this, the CR.null = TRUE option
#' should be specified. Finally, one may perform the comparison as either a one-tailed or a two-tailed (default) test.
#' 
#' @param ... saved analyses of class CR
#' @param CR.null A logical value to indicate whether a Null CR model (no modularity) should also be included in analysis.
#' @param two.tailed A logical value to indicate whether a two-tailed test (typical and default) should be performed.
#' @keywords analysis
#' @export
#' @author Dean Adams and Michael Collyer
#' @return An object of class compare.CR, returns a list of the following
#' \item{sample.z}{A vector of effect sizes for each sample.}
#' \item{sample.r.sd}{A vector of standard deviations for each sampling distribution.}
#' \item{pairwise.z}{A matrix of pairwise, two-sample z scores between all pairs of effect sizes.}
#' \item{pairwise.p}{A matrix of corresponding P-values.}
#' @references Adams, D.C. and M.L. Collyer. 2019.  Comparing the strength of modular signal, and evaluating alternative modular hypotheses,
#' using covariance ratio effect sizes for morphometric data. *Evolution*. (Accepted).
#' @examples
#' 
#' ### NOT RUN
#' # Example 1: Compare modular signal across datasets
#'
#' # data(pupfish) 
#' # Y.gpa<-gpagen(pupfish$coords, print.progress = FALSE)    #GPA-alignment    
#' ## landmarks on the body and operculum
#' # land.gps<-rep('a',56); land.gps[39:48]<-'b'
#' 
#' # group <- factor(paste(pupfish$Pop, pupfish$Sex, sep = "."))
#' # levels(group)
#' 
#' # coords.gp <- coords.subset(Y.gpa$coords, group)
#' 
#' # modularity.tests <- lapply(1:nlevels(group), function(j) modularity.test(coords.gp[[j]],
#' #  land.gps, iter = 499, print.progress = FALSE))
#' ## the lapply function performs the modularity test on each 3D array in the lists provided
#' 
#' # modularity.tests[[1]]
#' # modularity.tests[[2]]
#' # modularity.tests[[3]]
#' # modularity.tests[[4]]
#' 
#' # group.Z <- compare.CR(modularity.tests, CR.null = FALSE)
#' # group.Z  ## NOTE: need a summary function 
#' 
#' # Example 2: Compare alternative modular hypotheses
#' 
#' # land.gps3 <- rep('a',56); land.gps3[39:48]<-'b'; land.gps3[c(6:9,28:38)] <- 'c' 
#'    #3 module hypothesis (tail now a module)
#' # land.gps4 <- rep('a',56); land.gps4[39:48]<-'b'; land.gps4[c(6:9,28:38)] <- 'c'; 
#'  # land.gps4[c(10,49:56)] <- 'd'  #4 module hypothesis (eye now a module)
#' 
#' # m3.test <- modularity.test(coords.gp$Marsh.F,land.gps3, iter = 499, print.progress = FALSE)
#' # m4.test <- modularity.test(coords.gp$Marsh.F,land.gps4, iter = 499, print.progress = FALSE)
#' 
#' # model.Z <- compare.CR(modularity.tests[[1]],m3.test,m4.test, CR.null = TRUE)
#' # model.Z 
#' 
compare.CR <- function(...,CR.null = TRUE, two.tailed = TRUE){
   dots <- list(...)
   tails <- if(two.tailed) 2 else 1
   if(length(dots) == 1) n <- length(dots[[1]]) else n <- length(dots)
   if(n == 1) stop("At least two objects of class CR are needed")
   if(length(dots) == 1) {
     list.names <- names(dots[[1]]) 
     dots <- lapply(1:n, function(j) dots[[1]][[j]])
     names(dots) <- list.names
     } else list.names <- names(dots)
   
   if(length(dots) < 2) stop("At least two objects of class CR are needed")
     
   if(is.null(list.names)) list.names <- paste("CR", 1:n, sep = ".")
   names(dots) <- list.names
        
   is.CR <- function(x) inherits(x, "CR")
   sdn <- function(x) sqrt(sum((x-mean(x))^2)/length(x))
   list.check <- sapply(1:length(dots), function(j) any(is.CR(dots[[j]])))
   if(any(list.check == FALSE)) stop("Not all objects are class CR")
   k <- length(list.check)
   if(is.null(list.names)) list.names <- as.list(substitute(list(...)))[-1L]
   k.combn <- combn(k,2)
   list.drs <- sapply(1:k, function(j) dots[[j]]$random.CR[1] - mean(dots[[j]]$random.CR[-1])) 
   list.sds <- sapply(1:k, function(j) sdn(dots[[j]]$random.CR[-1]))
   list.zs <- sapply(1:k, function(j) effect.size(dots[[j]]$random.CR, center=TRUE))  
   if (CR.null == TRUE){
      k <- k + 1
      k.combn <- combn(k,2)
      list.drs <- c(0,list.drs) 
      list.sds <- c(0,list.sds)
      list.zs <- c(0,list.zs)
   }
   z12 <- sapply(1:ncol(k.combn), function(j){
     a <- k.combn[1,j]; b <- k.combn[2,j]
     r1 <- list.drs[a]; r2 <- list.drs[b] 
     sd1 <- list.sds[a]; sd2 <- list.sds[b]
     (r1-r2)/sqrt(sd1^2+sd2^2)
   })
   z12.p <- sapply(1:length(z12), function(j) pnorm(abs(z12[[j]]), lower.tail = FALSE)* tails)  
   d <- rep(0,k); names(d) <- list.names
   D <-dist(d)
   z12.pw <- p12.pw <- D
   for(i in 1:length(z12)) z12.pw[i] <-z12[i]
   for(i in 1:length(z12)) p12.pw[i] <-z12.p[i]
   names(list.zs) <- names(list.sds) <-list.names
   pairwise.z <- as.matrix(z12.pw)
   pairwise.P <- as.matrix(p12.pw)
   diag(pairwise.P) <- 1
   if (CR.null == TRUE){
      names(list.zs) <- names(list.sds) <-c("No_Modules",list.names)
      pairwise.z <- as.matrix(z12.pw); rownames(pairwise.z) <- colnames(pairwise.z)<-c("No_Modules",list.names)
      pairwise.P <- as.matrix(p12.pw); rownames(pairwise.P) <- colnames(pairwise.P)<-c("No_Modules",list.names)
      diag(pairwise.P) <- 1
   }
   comment <- c("NOTE: more negative effects represent stronger modular signal!")
   out <- list(comment=comment, sample.z = list.zs,
               sample.r.sd = list.sds,
               pairwise.z = abs(pairwise.z),
               pairwise.P = pairwise.P)
   class(out) <- "compare.CR"
   out
 }
 