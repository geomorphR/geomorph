#' Create a data frame with shape data
#' 
#' A list similar to a data frame to facilitate analysis of shape data.
#'
#' This function produces a list that can be used like a data frame in other analytical functions.
#' The purpose is similar to the function, \code{\link[base]{data.frame}}, but without the constraint that 
#' data must conform to an n (observations) x p (variables) matrix.  Rather, the list produced is 
#' constrained only by n.  List objects can be Procrustes residuals (coordinates) arrays, matrices, variables,
#' distance matrices, and phylogenetic trees.  Results from \code{\link{gpagen}} can be directly
#' imported into a geomorph.data.frame to utilize the coordinates and centroid size as variables. (See Examples)
#' @param ... a list of objects to include in the data frame.
#' @export
#' @author Michael Collyer
#' @keywords utilities
#' @examples
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land,PrinAxes=FALSE)
#' gdf <- geomorph.data.frame(Y.gpa)
#' attributes(gdf)
#' 
#' gdf <- geomorph.data.frame(Y.gpa, species = plethodon$species, site = plethodon$site)
#' attributes(gdf)
#' 
#' # Using geomorph.data.frame to facilitate analysis
#' procD.lm(coords ~ Csize + species * site, data = gdf)
geomorph.data.frame <- function(...) {
dots <- list(...)
list.check0 <- sapply(1:length(dots), function(j) any(is.geomorph.data.frame(dots[[j]])))
list.check00 <- sapply(1:length(dots), function(j) any(is.data.frame(dots[[j]])))
dots0 <- unlist(dots[list.check0], recursive=FALSE)
dots00 <- unlist(dots[list.check00], recursive=FALSE)
dots.updated <- dots[!list.check0 & !list.check00]
if(length(dots.updated) > 0) {
  list.check1 <- sapply(1:length(dots.updated), function(j) is.gpagen(dots.updated[[j]]))
  dots1 <- dots.updated[list.check1]
} else dots1 <- NULL
if(length(dots1) > 0){
  dots2 <- lapply(1:length(dots1), function(j){
    x <- unlist(dots1[j], recursive = FALSE)
    list(x$coords, x$Csize)
  })
  dots2 <- unlist(dots2, recursive = FALSE)
  dots2.names <- rep(c("coords", "Csize"), length(dots1))
  names(dots2) <- dots2.names
  dots.updated <- dots.updated[!list.check1]
} else {
  dots2 <- NULL
}
if(length(dots.updated) > 0){
  list.check2 <- sapply(1:length(dots.updated), function(j) is.phylo(dots.updated[[j]]))
  dots3 <- dots.updated[list.check2]
  dots4 <- dots.updated[!list.check2]
} else {
  dots3 <- NULL
  dots4 <- NULL
}
if(length(dots3) == 0) dots3 <- NULL
if(length(dots4) == 0) dots4 <- NULL
dots <- c(dots0,dots00,dots2, dots3,dots4)
N <- length(dots)
dots.ns <- array(NA,N)
for(i in 1:N){
  if(is.array(dots[[i]])) {
    if(length(dim(dots[[i]])) == 3) dots.ns[i] <- dim(dots[[i]])[[3]]
    if(length(dim(dots[[i]])) == 2) dots.ns[i] <- dim(dots[[i]])[[2]]
    if(length(dim(dots[[i]])) == 1) dots.ns[i] <- dim(dots[[i]])[[1]]
  }
  if(is.matrix(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[1]]
  if(class(dots[[i]]) == "dist") dots.ns[i] <- attr(dots[[i]], "Size")
  if(class(dots[[i]]) == "phylo") dots.ns[i] <- length(dots[[i]]$tip.label)
  if(is.data.frame(dots[[i]])) dots.ns[i] <- dim(dots[[i]])[[2]]
  if(is.vector(dots[[i]])) dots.ns[i] <- length(dots[[i]])
  if(is.factor(dots[[i]])) dots.ns[i] <- length(dots[[i]])
  if(is.logical(dots[[i]])) dots.ns[i] <- length(dots[[i]])
}
if(any(is.na(dots.ns))) stop("Some input is either dimensionless or inappropriate for data frames")
if(length(unique(dots.ns)) > 1) stop("Inputs have different numbers of observations")
class(dots) <- c("geomorph.data.frame")
dots
}

