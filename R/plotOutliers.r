#' Find potential outliers
#' 
#' Function plots all specimens ordered by distance from the mean.
#' 
#' The function creates a plot of all specimens ordered by their Procrustes distance from the mean shape. 
#' The median distance (unbroken line) and upper and lower quartiles (dashed lines) summarize the distances
#' from the mean shape. Specimens falling above the upper quartile are plotted in red. The addresses of all specimens are
#' returned in the order displayed in the plot for further inspection by \code{\link{plotRefToTarget}}.
#' 
#' If the data have strong group structure and there is reasonable belief that the whole sample mean should not be used,
#' then a factor defining the groups can be used.
#' 
#' @param A An array (p x k x n) containing GPA-aligned coordinates for a set of specimens
#' @param groups An optional factor defining groups
#' @export
#' @keywords utilities
#' @seealso  \code{\link{gpagen}}
#' @seealso  \code{\link{plotTangentSpace}}
#' @seealso  \code{\link{plotAllSpecimens}}
#' @author Emma Sherratt
#' @return Function returns the landmark addresses of all specimens ordered as in the plot. If groups are used, function returns 
#' a list structure and a plot for each level in groups.
#' @examples
#' data(plethodon)
#' # let's make some outliers
#' newland <- plethodon$land
#' newland[c(1,8),,2] <- newland[c(8,1),,2]
#' newland[c(3,11),,26] <- newland[c(11,3),,2]
#' Y<- gpagen(newland) # GPA
#' out <- plotOutliers(Y$coords) # function returns dimnames and address of all specimens ordered
#' plotRefToTarget(mshape(Y$coords), Y$coords[,,out[1]], method="vector", label=TRUE)
#' plotRefToTarget(mshape(Y$coords), Y$coords[,,out[2]], method="vector", label=TRUE)
#' 
plotOutliers <- function(A, groups = NULL){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(is.null(groups)){ groups = factor(rep("All Specimens",dim(A)[3]))}
  res <- lapply(levels(groups), function(j){
    mn <- matrix(t(mshape(A[,,which(groups==j)])), nrow=1) 
    A.d <- two.d.array(A[,,which(groups==j)])
    d <- NULL
    for(i in 1:nrow(A.d)){ d <- c(d, as.vector(dist(rbind(mn, A.d[i,])))) }
    if(is.null(dimnames(A.d)[[1]])) { dimnames(A.d)[[1]] <- as.character(seq(1:nrow(A.d)))}
    names(d) <- dimnames(A.d)[[1]] 
    D <- d[order(d, decreasing=TRUE)]
    Q <- summary(D)
    Med <- as.numeric(summary(D)[3])
    LL <- as.numeric(Q[2] - 1.5*(Q[5]-Q[2]))
    UL <- as.numeric(Q[5] + 1.5*(Q[5]-Q[2]))
    plot(D, type="p", ylab= "Procrustes Distance from Mean", pch=19, xlab="", xaxt='n', main = j)
      abline(a=LL, b=0,lty=2,col= "blue")
      abline(a=Med,b=0,col= "blue")
      abline(a=UL,b=0,lty=2,col= "blue")
      text(x= nrow(A.d), y=LL, labels= "lower quartile", col = "blue", cex=0.5)
      text(x= nrow(A.d), y=Med, labels= "median",col = "blue", cex=0.5)
      text(x= nrow(A.d), y=UL, labels= "upper quartile",col = "blue", cex=0.5)
    if(any(D >= UL)) { 
      points(D[which(D >= UL)], pch=19, col="red")
      text(D[which(D >= UL)], labels=names(D)[which(D >= UL)], col= "red", adj=0.8, pos=4, cex=0.5)
    } else { text(D, labels=names(D), adj=c(0.5, 0.1), pos=4, cex=0.5)}
    ordered<-match(D,d)        
    names(ordered) <- names(D)
    return(ordered)
  })
  names(res) <- levels(groups)
  if(length(levels(groups))==1){ res <- res$`All Specimens`}
  return(res)

}
