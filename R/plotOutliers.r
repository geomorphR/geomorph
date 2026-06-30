#' Find potential outliers
#' 
#' Function plots all specimens ordered by distance from the mean.
#' 
#' The function creates a plot of all specimens ordered by their distance from the mean shape. There
#' is an option for using specified principal components (PC) of shape.  If all shape dimensions are used
#' then distances are equal to Procrustes distances or distances in the tangent space of shape space that 
#' resemble Procrustes distances, depending on whether the projection was performed with generalized 
#' Procrustes analysis (GPA).  Once distances are calculated, a power-transformation is perfomed to 
#' normalize the distances.  From these data, an upper limit is estimated following Tukey's box-plot rule, 
#' as \eqn{Q_3 + 1.5 \times (Q_3 - Q_1)}, where \eqn{Q} refers to quartile.  This upper limit is back-transformed 
#' to distance, and any distances greater than this limit are colored red.  (Note: These shapes could 
#' be considered outliers but their red color does not mean they are necessarily outliers. )
#' 
#'  The user may optionally 
#' also inspect the shapes of identified configurations as compared to the consensus, in order
#' to identify potential digitization errors or other data issues. The addresses of all specimens are
#' returned in the order displayed in the plot for further inspection by \code{\link{plotRefToTarget}}.
#' 
#' If the data have strong group structure and there is reasonable belief that the whole sample mean should not be used,
#' then a factor defining the groups can be used.
#' 
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of specimens
#' @param groups An optional factor defining groups
#' @param PC A single number or range of principal components that can be used to check for outliers.  If
#' NULL, all shape dimensions are used.  This argument might be useful for investigating subtle but
#' important aspects of shape that would not elucidate outliers in all dimensions.
#' @param inspect.outliers A logical value indicating whether to plot outlier shape configurations as compared to the consensus
#' @export
#' @keywords utilities
#' @seealso  \code{\link{gpagen}}
#' @seealso  \code{\link{plotAllSpecimens}}
#' @author Emma Sherratt, Antigoni Kaliontzopoulou, & Michael Collyer
#' @return Function returns the landmark addresses of all specimens ordered as in the plot. If groups are used, function returns 
#' a list structure and a plot for each level in groups.
#' @examples
#' \dontrun{
#' data(plethodon)
#' # let's make some outliers
#' newland <- plethodon$land
#' newland[c(1,8),,2] <- newland[c(8,1),,2]
#' newland[c(3,11),,26] <- newland[c(11,3),,2]
#' Y <- gpagen(newland) # GPA
#' out <- plotOutliers(Y$coords) # function returns dimnames and address 
#' out
#' # of all specimens ordered
#' plotOutliers(Y$coords, inspect.outliers = TRUE) # function also produces 
#' # plots of identified outlier specimens compared to the mean shape
#' 
#' # example with groups
#' plotOutliers(Y$coords, groups = plethodon$species, 
#' inspect.outliers = TRUE)
#' 
#' # previous example using first three PCs of shape
#' plotOutliers(Y$coords, groups = plethodon$species, 
#' PC = 1:3, inspect.outliers = TRUE)
#' 
#' # previous example using just the first PC of shape
#' plotOutliers(Y$coords, groups = plethodon$species, 
#' PC = 1, inspect.outliers = TRUE)
#'  }
plotOutliers <- function(A, groups = NULL, 
                         PC = NULL,
                         inspect.outliers = FALSE){
  if(length(dim(A)) != 3)
    stop("\nData matrix not a 3D array (see 'arrayspecs').\n",
         call. = FALSE)  
  
  if(is.null(groups))
    groups <- factor(rep("All Specimens",dim(A)[3]))
  
  glen <- as.vector(table(groups))
  if(any(glen < 3))
    stop("\n One or more groups have fewer than 3 shapes.  Outliers cannot be deteced.\n",
         call. = FALSE)
  options(viewer = NULL)  # for 3D plotting
  Ymat <- gm.prcomp(A)$x
  if(is.null(PC)) PC <- 1:ncol(Ymat)
  if(any(PC > ncol(Ymat)))
    stop("\nChoose a number of PCs that does not exceed the total number of PCs\n.",
         call. = FALSE)
  Ymat <- as.matrix(Ymat[, PC])
  n <- NROW(Ymat)
  
  if(is.null(dimnames(A)[[3]]))
    dimnames(A)[[3]] <- rownames(Ymat) <- 1:n

  res <- lapply(levels(groups), function(j){
    y <- as.matrix(center(Ymat[which(groups == j), ]))
    D <- sqrt(rowSums(y^2))
    D <- D[order(D, decreasing=TRUE)]
    bc <- powerTrans(D)
    b <- bc$transformed
    UL <- quantile(b, 0.75) + 1.5 * IQR(b)
    UL <- (UL * bc$opt.lambda + 1)^(1/bc$opt.lambda)

    plot(D, type="p", ylab= "Distance from Mean", pch=19, xlab="", xaxt='n', main = j)
      abline(a=UL,b=0,lty=2,col= "blue")
      text(x= nrow(y), y=UL, labels= "upper limit",col = "blue", cex=0.5, adj=c(0.5, -0.5))
    if(any(D >= UL)) { 
      ol <- D[which(D >= UL)]
      p <- 1:length(ol)
      points(p, ol, pch=19, col="red")
      text(p, ol, labels = names(ol), col= "red", adj=0.8, pos=4, cex=0.5)
    }
      
    if(any(D >= UL)) {   
      if(inspect.outliers == TRUE){
        out.config <- names(ol)
        Alist <- lapply(1:n, function(x) as.matrix(A[,, x]))
        names(Alist) <- dimnames(A)[[3]]
        OLlist <- Alist[names(ol)]
        M <- mshape(A)
        
        if(dim(A)[2] == 2){
          for(i in seq_along(OLlist)){
            plotRefToTarget(M, OLlist[[i]], method="vector", label=TRUE)
            title(main = paste("group:", j, ", specimen:", out.config[i]))
          }
        }
        
        if(dim(A)[2] == 3){
          for(i in seq_along(OLlist)){
            fig <- plotRefToTarget(M, OLlist[[i]], 
                                   method="vector", label=TRUE)
            fig <- fig |>
              plotly::layout(title = paste("group:", j,
                                   "<br>specimen:", out.config[i]))
            print(fig)
          }
        }
      }
    } else { text(1:n, D, labels=names(D), adj=c(0.5, 0.1), pos=4, cex=0.5)}
    names(D)
  })
  names(res) <- levels(groups)
  if(length(glen) == 1) res <- res$`All Specimens`
  res
}
