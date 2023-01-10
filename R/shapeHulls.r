#' Update Plots with Convex Hulls for Groups
#' 
#' This function is used to update \code{\link{plot.procD.lm}} and \code{\link{plot.gm.prcomp}} ordination plot
#' objects with convex hulls for different groups.  If no groups are defined (groups is NULL) just a single
#' convex hull will be returned.  Groups do not need to be a factor in the original \code{\link{procD.lm}} fit.
#' 
#' This function is a wrapper for the \code{\link{points}} function. It is intentionally limited, so
#' as to not interfere with other plot parameter adjustments.
#' 
#' @param x A \code{\link{plot.procD.lm}} or \code{\link{plot.gm.prcomp}} plot object.
#' @param groups An optional vector or factor to define groups for hull.  If NULL, only one hull will be generated for all points.
#' @param group.cols An optional vector to define hull colors, arranged in the same order as factor levels.  If NULL and if multiple groups
#' exist, the general R color sequence (black, red, green, blue, etc.) will be used.
#' @param group.lwd An optional vector equal in length to the number of group levels, and arranged in the order of group levels,
#' to modify hull line width.
#' @param group.lty An optional vector equal in length to the number of group levels, and arranged in the order of group levels,
#' to modify hull line type.
#' @export
#' @author Michael Collyer
#' @seealso \code{\link{procD.lm}}
#' @keywords utilities
#' @examples
#' 
#' # Via procD.lm and plot.procD.lm
#' 
#' data("pupfish")
#' gdf <- geomorph.data.frame(coords = pupfish$coords, Sex = pupfish$Sex,
#' Pop = pupfish$Pop)
#' fit <- procD.lm(coords ~ Pop * Sex, data = gdf, print.progress = FALSE)
#' pc.plot <- plot(fit, type = "PC", pch = 19)
#' shapeHulls(pc.plot)
#' 
#' pc.plot <- plot(fit, type = "PC", pch = 19)
#' groups <- interaction(gdf$Pop, gdf$Sex)
#' 
#' shapeHulls(pc.plot, groups = groups, 
#' group.cols = c("dark red", "dark red", "dark blue", "dark blue"),
#' group.lwd = rep(2, 4), group.lty = c(2, 1, 2, 1))
#' 
#' legend("topright", levels(groups), 
#' col = c("dark red", "dark red", "dark blue", "dark blue"),
#' lwd = rep(2,4), lty = c(2, 1, 2, 1))
#' 
#' pc.plot <- plot(fit, type = "PC", pch = 19)
#' shapeHulls(pc.plot, groups = gdf$Sex, group.cols = c("black", "black"), 
#' group.lwd = rep(2, 2), group.lty = c(2, 1))
#' legend("topright", levels(gdf$Sex), lwd = 2, lty = c(2, 1))
#' 
#' # Via gm.prcomp and plot.gm.prcomp
#' 
#' data(plethspecies) 
#' Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
#' pleth.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#' summary(pleth.phylo)
#' 
#' pc.plot <- plot(pleth.phylo, phylo = TRUE)
#' gp <- factor(c(rep(1, 5), rep(2, 4)))
#' shapeHulls(pc.plot, groups = gp, group.cols = 1:2, 
#' group.lwd = rep(2, 2), group.lty = c(2, 1))
#' legend("topright", c("P. cinereus clade", "P. hubrichti clade"), 
#' col = 1:2, lwd = 2, lty = c(2, 1))

shapeHulls <- function(x, groups = NULL, group.cols = NULL, 
                       group.lwd = NULL, group.lty = NULL){
  y <- as.matrix(x$PC.points)
  if(NCOL(y) < 2) stop("Cannot generate hulls in fewer than 2 dimensions")
  if(NCOL(y) > 2) y <- y[,1:2]
  n <- NROW(y)
  if(!is.null(groups) && length(groups) != n ) stop("Different number of observations between groups factor and PC plot.\n",
                                                    call. = FALSE)

  if(is.null(groups)) groups <- rep(1, n)
  groups <- as.factor(groups)
  if(length(unique(groups)) != length(levels(groups)))
    warning("The levels in the grouping factor do not match the number of unique factor levels.\n",
            call. = FALSE, immediate. = TRUE)
  ug <- unique(groups)
  g <- length(ug)
  if(is.null(group.cols)) group.cols <- 1:g
  if(is.null(group.lwd)) group.lwd <- rep(1, g)
  if(is.null(group.lty)) group.lty <- rep(1, g)
  if(length(group.cols) != g) stop("Number of requested group colors does not match the number of groups")
  if(length(group.lwd) != g) stop("Number of requested group widths does not match the number of groups")
  if(length(group.lty) != g) stop("Number of requested group line types does not match the number of groups")
  
  for(i in 1:g){
    yy <- y[groups == ug[i],]
    chp <- chull(yy)
    chp <- c(chp, chp[1])
    points(yy[chp,], type = "l", lty = group.lty[i],
           lwd = group.lwd[i], col = group.cols[i])
  }
}
