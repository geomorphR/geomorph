#' Convert geomorph plots to ggplot objects
#'
#' Function attempts to coerce plot information from a geomorph plot object to an 
#' amenable ggplot object.  
#'
#' This function will attempt to use the plot arguments from an geomorph plot object 
#' to make a ggplot that can be additionally updated, as desired.  Not all plot 
#' characteristics might be converted.  Nonetheless, a ggplot will be coerced and could
#' be updated, according to user preference.
#' 
#' This function assumes no responsibility for arguments made by \code{\link[ggplot2]{ggplot}}.
#' It merely produces a ggplot object that should resemble a geomorph plot default.  Any 
#' augmentation of ggplot objects can be done either by direct intervention of the ggplot 
#' produced or reformatting the initial geomorph plot produced.  One should not expect direct
#' correspondence between R base plot parameters and ggplot parameters.  
#' 
#' @param object A plot object produced from \code{\link{plot.gm.prcomp}}, 
#' \code{\link{plot.pls}},\code{\link{plot.procD.lm}}, or \code{\link{plotAllometry}}.
#' For \code{\link{plot.procD.lm}} objects, only types "PC" or "regression"
#' should work.
#' @keywords utilities
#' @export
#' @author Michael Collyer
#' @examples
#' \dontrun{
#' 
#' ### PLS Example
#'  data(plethodon) 
#'  Y.gpa <- gpagen(plethodon$land)    #GPA-alignment    
#'  landmarks on the skull and mandible assigned to partitions
#'  land.gps <- c("A","A","A","A","A","B","B","B","B","B","B","B") 
#'  IT <- integration.test(Y.gpa$coords, partition.gp = land.gps, iter = 999)
#'  summary(IT) # Test summary
#'  P <- plot(IT) # PLS plot
#'  make_ggplot(P) # same plot in ggplot
#' 
#' ### Allometry example
#' 
#'  data(plethodon) 
#'  Y.gpa <- gpagen(plethodon$land, print.progress = FALSE)    #GPA-alignment  
#'
#'  gdf <- geomorph.data.frame(Y.gpa, site = plethodon$site, 
#'                            species = plethodon$species) 
#' 
#'  fit <- procD.lm(coords ~ Csize * species * site, data=gdf, iter=0, 
#'                  print.progress = FALSE)
#' 
#'  P <- plotAllometry(fit, size = gdf$Csize, logsz = TRUE, method = "PredLine", 
#'                      pch = 19, col = as.numeric(interaction(gdf$species, gdf$site)))
#'
#'  make_ggplot(P)
#' 
#' ### Tangent Space plot
#' 
#'  data(plethspecies) 
#'  Y.gpa <- gpagen(plethspecies$land)    #GPA-alignment
#' 
#'  PCA.w.phylo <- gm.prcomp(Y.gpa$coords, phy = plethspecies$phy)
#'  P <- plot(PCA.w.phylo, phylo = TRUE, main = "PCA.w.phylo")
#'  make_ggplot(P)
#' }
make_ggplot <- function(object){
  
  x <- y <- lbl <- NULL
  pa <- object$plot_args
  
  if(is.null(pa)) 
    stop("Plot object does not appear to be a type that can be converted to ggplot.\n",
         call. = FALSE)
  
  df <- data.frame(x = pa$x, y = pa$y)
  
  if(is.null(pa$cex)) pa$cex <- 1
  if(is.null(pa$pch)) pa$pch <- 19
  if(is.null(pa$bg)) pa$bg <- NA
  if(is.null(pa$col)) pa$col <- 1
  if(is.null(pa$lty)) pa$lty <- 1
  if(is.null(pa$lwd)) pa$lwd <- 1
  
  if(!is.null(pa$xlim)) xLim <- pa$xlim
  if(!is.null(pa$ylim)) yLim <- pa$ylim
  if(is.null(pa$xlim)) xLim <- c(min(df$x), max(df$x))
  if(is.null(pa$ylim)) yLim <- c(min(df$y), max(df$y))
  
  # 5% span boost for graphics
  xboost <- 0.05 *(xLim[2] - xLim[1])
  yboost <- 0.05 *(yLim[2] - yLim[1])
  
  xLim[1] <- xLim[1] - xboost * abs(xLim[1])
  xLim[2] <- xLim[2] + xboost * abs(xLim[2])
  yLim[1] <- yLim[1] - yboost * abs(yLim[1])
  yLim[2] <- yLim[2] + yboost * abs(yLim[2])
  
  if(!is.null(pa$asp)){
    g <- ggplot(df, aes(x, y)) + 
      labs(title = pa$main, x = pa$xlab, y = pa$ylab) +
      xlim(xLim) + ylim(yLim) + coord_fixed(ratio = pa$asp)
  } else {
    g <- ggplot(df, aes(x, y)) + 
      labs(title = pa$main, x = pa$xlab, y = pa$ylab) +
      xlim(xLim) + ylim(yLim) 
    
  }

  if(!is.null(object$phylo)) {
    ppa <- object$phylo$phylo.par
    dfp <- as.data.frame(object$phylo$phy.pcdata)
    colnames(dfp) <- c("x", "y")
    tree <- object$phylo$phy
    N <- length(tree$tip.label)
    edges <- as.matrix(tree$edge)
    tip.labeled <- ppa$tip.labels
    dfp$lbl <- rownames(dfp)
    node.labeled <- ppa$node.labels
    
    for(i in 1:nrow(edges)) {
      pts <- dfp[edges[i,], ]
      g <- g + geom_path(data = pts, aes(x, y), color = ppa$edge.color) +
        geom_point(data = pts, shape = ppa$node.pch, size = ppa$node.cex, 
                   fill = ppa$node.bg)
    }
    
    if(tip.labeled)
      g <- g + geom_text(data = dfp[1:N,], 
                         aes(x = x, y = y, 
                             label=lbl), 
                         size = ppa$tip.txt.cex * 2,
                         color = ppa$tip.txt.col,
                         hjust = ppa$tip.txt.adj[1] * 2, 
                         vjust = ppa$tip.txt.adj[2] * 2)
    
    if(node.labeled)
      g <- g + geom_text(data = dfp[-(1:N),], 
                         aes(x = x, y = y, 
                             label=lbl), 
                         size = ppa$node.txt.cex * 2,
                         color = ppa$node.txt.col,
                         hjust = ppa$node.txt.adj[1] * 2, 
                         vjust = ppa$node.txt.adj[2] * 2)
    
  }
  
  pts <- data.frame(x = pa$x, y = pa$y)
  
  g <- g + geom_point(data = pts, aes(x, y), 
                      shape = pa$pch, size = pa$cex * 5, col = pa$col, fill = pa$bg) 
  g
}