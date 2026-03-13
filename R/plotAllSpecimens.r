#' Plot landmark coordinates for all specimens
#'
#' Function plots landmark coordinates for a set of specimens
#'
#' The function creates a plot of the landmark coordinates for all specimens. This is useful for examining 
#'  patterns of variation in Procrustes shape variables, after a GPA has been performed. If "mean = TRUE", the mean shape will be calculated and added to the plot.
#'  Additionally, if a matrix of links is provided, the landmarks of the mean shape will be connected by lines.  
#'  The link matrix is an m x 2 matrix, where m is the desired number of links. Each row of the link matrix 
#'  designates the two landmarks to be connected by that link. The function will plot either two- or 
#'  three-dimensional data.
#'
#' @param A A 3D array (p x k x n) containing Procrustes shape variables for a set of specimens
#' @param mean A logical value indicating whether the mean shape should be included in the plot
#' @param links An optional matrix defining for links between landmarks (only if mean=TRUE)
#' @param label A logical value indicating whether landmark numbers will be plotted (only if mean=TRUE)
#' @param plot_param A list of plot parameters for the points (pt.bg, pt.cex), mean (mean.bg, mean.cex), links (link.col, link.lwd, link.lty) and landmark labels (txt.cex, txt.adj, txt.pos, txt.col)
#' @export
#' @keywords visualization
#' @author Dean Adams
#' @examples
#' \dontrun{
#' data(plethodon) 
#' Y.gpa <- gpagen(plethodon$land)    #GPA-alignment
#'
#' plotAllSpecimens(Y.gpa$coords, links = plethodon$links)
#' }
plotAllSpecimens<-function(A,mean=TRUE,links=NULL,label=FALSE,plot_param = list()){
  if (length(dim(A))!=3){
    stop("Data matrix not a 3D array (see 'arrayspecs').")  }
  if(any(is.na(A))==T){
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")  }
  k <- dim(A)[2]
  p <- dim(A)[1]
  if(mean==TRUE){ mn<-mshape(A) }
  
  p.p <- plot_param
  if(is.null(p.p$pt.bg)) p.p$pt.bg="gray" ; if(is.null(p.p$pt.cex)) p.p$pt.cex=1 ; 
  if(is.null(p.p$mean.bg)) p.p$mean.bg="black" ; if(is.null(p.p$mean.cex)) p.p$mean.cex=2
  if(is.null(p.p$link.col)) p.p$link.col="black" ; if(is.null(p.p$link.lwd)) p.p$link.lwd=2
  if(is.null(p.p$link.lty)) p.p$link.lty=1 ; if(is.null(p.p$txt.adj)) p.p$txt.adj=c(-.1,-.1)
  if(is.null(p.p$txt.col)) p.p$txt.col="black" ; if(is.null(p.p$txt.cex)) p.p$txt.cex=0.8
  if(is.null(p.p$txt.pos)) p.p$txt.pos=1
  
  if(k==2){
    plot(A[,1,],A[,2,],asp=1, pch=21,bg=p.p$pt.bg,cex=p.p$pt.cex*1,xlab="x",ylab="y") 
    if(mean==TRUE){ 
      if(is.null(links)==FALSE){
        linkcol <- rep(p.p$link.col,nrow(links))[1:nrow(links)]
        linklwd <- rep(p.p$link.lwd,nrow(links))[1:nrow(links)]
        linklty <- rep(p.p$link.lty,nrow(links))[1:nrow(links)]
        for (i in 1:nrow(links)){
          segments(mn[links[i,1],1],mn[links[i,1],2],mn[links[i,2],1],mn[links[i,2],2],
                   col=linkcol[i],lty=linklty[i],lwd=linklwd[i])
        }
      }
      points(mn,pch=21,bg=p.p$mean.bg,cex=p.p$mean.cex)
      if(label == TRUE){text(mn, label=paste(1:dim(mn)[1]),adj=(p.p$txt.adj+p.p$mean.cex),
                             pos=p.p$txt.pos,cex=p.p$txt.cex,col=p.p$txt.col)}
    }
  }
  if(k==3){
    A3d<-NULL
    for (i in 1:dim(A)[[3]]){
      A3d<-rbind(A3d,A[,,i])
    }
    A3d <- data.frame(A3d)
    
    fig <- plot_ly()
    fig <- fig |>
      plotly::layout(scene = list(aspectmode = "data",
          xaxis = list(title = '', 
              showgrid = F, visible = F,showticklabels = F, 
              zeroline = F, showbackground = F),
          yaxis = list(title = '', showgrid = F, 
              visible = F, showticklabels = F, zeroline = F, 
              showbackground = F),
          zaxis = list(title = '', showgrid = F, 
              visible = F, showticklabels = F, zeroline = F, 
              showbackground = F)),showlegend = FALSE) |>
       add_trace(x = ~A3d$X, y = ~A3d$Y, z = ~A3d$Z, type = "scatter3d",
                 mode = "markers", name = "LM",
                 marker = list(color = p.p$pt.bg, 
                 size = p.p$pt.cex*1.5),
              showlegend = FALSE,
              inherit = FALSE) 
    if(mean==TRUE){ 
      mn_df <- as.data.frame(mn)
      colnames(mn_df) <- c("X","Y","Z")
      fig <- fig |>
      add_trace(x = ~mn_df$X, y = ~mn_df$Y, z = ~mn_df$Z, 
        type = "scatter3d", mode = "markers", 
        name = "Mean", marker = list(color = p.p$mean.bg, 
                                     size = p.p$mean.cex*2),
        showlegend = FALSE,
        inherit = FALSE) 
    } 
    if(is.null(links)==FALSE){ 
      dash_map <- c("solid", "dash", "dot", "dashdot", 
                    "longdash",  "longdashdot")
      line_df <- data.frame()
    for (i in 1:nrow(links)) {
      line_df <- rbind(
        line_df,
        mn_df[links[i, 1], ],
        mn_df[links[i, 2], ],
        data.frame(X = NA, Y = NA, Z = NA)
      )
    }
      fig <- fig |>
      add_trace(x = ~line_df$X, y = ~line_df$Y, z = ~line_df$Z,
          type = "scatter3d", mode = "lines",      
          line = list(color = p.p$link.col, width = p.p$link.lwd,
                      dash = dash_map[p.p$link.lty]),
          showlegend = FALSE,
          inherit = FALSE)    
    }  
    if(label == TRUE){ 
      fig <- fig |>
      add_trace(x = ~mn_df$X, y = ~mn_df$Y, z = ~mn_df$Z,
          type = "scatter3d", mode = "text",
          text = seq_len(p), 
          textfont = list(color = p.p$txt.col, size = p.p$txt.cex),
          textposition = "top center",
          showlegend = FALSE,
          inherit = FALSE)
    } 
    fig
  }
}
