#' Plots results from select
#'
#' This function plots results (species probabilities/optimum solutions) from the select() function
#'
#' @param result A saved object from function select()
#' @param traits A 2D matrix of trait values
#' @return 3D barplot that illustrates probabilities of species located within a 2D trait space
#' @export
#' @examples
#' ##### 2 traits: Constrain trait X to value 2.5, diversify trait Y
#' trait.matrix <- as.matrix(cbind(traitX=c(rep(1,3),rep(2,3),rep(3,3)),
#'                 traitY=c(rep(c(1,2,3),3))))
#' rownames(trait.matrix)=c(letters[1:9])
#' result3 = select(t2c=as.matrix(trait.matrix[,1]),constraints=c(2.5),
#'           t2d=as.matrix(trait.matrix[,2]),capd=TRUE,obj="FDandE")
#' plotprobs(result3,trait.matrix)

plotprobs = function(result,traits){
  cols<-function(n) {grDevices::colorRampPalette(c("lightblue", "blue"))(20) } # 20 distinct colors
  lattice::cloud(result~as.vector(t(traits[,1]))+as.vector(t(traits[,2])),
              panel.3d.cloud=latticeExtra::panel.3dbars, perspective=TRUE,
              xbase=0.5, ybase=0.5, scales=list(arrows=FALSE, col=1),  distance=0.3,
              xlim=c(min(traits[,1])-1,max(traits[,1])+1),
              ylim=c(min(traits[,2])-1,max(traits[,2])+1), zlim=c(0,max(result)),
              par.settings = list(axis.line = list(col = "transparent")),
              xlab=list("traitX",rot=50,cex=1.5),ylab=list("traitY",rot=-30,cex=1.5),
              zlab=list("Abundance",rot=90,cex=1.5), screen = list(z = 55, x = -55),
              col.facet = lattice::level.colors(result, at = lattice::do.breaks(range(result), 20),
              col.regions = cols, colors = TRUE))
}
