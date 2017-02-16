#' Plots results from the select function
#'
#' This function plots results (species probabilities/optimum solutions) from the select() function
#'
#' @param result A saved object from function select()
#' @param traits A matrix of trait values where traits are columns and rows are species. Maximum of 2 traits.
#' @return 2D barplot of probabilities for each species or 3D barplot that illustrates probabilities of species located within a 2D trait space
#' @export
#' @examples
#' ### 1 trait constraint with maximum functional diversity and entropy
#' Spp=5 #S = number of species
#' trait <- as.matrix(data.frame(trait=c(1:Spp)))
#' rownames(trait)=c(letters[1:nrow(trait)])
#' result1 = select(t2c=trait, constraints=c(3.5), t2d=trait, obj="QH", capd=FALSE)
#' plotprobs(result1,trait)
#'
#' ##### 2 traits: Constrain trait X to value 2.5, diversify trait Y
#' trait.matrix <- as.matrix(cbind(traitX=c(rep(1,3),rep(2,3),rep(3,3)),
#'                 traitY=c(rep(c(1,2,3),3))))
#' rownames(trait.matrix)=c(letters[1:9])
#' result2 = select(t2c=as.matrix(trait.matrix[,1]),constraints=c(2.5),
#'           t2d=as.matrix(trait.matrix[,2]),capd=TRUE,obj="QH")
#' plotprobs(result2,trait.matrix)

plotprobs = function(result,traits){

  res <- result$prob

  cols<-function(n) {grDevices::colorRampPalette(c("lightblue", "blue"))(20) } # 20 distinct colors

  if(ncol(traits)=="1"){
    graphics::barplot(t(res),ylim=c(0,max(res)),col="blue",names=rownames(res),ylab="Probability")
  }

  if(ncol(traits)=="2"){
  lattice::cloud(res~as.vector(t(traits[,1]))+as.vector(t(traits[,2])),
              panel.3d.cloud=latticeExtra::panel.3dbars, perspective=TRUE,
              xbase=0.5, ybase=0.5, scales=list(arrows=FALSE, col=1),  distance=0.3,
              xlim=c(min(traits[,1])-1,max(traits[,1])+1),
              ylim=c(min(traits[,2])-1,max(traits[,2])+1), zlim=c(0,max(res)),
              par.settings = list(axis.line = list(col = "transparent")),
              xlab=list("Trait X",rot=50,cex=1.5),ylab=list("Trait Y",rot=-30,cex=1.5),
              zlab=list("Probability",rot=90,cex=1.5), screen = list(z = 55, x = -55),
              col.facet = lattice::level.colors(res, at = lattice::do.breaks(range(res), 20),
              col.regions = cols, colors = TRUE))
  }
}
