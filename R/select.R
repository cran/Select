#' Select species based on traits
#'
#' This function returns a probability distribution for a species pool based on their traits and a desired trait profile (Laughlin 2014).
#'
#' @param t2c Traits to constrain: A matrix of species trait values. Can be up to 2 dimensions. Organize species as rows and traits as columns.
#' @param constraints Trait constraints: A vector of trait values that serve as constants in the constraint equations
#' @param t2d Traits to diversify: A matrix of species trait values to diversify. Organize species as rows, traits as columns. Can be any dimension.
#' @param obj Objective function: The objective function to optimise, one of four possibilities = c("FDandE", "FD", "E", "FE").
#' FDandE = Rao functional diversity plus entropy; FD = Rao functional diversity only; E = entropy only; FE = functional entropy.
#' @param phi A parameter bounded between 0 and 1 that weights the importance of either entropy or Rao (default = 0.5)
#' @param traitConstraint A logical stating whether solutions should be constrained to a trait mean (default = TRUE),
#' when TRUE, a vector of constraints must be provided as an argument.
#' @param capd A logical stating whether the distance matrix should be capped at the mean distance among species (default = FALSE).
#' @return Optimal solutions (i.e. probabilities) for each species.
#' @examples
#' ##### 1 trait with a constraint
#' trait <- as.matrix(data.frame(trait=c(1,2,3)))
#' rownames(trait)=c(letters[1:nrow(trait)])
#' result1 = select(t2c=trait, constraints=c(2.5), t2d=trait, obj="FDandE", capd=TRUE)
#' barplot(t(result1),ylim=c(0,0.7),col="blue",names=c(c(1,2,3)),ylab="Relative abundance")
#'
#' ##### 1 trait and no trait constraint
#' result2 = select(t2c=trait, t2d=trait, obj="FDandE", traitConstraint=FALSE, capd=TRUE)
#' barplot(t(result2),ylim=c(0,0.7),col="blue",names=c(c(1,2,3)),ylab="Relative abundance")
#'
#' ##### 2 traits: Constrain trait X to value 2.5, diversify trait Y
#' trait.matrix <- as.matrix(cbind(traitX=c(rep(1,3),rep(2,3),rep(3,3)),
#'                 traitY=c(rep(c(1,2,3),3))))
#' rownames(trait.matrix)=c(letters[1:9])
#' result3 = select(t2c=as.matrix(trait.matrix[,1]),constraints=c(2.5),
#'           t2d=as.matrix(trait.matrix[,2]),capd=TRUE,obj="FDandE")
#' plotprobs(result3,trait.matrix)
#'
#' ##### 3 traits: Constrain trait Z to value 2.5, diversify trait X and Y
#' traitZ <- as.matrix(data.frame(c(1,3,2,2,3,1,2,3,1)))
#' result4 = select(t2c=traitZ,constraints=c(2.5),t2d=trait.matrix, capd=TRUE, obj="FDandE")
#' plotprobs(result4,trait.matrix)
#'
#' @references Laughlin, D.C. (2014) Applying trait-based models to achieve functional targets for theory-driven ecological restoration. Ecology Letters, 17, 771-784.
#' @export

select = function(t2c, constraints, t2d, obj="FDandE", phi=0.5, traitConstraint=TRUE, capd=FALSE){
  if(ncol(t2d)>1 & obj=="FE") stop("function cannot optimise fe with multidimensional traits, choose different obj")

  d = as.matrix(stats::dist(t2d,upper=TRUE, diag=TRUE))
  N = nrow(d)
  mean.dist <- mean(d[upper.tri(d,diag=FALSE)])

  dcap.fun = function(d){ d[d>mean.dist]=mean.dist; return(d) }
  if(capd==TRUE){d=dcap.fun(d)}

  mean.t2d <- mean(t2d)
  tdist <- abs(t2d-mean.t2d)
  for(i in 1:length(tdist)){if(tdist[i]==0){tdist[i]<-1e-3}}

  fde=function(p){ -(t(p)%*%(d/2)%*%p*phi + -(t(p)%*%log(p))*(1-phi)) }  #Rao+Entropy
  fd= function(p){ -(t(p)%*%(d/2)%*%p) }                                 #Rao only
  e = function(p){ -(-(t(p)%*%log(p))) }                                 #Entropy only
  fe= function(p){ -(-((t(p)*t(tdist))%*%log(p*tdist))) }                #functional entropy, not great

  if(traitConstraint==TRUE){
  if(ncol(t2c)==1){
    eqfun=function(p){
      z1=sum(p)
      z2=t2c[,1]%*%p
      return(c(z1,z2))}}
  else if(ncol(t2c)==2){
    eqfun=function(p){
      z1=sum(p)
      z2=t2c[,1]%*%p
      z3=t2c[,2]%*%p
      return(c(z1,z2,z3))}}
    all.constraints = c(1,constraints)}

  if(traitConstraint==FALSE){
    eqfun=function(p){
      z1=sum(p)
      return(c(z1))}
    all.constraints=c(1)
    t2c=c(rep(1,N))}

  if(obj=="FDandE")
  {prob=Rsolnp::solnp(pars=c(rep(1/N,N)),fun=fde,eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))$pars}
  else if(obj=="FD")
  {prob=Rsolnp::solnp(pars=c(rep(1/N,N)),fun=fd, eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))$pars}
  else if(obj=="E")
  {prob=Rsolnp::solnp(pars=c(rep(1/N,N)), fun=e, eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))$pars}
  else if(obj=="FE")
  {prob=Rsolnp::solnp(pars=c(rep(1/N,N)), fun=fe, eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))$pars}

  prob=as.matrix(prob)
  rownames(prob)=rownames(t2c)

  return(prob)
}
