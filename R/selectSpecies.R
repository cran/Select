#' Select species based on traits
#'
#' This function returns a probability distribution for a species pool based on their traits and a desired trait profile (Laughlin 2014). It can simultaneously constrain specific trait value(s) and optimise functional diversity.
#'
#' @param t2c Traits to constrain: A matrix of species trait values. Organize species as rows and traits as columns. 5 traits maximum.
#' @param constraints Trait constraints: A vector of trait values that serve as constants in the constraint equations. 5 constraints maximum. Must be listed in same order as columns in t2c.
#' @param t2d Traits to diversify: A matrix of species trait values to diversify. Organize species as rows, traits as columns. Can be any dimension (there is no upper limit to the number of traits to diversify).
#' @param obj Objective function: The objective function to optimise, one of three possibilities = c("QH", "Q", "H"). QH = Quadratic entropy (Q) plus entropy (H'); Q = Quadratic entropy; H = entropy.
#' @param phi A parameter bounded between 0 and 1 that weights the importance of either quadratic entropy or entropy (default = 0.5). Phi of 1 corresponds to 100 percent Q, phi of 0.5 corresponds to 50 percent Q and 50 perfect H', phi of 0 corresponds to 100 percent H'.
#' @param traitConstraint A logical stating whether solutions should be constrained to a trait mean (default = TRUE), when TRUE, a vector of constraints must be provided as an argument.
#' @param capd A logical stating whether the distance matrix should be capped at the mean distance among species (default = FALSE). Mean distance is calculated as the average of all upper triangular entries in the distance matrix calculated from t2d.
#' @return A list with the elements: \item{prob}{Probabilities, i.e. optimal solutions} \item{cwm}{Final moment of constraint computed as prob x t2c using matrix multiplication.} \item{objval}{Values of the objective function being maximized. The last value is the maximum.} \item{lagrange}{Lagrange multipliers.} \item{hessian}{The Hessian at the optimal solution.}
#' @examples
#' ### 1 trait constraint with maximum entropy
#' Spp=5 #S = number of species
#' trait <- as.matrix(data.frame(trait=c(1:Spp)))
#' rownames(trait)=c(letters[1:nrow(trait)])
#' result1 = selectSpecies(t2c=trait, constraints=c(3.5), t2d=trait, obj="H", capd=FALSE)

#' ### compare result1 with virtually identical maxent output from FD package
#' #FD::maxent(constr=c(3.5),states=trait)$prob

#' ### 1 trait constraint with maximum functional diversity
#' result2 = selectSpecies(t2c=trait, constraints=c(3.5), t2d=trait, obj="Q", capd=FALSE)

#' ### 1 trait constraint with maximum functional diversity and entropy
#' result3 = selectSpecies(t2c=trait, constraints=c(3.5), t2d=trait, obj="QH", capd=FALSE)
#'
#' ### Plot results
#' plotProbs(result1,trait)
#' plotProbs(result2,trait)
#' plotProbs(result3,trait)
#'
#' ### 1 trait and no trait constraint
#' result4 = selectSpecies(t2c=trait, t2d=trait, obj="QH", traitConstraint=FALSE, capd=FALSE)
#' plotProbs(result4,trait)
#'
#' ##### 2 traits: Constrain trait X at X=3, diversify trait Y
#' trait.matrix <- as.matrix(cbind(traitX=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)),
#'   traitY=c(rep(c(1,2,3,4),4))))
#' rownames(trait.matrix)=c(letters[1:16])
#' result5 = selectSpecies(t2c=as.matrix(trait.matrix[,1]),constraints=c(3),
#'   t2d=as.matrix(trait.matrix[,2]),obj="Q",capd=FALSE)
#' result6 = selectSpecies(t2c=as.matrix(trait.matrix[,1]),constraints=c(3),
#'   t2d=as.matrix(trait.matrix[,2]),obj="QH",capd=TRUE)
#' plotProbs(result5,trait.matrix)
#' plotProbs(result6,trait.matrix)
#'
#' ##### 3 traits: Constrain trait Z to value 2.5, diversify trait X and Y
#' traitZ <- as.matrix(data.frame(c(1,3,2,2,3,1,2,3,1,2,1,3,2,3,2,2)))
#' result7 = selectSpecies(t2c=traitZ,constraints=c(2.5),t2d=trait.matrix, capd=TRUE, obj="QH")
#' plotProbs(result7,trait.matrix)
#'
#' @references Laughlin, D.C. (2014) Applying trait-based models to achieve functional targets for theory-driven ecological restoration. Ecology Letters, 17, 771-784.
#' @export

selectSpecies = function(t2c, constraints, t2d, obj="QH", phi=0.5, traitConstraint=TRUE, capd=FALSE){

  if(ncol(t2c) >= nrow(t2c)-1) stop("There are more traits than species (i.e. T >= S-1). Use fewer traits or more species.")

  d = as.matrix(stats::dist(t2d, upper=TRUE, diag=TRUE))
  N = nrow(d)
  mean.dist <- mean(d[upper.tri(d,diag=FALSE)])

  dcap.fun = function(d){ d[d>mean.dist]=mean.dist; return(d) }
  if(capd==TRUE){d=dcap.fun(d)}

  mean.t2d <- mean(t2d)
  tdist <- abs(t2d-mean.t2d)
  for(i in 1:length(tdist)){if(tdist[i]==0){tdist[i]<-1e-3}}

  qh =function(p){ -(t(p)%*%(d/2)%*%p*phi + -(t(p)%*%log(p))*(1-phi)) }  #Q+H'
  q = function(p){ -(t(p)%*%(d/2)%*%p) }                                 #Q only
  h = function(p){ -(-(t(p)%*%log(p))) }                                 #H' only

  if(traitConstraint==TRUE){
    if(length(constraints)!=ncol(t2c)) stop("The number of constraints must be equal to the number of traits to constrain.")

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

  else if(ncol(t2c)==3){
    eqfun=function(p){
      z1=sum(p)
      z2=t2c[,1]%*%p
      z3=t2c[,2]%*%p
      z4=t2c[,3]%*%p
      return(c(z1,z2,z3,z4))}}

  else if(ncol(t2c)==4){
    eqfun=function(p){
      z1=sum(p)
      z2=t2c[,1]%*%p
      z3=t2c[,2]%*%p
      z4=t2c[,3]%*%p
      z5=t2c[,4]%*%p
      return(c(z1,z2,z3,z4,z5))}}

  else if(ncol(t2c)==5){
    eqfun=function(p){
      z1=sum(p)
      z2=t2c[,1]%*%p
      z3=t2c[,2]%*%p
      z4=t2c[,3]%*%p
      z5=t2c[,4]%*%p
      z6=t2c[,5]%*%p
      return(c(z1,z2,z3,z4,z5,z6))}}

  all.constraints = c(1,constraints)
  }

  if(traitConstraint==FALSE){
    eqfun=function(p){
      z1=sum(p)
      return(c(z1))}
    all.constraints=c(1)
    t2c=c(rep(1,N))}

  oldw <- getOption("warn") # save default warnings
  options(warn = -1)        # Rsolnp::solnp tries to create object "tempdf", this fix overrides that

  if(obj=="QH")
  {res=Rsolnp::solnp(pars=c(rep(1/N,N)),fun=qh,eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))}
  else if(obj=="Q")
  {res=Rsolnp::solnp(pars=c(rep(1/N,N)),fun=q, eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))}
  else if(obj=="H")
  {res=Rsolnp::solnp(pars=c(rep(1/N,N)), fun=h, eqfun, eqB=all.constraints, LB = c(rep(0,N)), UB = c(rep(1,N)))}

  result = list()
    result$prob <- as.matrix(res$pars); rownames(result$prob) <- rownames(t2c)
    result$cwm <- as.matrix(res$pars%*%t2c)
    result$objval <- as.matrix(res$values)
    result$lagrange <- as.matrix(res$lagrange)
    result$hessian <- as.matrix(res$hessian)

  return(result)
  options(warn = oldw) # reinstate warnings
}
