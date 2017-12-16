% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectSpecies.R
\name{selectSpecies}
\alias{selectSpecies}
\title{Select species based on traits}
\usage{
selectSpecies(t2c, constraints, t2d, obj = "QH", phi = 0.5,
  traitConstraint = TRUE, capd = FALSE)
}
\arguments{
\item{t2c}{Traits to constrain: A matrix of species trait values. Organize species as rows and traits as columns. 5 traits maximum.}

\item{constraints}{Trait constraints: A vector of trait values that serve as constants in the constraint equations. 5 constraints maximum. Must be listed in same order as columns in t2c.}

\item{t2d}{Traits to diversify: A matrix of species trait values to diversify. Organize species as rows, traits as columns. Can be any dimension (there is no upper limit to the number of traits to diversify).}

\item{obj}{Objective function: The objective function to optimise, one of three possibilities = c("QH", "Q", "H"). QH = Quadratic entropy (Q) plus entropy (H'); Q = Quadratic entropy; H = entropy.}

\item{phi}{A parameter bounded between 0 and 1 that weights the importance of either quadratic entropy or entropy (default = 0.5). Phi of 1 corresponds to 100 percent Q, phi of 0.5 corresponds to 50 percent Q and 50 perfect H', phi of 0 corresponds to 100 percent H'.}

\item{traitConstraint}{A logical stating whether solutions should be constrained to a trait mean (default = TRUE), when TRUE, a vector of constraints must be provided as an argument.}

\item{capd}{A logical stating whether the distance matrix should be capped at the mean distance among species (default = FALSE). Mean distance is calculated as the average of all upper triangular entries in the distance matrix calculated from t2d.}
}
\value{
A list with the elements: \item{prob}{Probabilities, i.e. optimal solutions} \item{cwm}{Final moment of constraint computed as prob x t2c using matrix multiplication.} \item{objval}{Values of the objective function being maximized. The last value is the maximum.} \item{lagrange}{Lagrange multipliers.} \item{hessian}{The Hessian at the optimal solution.}
}
\description{
This function returns a probability distribution for a species pool based on their traits and a desired trait profile (Laughlin 2014). It can simultaneously constrain specific trait value(s) and optimise functional diversity.
}
\examples{
### 1 trait constraint with maximum entropy
Spp=5 #S = number of species
trait <- as.matrix(data.frame(trait=c(1:Spp)))
rownames(trait)=c(letters[1:nrow(trait)])
result1 = selectSpecies(t2c=trait, constraints=c(3.5), t2d=trait, obj="H", capd=FALSE)
### compare result1 with virtually identical maxent output from FD package
#FD::maxent(constr=c(3.5),states=trait)$prob
### 1 trait constraint with maximum functional diversity
result2 = selectSpecies(t2c=trait, constraints=c(3.5), t2d=trait, obj="Q", capd=FALSE)
### 1 trait constraint with maximum functional diversity and entropy
result3 = selectSpecies(t2c=trait, constraints=c(3.5), t2d=trait, obj="QH", capd=FALSE)

### Plot results
plotProbs(result1,trait)
plotProbs(result2,trait)
plotProbs(result3,trait)

### 1 trait and no trait constraint
result4 = selectSpecies(t2c=trait, t2d=trait, obj="QH", traitConstraint=FALSE, capd=FALSE)
plotProbs(result4,trait)

##### 2 traits: Constrain trait X at X=3, diversify trait Y
trait.matrix <- as.matrix(cbind(traitX=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)),
  traitY=c(rep(c(1,2,3,4),4))))
rownames(trait.matrix)=c(letters[1:16])
result5 = selectSpecies(t2c=as.matrix(trait.matrix[,1]),constraints=c(3),
  t2d=as.matrix(trait.matrix[,2]),obj="Q",capd=FALSE)
result6 = selectSpecies(t2c=as.matrix(trait.matrix[,1]),constraints=c(3),
  t2d=as.matrix(trait.matrix[,2]),obj="QH",capd=TRUE)
plotProbs(result5,trait.matrix)
plotProbs(result6,trait.matrix)

##### 3 traits: Constrain trait Z to value 2.5, diversify trait X and Y
traitZ <- as.matrix(data.frame(c(1,3,2,2,3,1,2,3,1,2,1,3,2,3,2,2)))
result7 = selectSpecies(t2c=traitZ,constraints=c(2.5),t2d=trait.matrix, capd=TRUE, obj="QH")
plotProbs(result7,trait.matrix)

}
\references{
Laughlin, D.C. (2014) Applying trait-based models to achieve functional targets for theory-driven ecological restoration. Ecology Letters, 17, 771-784.
}