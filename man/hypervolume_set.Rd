\name{hypervolume_set}
\alias{hypervolume_set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set operations (intersection / union / unique components)
}
\description{
Computes the intersection, union, and unique components of two hypervolumes.
}
\usage{
hypervolume_set(hv1, hv2, reduction_factor = 1, verbose = T, check_memory = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv1}{
A n-dimensional hypervolume
}
  \item{hv2}{
A n-dimensional hypervolume
}
  \item{reduction_factor}{
A number in (0,1] that represents the fraction of random points sampled from each hypervolume for the stochastic inclusion test. Larger values are more accurate but computationally slower.
}
  \item{verbose}{
Logical value; print diagnostic output if true.
}
  \item{check_memory}{
Logical value; returns information about expected memory usage if true.
}
}
\details{
Uses the inclusion test approach to identify points in the first hypervolume that are or are not within the second hypervolume and vice-versa. The intersection is the points in both hypervolumes, the union those in either hypervolume, and the unique components the points in one hypervolume but not the other.

By default, the function uses {check_memory=TRUE} which will provide an estimate of the computational cost of the set operations. The function should then be re-run with {check_memory=FALSE} if the cost is acceptable. This algorithm's memory and time cost scale quadratically with the number of input points, so large datasets can have disproportionately high costs. This error-checking is intended to prevent the user from large accidental memory allocation. 

If the computational cost is too high, a small value of {reduction_factor} can be specified to reduce the number of points randomly sampled for this stochastic analysis. Linear decreases in this parameter produce quadratic improvements in cost. 
}
\value{
If {check_memory} is false, returns a HypervolumeList object, with four items in its HVList slot:

  \item{Intersection }{The intersection of hv1 and hv2}
  \item{Union }{The union of hv1 and hv2}
  \item{Unique_1 }{The unique component of hv1 relative to hv2}
  \item{Unique_2 }{The unique component of hv2 relative to hv1}
  
  Note that the output hypervolumes will have lower random point densities than the input hypervolumes. 
  
You may find it useful to define a Jaccard-type fractional overlap between hv1 and hv2 as {hv_set@HVList$Intersection@Volume / hv_set@HVList$Union@Volume}.
  
If {check_memory} is true, instead returns a scalar with the expected number of pairwise comparisons.
}

\examples{
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)
hv2 = hypervolume(subset(iris, Species=="virginica")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)

hv_set = hypervolume_set(hv1, hv2, reduction_factor = 0.5, check_memory=FALSE)

# no overlap found
sapply(hv_set@HVList, function(x) {x@Volume})
}
