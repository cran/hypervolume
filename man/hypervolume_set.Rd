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
hypervolume_set(hv1, hv2, reduction_factor = 1, verbose = T)
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
}
\details{
Uses the inclusion test approach to identify points in the first hypervolume that are or are not within the second hypervolume and vice-versa. The intersection is the points in both hypervolumes, the union those in either hypervolume, and the unique components the points in one hypervolume but not the other.
}
\value{
A HypervolumeList object, with four items in its HVList slot:

  \item{Intersection }{The intersection of hv1 and hv2}
  \item{Union }{The union of hv1 and hv2}
  \item{Unique_1 }{The unique component of hv1 relative to hv2}
  \item{Unique_2 }{The unique component of hv2 relative to hv1}
  
  Note that the output hypervolumes will have lower random point densities than the input hypervolumes.
}

\examples{
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)
hv2 = hypervolume(subset(iris, Species=="virginica")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)

hv_set = hypervolume_set(hv1, hv2, reduction_factor = 0.5)

# no overlap found
sapply(hv_set@HVList, function(x) {x@Volume})
}
