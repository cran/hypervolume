\name{nonconvexfeatures}
\alias{nonconvexfeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-convex feature detection (dents and holes)
}
\description{
Infers the presence of non-convex features (dents and holes) within a hypervolume. The algorithm's approach is to 1) compute a convex hull that encloses the data points, representing the maximal hypervolume and simplest linear inequality constraints that could describe the data; 2) sample a random set of points from this hull, 3) compute a thresholded kernel density estimate of the data points, representing the multidimensional hypervolume that is the niche;  4) sample a random set of points from this hypervolume; 4) determine the subset of random points that are within the convex hull and are not within the hypervolume; 5) identify non-convex features by determining the connectivity of these random points; and 6) segment these random points into holes or dents by determining if they touch the convex boundary. 



}
\usage{
nonconvexfeatures(hv, 
  numconvexhullpoints = NULL, memory_reduction_factor = 1,
  check_convexhull = TRUE, check_memory = TRUE, 
  minnumpoints = 5, npointsdesired = 1000, 
  inflationfactors = seq(1, 2, by = 0.2), 
  plotdiag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
The hypervolume to be analyzed. Must have a non-empty \code{Data} slot, i.e. cannot be the output of set operations. The convex hull must be drawn around data points rather than uniformly random points. The kernel density estimate is not smooth due to the underlying hyperbox kernel. Thus delineating a convex hull around the hypervolume's random points necessarily will cause the algorithm to detect non-convex features where there are jagged transitions.  
}
  \item{numconvexhullpoints}{
Number of points in \code{hv@Data} to use to delineate the convex hull. Smaller numbers will produce exponentially faster results but also much less accurate delineation of the convex hull and potential omission of non-convex features. The algorithm does preferentially sample points far away from the data centroid to minimize this effect. Leaving this parameter as \code{NULL} results in the algorithm using all data points to delineate the convex hull.
}
  \item{memory_reduction_factor}{
A number between 0 and 1 indicating the fraction of random points in the convex hull and hypervolume to use for set operations. Larger values are more accurate but use more computational resources.
}
  \item{check_convexhull}{
If \code{TRUE}, calculates convex hull and estimates memory usage necessary to perform hit-and-run sampling for generating uniformly random points, then exits. If \code{FALSE}, prints resource usage and continues algorithm. It is very useful for preventing crashes to check the estimated memory usage on large or high dimensional datasets before running the full algorithm.
}
  \item{check_memory}{
If \code{TRUE}, estimates the memory usage required to perform set operations between the convex hull and hypervolume, then exits. If \code{FALSE}, prints resource usage and continues algorithm. It is very useful for preventing crashes to check the estimated memory usage on large or high dimensional datasets before running the full algorithm.

}
  \item{minnumpoints}{
Number of random points in the smallest non-convex feature to be retained by the algorithm. Very small non-convex features are assumed to be false positives and pruned from the algorithm's output. 
}
  \item{npointsdesired}{
Number of random points to sample from the convex hull. Larger values are more accurate.
}
  \item{inflationfactors}{
A vector of numbers representing multiplicative factors used when segmenting non-convex features. Features separated by at least a critical distance multiplied by each inflation factor are inferred to be separate objects. Larger values produce fewer but larger features. The algorithm returns a \code{HypervolumeList} of non-convex features for each value in \code{inflationfactors}, so that the user can see how the choice of critical distance affects inference of the number and size of non-convex features.  A value of approximately two seem to work well in practice.
}
  \item{plotdiag}{
Generate diagnostic plot of main steps in algorithm. If true, will create a file called 'diag.pdf'. Only works in n=2 dimensions.
}
}
\details{

The algorithm's computational complexity can scale exponentially with the number of data points in the hypervolume as well as the dimensionality. For this reason, the algorithm prints extensive diagnostic output describing the resources allocated at each step. The algorithm also terminates by default before computationally intensive steps to warn the user of the expected resources needed; this behavior can be disabled by setting the appropriate flags (\code{check_convexhull} and \code{check_memory} to false).
}
\value{
A list of \code{HypervolumeList} objects, suitable for plotting. Each item in the list corresponds to the non-convex features inferred for a given value of \code{inflationfactors} and is named accordingly. Each HypervolumeList has a \code{HVList} slot. In this slot, the first list item is the original hypervolume. All subsequent list items, if any, are individual non-convex features. They are named 'Dent' if a dent and 'Hole' if a hole.

Typically the user will want to access the last item in the list, as segmentation of non-convex features tends to be most accurate for larger inflation factors.
}
\references{
Blonder, B. Is niche geometry non-convex, and how would we know? (in preparation)  
}
\author{
Benjamin Blonder
}

\seealso{
\code{demo(iris_nonconvexity)}
}
\examples{
# generate annulus data
data_annulus <- data.frame(matrix(data=runif(4000),ncol=2))
names(data_annulus) <- c("x","y")
data_annulus  <- subset(data_annulus, 
sqrt((x-0.5)^2+(y-0.5)^2) > 0.4 & sqrt((x-0.5)^2+(y-0.5)^2) < 0.5)

# MAKE HYPERVOLUME (low reps for fast execution)
hv_annulus <- hypervolume(data_annulus,bandwidth=0.1,name='annulus',reps=500)

# DETECT FEATURES (low npoints for fast execution)
features_annulus <- nonconvexfeatures(hv_annulus, 
check_memory = FALSE, check_convexhull=FALSE, 
inflationfactors = 2, minnumpoints=10)

# ANALYZE FOR THE LARGEST INFLATION FACTOR
features_annulus_largest <- features_annulus[[length(features_annulus)]]

# PLOT RESULTS
plot(features_annulus_largest)

}