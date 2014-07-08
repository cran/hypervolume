\name{hypervolume-package}
\alias{hypervolume-package}
\docType{package}
\title{
N-dimensional hypervolume operations
}
\description{
Estimates the shape and volume of high-dimensional objects and performs set operations: intersection / overlap, union, unique components, and inclusion test. Can measure the n-dimensional ecological hypervolume and perform species distribution modeling.
}

\details{
Frequently asked questions

{What quantile threshold should I use?} 

A value of 0.0 will ensure that the final hypervolume encloses all observed points. A value of 0.1 is more robust to outliers.

{How do I know the hypervolume algorithm has converged on the correct answer?} 

You can generate multiple hypervolumes, then plot computed volume against the number of Monte Carlo samples. You should choose a large enough number of Monte Carlo samples to constrain variation in volume to within your desired tolerance.

{What fraction of points should I use for set operations? If I use all the points the algorithm is too slow.} 

The set operations scale approximately quadratically with the number of points, so choosing a smaller point fraction will improve runtimes quadratically. However, using fewer points decreases the accuracy of the computation and may result in some points being identified as `out' when they should be `in', and could potentially bias estimates of computed volumes. Choose a large enough point fraction that the computed volume converges to a value within your desired tolerance.

{How can I export hypervolume data to another program for plotting?} 

Each hypervolume object stores a slot called \code{RandomUniformPointsThresholded}. Points in this field are a random uniform sample of the hypervolume, so you can export these points for graphing or inference in other programs.

{Why does the hypervolume look blocky, and not smooth around the edges?} 

The algorithms use a hyperbox kernel for computational efficiency. For small sample sizes or widely spaced points, the kernel bandwidth may be large relative to the spread of points. As a result, the hypervolume will appear to have jagged edges. This effect cannot be avoided in the current implementation.

{I have a multi-core machine. How can I speed up the algorithms?} 

You can run multiple hypervolume commands in parallel using several possible packages. See http://cran.r-project.org/web/views/HighPerformanceComputing.html for details.

{I don't understand what the `uniformly random' output means.} 

The hypervolume algorithm proceeds by 1) computing a kernel density estimate, 2) thresholding this estimate to include some quantile of the total probability density, 3) defining the volume as all points whose density value exceeds this threshold, and 4) sampling random points from within this volume, such that the point density of the final object is constant within the volume and zero outside the volume. This is the stochastic equivalent of a binary definition of an object (points are either in or out) and is necessary to do set operations. Thus the object is completely defined by the `edge' of these points and no further contouring should be done. 

{Why do some parts of the hypervolume appear to have higher point density?} 

This is because each pairplot is a projection of a high-dimensional object to two dimensions - so `thicker' regions of the object will effectively show more points.

{Why do I not get the same answer if I run the same code repeatedly?} 

The algorithms are stochastic and depend on the state of the random number generator. If results are unreliable, increase the number of Monte Carlo samples. Alternative you can make results repeatable by fixing the random number generator seed in your code, e.g. \code{set.seed(3)}.



}

\author{
Benjamin Blonder <bblonder@gmail.com>
}
\references{
Blonder, B., Lamanna, C., Violle, C. and Enquist, B. J. (2014), The n-dimensional hypervolume. Global Ecology and Biogeography, 23: 595-609. doi: 10.1111/geb.12146
}