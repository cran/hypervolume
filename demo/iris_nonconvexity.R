# Load the built-in Fisher/Anderson dataset
data(iris)

# extract trait data for just Iris versicolor, then rescale all axes via z-transformation
data_versicolor <- scale(subset(iris, Species=='versicolor')[,1:4])

# generate a hypervolume for the floral trait data
hv_versicolor <- hypervolume(data_versicolor, bandwidth=0.5, name='versicolor')

# find non-convex features (run without error-checking first)
features_versicolor <- nonconvexfeatures(hv_versicolor, check_memory=FALSE, check_convexhull=FALSE)

# plot the non-convex features inferred for the largest inflation factor
plot(features_versicolor[[length(features_versicolor)]])