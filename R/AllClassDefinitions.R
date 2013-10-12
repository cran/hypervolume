setClass("Hypervolume", slots=c(
    Name="character",
    Data="matrix",
    Dimensionality="numeric",
    Volume="numeric",
    PointDensity="numeric",
    Bandwidth="numeric",
    QuantileThresholdDesired="numeric",
    QuantileThresholdObtained="numeric",
    RandomUniformPointsThresholded="matrix",
    ProbabilityDensityAtRandomUniformPoints="numeric"
    ))

setClass("HypervolumeList", slots=c(
    HVList="list"
  ))

summary.Hypervolume <- function(object, ...)
{
  cat(sprintf("Hypervolume\n\tName: %s\n\tNr. of observations: %d\n\tDimensionality: %d\n\tVolume: %f\n\tBandwidth: %s\n\tQuantile desired: %f\n\tQuantile obtained: %f\n\tNumber of random points: %d\n", 
              object@Name, ifelse(all(is.nan(object@Data)), NA, nrow(object@Data)), object@Dimensionality, object@Volume, paste(format(object@Bandwidth,digits=2),collapse=' '), object@QuantileThresholdDesired, object@QuantileThresholdObtained, nrow(object@RandomUniformPointsThresholded)))
  
}

summary.HypervolumeList <- function(object, ...)
{
  sapply(object@HVList, summary, ...)
}

