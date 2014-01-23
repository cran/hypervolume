setClass("Hypervolume", slots=c(
    Name="character",
    Data="matrix",
    Dimensionality="numeric",
    Volume="numeric",
    PointDensity="numeric",
    Bandwidth="numeric",
    RepsPerPoint="numeric",
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
  cat(sprintf("Hypervolume\n\tName: %s\n\tNr. of observations: %.0f\n\tDimensionality: %.0f\n\tVolume: %f\n\tBandwidth: %s\n\tQuantile desired: %f\n\tQuantile obtained: %f\n\tNr. of repetitions per point: %.0f\n\tNumber of random points: %.0f\n", 
              object@Name, ifelse(all(is.nan(object@Data)), NA, nrow(object@Data)), object@Dimensionality, object@Volume, paste(format(object@Bandwidth,digits=2),collapse=' '), object@QuantileThresholdDesired, object@QuantileThresholdObtained, object@RepsPerPoint, nrow(object@RandomUniformPointsThresholded)))

}

summary.HypervolumeList <- function(object, ...)
{
  for (i in 1:length(object@HVList))
  {
    summary(object@HVList[[i]])
  }
}

