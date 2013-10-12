hypervolume_set <- function(hv1, hv2, reduction_factor=1, verbose=T)
{
  # determine dataset sizes and dimensionality
  np1 = nrow(hv1@RandomUniformPointsThresholded)
  np2 = nrow(hv2@RandomUniformPointsThresholded)
  
  hv1_point_density = np1 / hv1@Volume
  hv2_point_density = np2 / hv2@Volume
  
  dimhv1 = ncol(hv1@RandomUniformPointsThresholded)
  dimhv2 = ncol(hv2@RandomUniformPointsThresholded)
  
  if (dimhv1 != dimhv2)
  {
    stop('Dimensionality of hypervolumes is not the same.')
  }
  
  if (reduction_factor <= 0 | reduction_factor > 1)
  {
    stop('Reduction factor is not in (0,1].')
  }
  
  # define the consistent dimensionality
  dim = dimhv1
  
  # sample both hypervolumes down to the minimum point density
  mindensity = min(c(hv1_point_density, hv2_point_density))
  
  numpointstokeep_hv1 = floor(reduction_factor * mindensity * hv1@Volume)
  numpointstokeep_hv2 = floor(reduction_factor * mindensity * hv2@Volume)
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d points in hv1 and %d points in hv2.\n', numpointstokeep_hv1, numpointstokeep_hv2))
  }
  
  if (reduction_factor < 1)
  {
    hv1_points_ss = hv1@RandomUniformPointsThresholded[sample(1:np1,size=numpointstokeep_hv1),]
    hv2_points_ss = hv2@RandomUniformPointsThresholded[sample(1:np2,size=numpointstokeep_hv2),]
  }
  else
  {
    hv1_points_ss = hv1@RandomUniformPointsThresholded
    hv2_points_ss = hv2@RandomUniformPointsThresholded
  }

  
  point_density = nrow(hv1_points_ss) / hv1@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dim) 
  
  
  # figure out which points are 'close enough' to other points
  # (within a n-ball of the critical distance)
  if (verbose==TRUE)
  {
    cat('Beginning ball queries... \n')
  }
  p2_in_1_all = evalfspherical(hv1_points_ss, cutoff_dist, hv2_points_ss,verbose=verbose)
  p1_in_2_all = evalfspherical(hv2_points_ss, cutoff_dist, hv1_points_ss,verbose=verbose)
  if (verbose==TRUE)
  {
    cat('Finished ball queries. \n')
  }
  
  # subset to retain only those 'close enough' points
  p2_in_1 = hv2_points_ss[p2_in_1_all > 0,]
  p1_in_2 = hv1_points_ss[p1_in_2_all > 0,]
  
  # the final volume is proportional to the fraction 
  # of points in hv1 in hv2, and vice versa
  v1 = nrow(p1_in_2) / nrow(hv1_points_ss) * hv1@Volume
  v2 = nrow(p2_in_1) / nrow(hv2_points_ss) * hv2@Volume
  
  # take the lower estimate as a conservative estimate
  final_volume_intersection = min(c(v1,v2))
  
  # create the intersection point cloud by merging both sets of sampled points.
  final_points_intersection = unique(rbind(p1_in_2, p2_in_1))
  final_density_intersection = nrow(final_points_intersection) / final_volume_intersection
  
  # now find the union point cloud
  p1_not_in_2 = hv1_points_ss[p1_in_2_all == 0,] # the points only in the first hypervolume
  p2_not_in_1 = hv2_points_ss[p2_in_1_all == 0,] # the points only in the second hypervolume
  
  num_points_to_sample_in_intersection = floor(point_density * final_volume_intersection) # choose the right number of points to keep the point density constant
  p_in_1_and_2 = final_points_intersection[sample(1:nrow(final_points_intersection), size=num_points_to_sample_in_intersection),] # randomly sample the intersection to grab
  
  final_volume_union = hv1@Volume + hv2@Volume - final_volume_intersection # union is sum minus intersection 
  final_points_union = unique(rbind(p1_not_in_2, p2_not_in_1, p_in_1_and_2)) 
  final_density_union = nrow(final_points_union) / final_volume_union
  
  # calculate the unique components for hv1 
  # (these will occasionally be too permissive and show some outlying points)
  final_volume_unique_hv1 = hv1@Volume - final_volume_intersection
  final_points_in_unique_1 = unique(p1_not_in_2)
  final_density_unique_1 = nrow(final_points_in_unique_1) / final_volume_unique_hv1
  
  # calculate the unique components for hv2
  final_volume_unique_hv2 = hv2@Volume - final_volume_intersection
  final_points_in_unique_2 = unique(p2_not_in_1)
  final_density_unique_2 = nrow(final_points_in_unique_2) / final_volume_unique_hv2
  
  # clean up memory
  gc()
  
  # prepare final hypervolumes
  result_intersection = new("Hypervolume")
  result_intersection@Name = sprintf("Intersection of (%s, %s)", hv1@Name, hv2@Name)
  result_intersection@Data = matrix(NaN,nrow=1,ncol=dim)
  result_intersection@Dimensionality = dim
  result_intersection@Volume = final_volume_intersection
  result_intersection@PointDensity = final_density_intersection
  result_intersection@Bandwidth = rep(NaN,dim)
  result_intersection@QuantileThresholdDesired = 0
  result_intersection@QuantileThresholdObtained = 0
  result_intersection@RandomUniformPointsThresholded = as.matrix(final_points_intersection); dimnames(result_intersection@RandomUniformPointsThresholded) = dimnames(hv1@RandomUniformPointsThresholded);
  result_intersection@ProbabilityDensityAtRandomUniformPoints = rep(1,nrow(result_intersection@RandomUniformPointsThresholded))
 
  result_union = new("Hypervolume")
  result_union@Name = sprintf("Union of (%s, %s)", hv1@Name, hv2@Name)
  result_union@Data = matrix(NaN,nrow=1,ncol=dim)
  result_union@Dimensionality = dim
  result_union@Volume = final_volume_union
  result_union@PointDensity = final_density_union
  result_union@Bandwidth = rep(NaN,dim)
  result_union@QuantileThresholdDesired = 0
  result_union@QuantileThresholdObtained = 0
  result_union@RandomUniformPointsThresholded = as.matrix(final_points_union); dimnames(result_union@RandomUniformPointsThresholded) = dimnames(hv1@RandomUniformPointsThresholded);
  result_union@ProbabilityDensityAtRandomUniformPoints = rep(1,nrow(result_union@RandomUniformPointsThresholded))
  
  result_unique_hv1 = new("Hypervolume")
  result_unique_hv1@Name = sprintf("Unique component of (%s) relative to (%s)", hv1@Name, hv2@Name)
  result_unique_hv1@Data = matrix(NaN,nrow=1,ncol=dim)
  result_unique_hv1@Dimensionality = dim
  result_unique_hv1@Volume = final_volume_unique_hv1
  result_unique_hv1@PointDensity = final_density_unique_1
  result_unique_hv1@Bandwidth = rep(NaN,dim)
  result_unique_hv1@QuantileThresholdDesired = 0
  result_unique_hv1@QuantileThresholdObtained = 0
  result_unique_hv1@RandomUniformPointsThresholded = as.matrix(final_points_in_unique_1); dimnames(result_unique_hv1@RandomUniformPointsThresholded) = dimnames(hv1@RandomUniformPointsThresholded);
  result_unique_hv1@ProbabilityDensityAtRandomUniformPoints = rep(1,nrow(result_unique_hv1@RandomUniformPointsThresholded))

  result_unique_hv2 = new("Hypervolume")
  result_unique_hv2@Name = sprintf("Unique component of (%s) relative to (%s)", hv2@Name, hv1@Name)
  result_unique_hv2@Data = matrix(NaN,nrow=1,ncol=dim)
  result_unique_hv2@Dimensionality = dim
  result_unique_hv2@Volume = final_volume_unique_hv2
  result_unique_hv2@PointDensity = final_density_unique_2
  result_unique_hv2@Bandwidth = rep(NaN,dim)
  result_unique_hv2@QuantileThresholdDesired = 0
  result_unique_hv2@QuantileThresholdObtained = 0
  result_unique_hv2@RandomUniformPointsThresholded = as.matrix(final_points_in_unique_2); dimnames(result_unique_hv2@RandomUniformPointsThresholded) = dimnames(hv2@RandomUniformPointsThresholded);
  result_unique_hv2@ProbabilityDensityAtRandomUniformPoints = rep(1,nrow(result_unique_hv2@RandomUniformPointsThresholded))
  
  # assemble final results into a list
  result = new("HypervolumeList")
  result@HVList = list(Intersection = result_intersection, Union = result_union, Unique_1 = result_unique_hv1, Unique_2 = result_unique_hv2)
  
  return(result)
}