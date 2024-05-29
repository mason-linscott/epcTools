#' Environmental spline function
#'
#' Given a dataframe of environment values at a given time, return a vector of length b that corresponds to environmental values at a time slice. This can be done through a loess smoothed function or from taking the average of values at a given time interval (binning).
#'
#' @param env A dataframe of environment data and time
#'
#' @param b Number of time slices
#'
#' @param R Resolution of spline
#'
#' @return A numeric likelihood value
#'
#' @export
getSliceEnv <- function(env,b){
    rootfirst_env<-env[order(-env[,1]), ]
    spline_fun<-loess(rootfirst_env[,2] ~ rootfirst_env[,1])
    env_range<-range(env[,1])
    spline_data<-as.data.frame(predict(spline_fun,seq(env_range[1],env_range[2],env_range[2]/(10*b))))
    n=nrow(spline_data)/b
    colnames(spline_data)<-"slice_averages"
    sliced_data<-aggregate(spline_data, list(rep(1:(nrow(spline_data) %/% n), each = n, len = nrow(spline_data))), mean)[-1];
    sliced_data<-rev(sliced_data[,1])
    return(sliced_data)

}

getBinEnv <- function(env,b){
  # Sort the data frame by time
  rootfirst_env<-env[order(-env[,1]), ]

  # Determine the range of the total time period
  time_range <- range(rootfirst_env[,1], na.rm = TRUE)

  # Calculate the width of each time slice
  slice_width <- abs(diff(time_range)) / b

  # Determine the breakpoints for time slices
  breakpoints <- seq(time_range[1], time_range[2], by = slice_width)

  # Initialize an empty vector to store the mean values
  mean_values <- numeric(b)

  # Calculate mean for each time slice
  for (i in 1:b) {
    # Extract the indices of rows within the current time slice
    slice_indices <- which(rootfirst_env[,1] >= breakpoints[i] & rootfirst_env[,1] <= breakpoints[i + 1])

    # Extract the environmental values for the current time slice
    slice_data <- rootfirst_env[,2][slice_indices]

    # Calculate the mean for the current time slice
    mean_values[i] <- mean(slice_data, na.rm = TRUE)
  }
  mean_values<-rev(mean_values)
  return(mean_values)

}
