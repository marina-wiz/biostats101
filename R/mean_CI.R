#' Construct the Confidence Interval for the Mean
#'
#' Constructs a confidence interval (CI) for the mean of a numeric vector using the
#' Student's t-distribution. The CI is calculated based on the specified degrees of
#' freedom, confidence level, and alternative hypothesis.
#'
#' @param data A numeric vector from which the mean and confidence interval
#' will be calculated.
#' @param conf.level A numeric value representing the confidence level for the
#' confidence interval estimation. Default is 0.95.
#' @param alternative The alternative hypothesis to be considered on the
#' confidence interval estimation.
#' Options are 'two.sided' (default), 'greater', or 'less'.
#' @return A named numeric vector with the mean and the lower and upper bounds
#' of the confidence interval.
#' @examples
#' # Example data
#' values = c(5.2, 4.8, 6.3, 6.1, 7.2, 3.5, 4.9, 2.2, 3.7, 3.5, 8.9)
#'
#' # Construct a 95% confidence interval for the mean
#' mean_CI(values)
#' @export
mean_CI <- function(data, conf.level = 0.95, alternative = "two.sided") {

  # Remove NA values
  data <- na.omit(data)

  # Sample size
  n <- length(data)

  # Sample mean
  sample_mean <- mean(data)

  # Sample standard deviation
  std_dev <- sd(data)

  # Calculate alpha
  alpha <- 1 - conf.level

  # Calculate confidence interval based on the specified type
  if (alternative == "two.sided") {
    lower_bound <-
      sample_mean - qt(p = 1 - alpha/2, df = n - 1) * std_dev / sqrt(n)
    upper_bound <-
      sample_mean + qt(p = 1 - alpha/2, df = n - 1) * std_dev / sqrt(n)
  } else if (alternative == "greater") {
    lower_bound <-
      sample_mean - qt(p = 1 - alpha, df = n - 1) * std_dev / sqrt(n)
    upper_bound <- Inf
  } else if (alternative == "less") {
    lower_bound <- -Inf
    upper_bound <-
      sample_mean + qt(p = 1 - alpha, df = n - 1) * std_dev / sqrt(n)
  } else {
    stop("Alternative must be one of 'two.sided', 'greater', or 'less'.")
  }

  # Create a named vector with mean and confidence interval bounds
  results <- c(sample_mean, lower_bound, upper_bound)
  names(results) <- c("mean", paste0("lwr.", 100 * conf.level),
                      paste0("upr.", 100 * conf.level))

  return(results)
}
