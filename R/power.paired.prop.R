#' Calculate the Power and Sample Size for Paired Proportions.
#'
#' This function calculates either the power given the sample size or the
#' sample size given the power for paired proportions p1 and p2.
#'
#' @param p1 Numeric, the proportion for the first group
#' @param p2 Numeric, the proportion for the second group.
#' @param n Numeric, the sample size. If `n` is provided, the function will
#' calculate the power.
#' @param power Numeric, the desired power (1 - beta). If `power` is provided,
#' the function will calculate the required sample size.
#' @param conf.level Numeric, the confidence level (1 - alpha). Default is 0.95.
#' @param alternative Character, the type of alternative hypothesis.
#' Options are 'two.sided' (default) or 'one.sided'.
#' @return A list containing the sample size, power, confidence level, and
#' alternative hypothesis.
#' @references
#' McNemar, Q. (1947). Note on the sampling error of the difference between
#' correlated proportions or percentages. *Psychometrika*, 12(2), 153-157.
#' https://doi.org/10.1007/BF02295996.
#' Connor, R. J. (1987). Sample size for testing differences in proportions for
#' the paired-sample design. *Biometrics*, 207-211.
#' https://doi.org/10.2307/2531961.
#' @examples
#' # Calculate the power given the sample size for paired proportions
#' power.paired.prop(p1 = 0.1, p2 = 0.15, n = 900)
#'
#' # Calculate the sample size given the power for paired proportions
#' power.paired.prop(p1 = 0.15, p2 = 0.1, n = 900, alternative='one.sided')
#' @export
power.paired.prop <- function(
    p1, p2, n = NULL, power = NULL,
    conf.level = 0.95, alternative = "two.sided") {

  if (is.null(n) && is.null(power)) {
    stop("Either sample size 'n' or 'power' must be provided.")
  }

  if (!is.null(n) && !is.null(power)) {
    stop("Only one of 'n' or 'power' should be provided.")
  }

  if (!alternative %in% c("two.sided", "one.sided")) {
    stop("Alternative must be one of 'two.sided' or 'one.sided'.")
  }

  # Function to compute power given sample size
  compute_power <- function(n, p1, p2, alpha, alternative) {

    alpha <- 1 - conf.level
    pdiff <- p2 - p1
    pdisc <- p1 + p2

    if (alternative == "two.sided") {
      z_alpha <- qnorm(1 - alpha / 2)
      term1 <- (pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) / sqrt(pdisc - pdiff^2)
      term2 <- (-pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) / sqrt(pdisc - pdiff^2)
      power <- pnorm(term1) + pnorm(term2)
    } else {
      z_alpha <- qnorm(1 - alpha)
      if (p1 > p2) {
        term <- (-pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) / sqrt(pdisc - pdiff^2)
      } else {
        term <- (pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) / sqrt(pdisc - pdiff^2)
      }
      power <- pnorm(term)
    }

    return(power)
  }

  # Function to compute sample size given power
  compute_sample_size <- function(power, p1, p2, alpha, alternative) {
    pdiff <- p2 - p1
    pdisc <- p1 + p2

    if (alternative == "two.sided") {
      z_alpha <- qnorm(1 - alpha / 2)
      z_beta <- qnorm(power)
      num <- (z_alpha * sqrt(pdisc) + z_beta * sqrt(pdisc - pdiff^2))^2
    } else {
      z_alpha <- qnorm(1 - alpha)
      z_beta <- qnorm(power)
      num <- (z_alpha * sqrt(pdisc) + z_beta * sqrt(pdisc - pdiff^2))^2
    }

    denom <- pdiff^2
    n <- num / denom
    return(ceiling(n))
  }

  if (!is.null(n)) {
    # Compute power given sample size
    result_power <- compute_power(n, p1, p2, alpha, alternative)
    return(list(sample_size = n, power = result_power, conf.level = conf.level,
                alternative = alternative))
  } else {
    # Compute sample size given power
    result_n <- compute_sample_size(power, p1, p2, alpha, alternative)
    return(list(sample_size = result_n, power = power, conf.level = conf.level,
                alternative = alternative))
  }
}
