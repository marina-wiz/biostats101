#' Calculate the Power and Sample Size for Paired Proportions.
#'
#' This function calculates either the power given the sample size or the
#' sample size given the power for paired proportions p1 and p2.
#'
#' @param p1 Numeric, the proportion at the first occasion.
#' @param p2 Numeric, the proportion at the second occasion.
#' @param n Numeric, the sample size.
#' @param power Numeric, the desired power (1 - beta). Default is 0.8 when
#' calculating sample size.
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
#' power.paired.prop(p1 = 0.15, p2 = 0.1, power = 0.8)
#' @export
power.paired.prop <- function(
    p1, p2, n = NULL, power = NULL,
    conf.level = 0.95, alternative = "two.sided") {

  if (p1 == p2) {
    stop("'p1' and 'p2' proportions are equal; this is not allowed")
  }
  if (is.null(p1) || is.null(p2)) {
    stop("Both 'p1' and 'p2' must be provided.")
  }

  if (!is.numeric(p1) || p1 <= 0 || p1 >= 1) {
    stop("'p1' must be numeric in [0, 1].")
  }
  if (!is.numeric(p2) || p2 <= 0 || p2 >= 1) {
    stop("'p2' must be numeric in [0, 1].")
  }

  if (!is.null(n) && !is.null(power)) {
    stop("Only one of 'n' or 'power' should be provided.")
  }

  if (!alternative %in% c("two.sided", "one.sided")) {
    stop("Alternative must be one of 'two.sided' or 'one.sided'.")
  }

  # define terms
  alpha <- 1 - conf.level
  pdiff <- p2 - p1
  pdisc <- p1 + p2

  if (!is.null(n)) { # If sample size is provided, compute power

    # Power computation
    if (alternative == "two.sided") {
      z_alpha <- qnorm(1 - alpha / 2)
      term1 <- (pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
        sqrt(pdisc - pdiff^2)
      term2 <- (-pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
        sqrt(pdisc - pdiff^2)
      power <- pnorm(term1) + pnorm(term2)
    } else { # one.sided
      z_alpha <- qnorm(1 - alpha)
      if (p1 > p2) {
        term <- (-pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
          sqrt(pdisc - pdiff^2)
      } else {
        term <- (pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
          sqrt(pdisc - pdiff^2)
      }
      power <- pnorm(term)
    }
    # Results for power calculation
    result_power <- list(sample_size = n, power = round(power, 4),
                         conf.level = conf.level, alternative = alternative)
    return(result_power)
  } else { # If sample size is not provided, compute sample size

    # If sample size and power are not provided, default power = 0.8
    if (is.null(power)) {
      power = 0.8
    }

    # Sample size computation
    if (alternative == "one.sided") {

      # Define power function for one-sided test
      power_one_sided <- function(p1, p2, n, alpha) {
        z_alpha <- qnorm(1 - alpha)
        if (p1 > p2) {
          term <- (-pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
            sqrt(pdisc - pdiff^2)
        } else {
          term <- (pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
            sqrt(pdisc - pdiff^2)
        }
        power <- pnorm(term)
        return(power)
      }

      # Compute the sample size for a one-sided test
      z_alpha <- qnorm(1 - alpha)
      z_beta <- qnorm(power)
      num <- (z_alpha * sqrt(pdisc) + z_beta * sqrt(pdisc - pdiff^2))^2
      denom <- pdiff^2
      n <- ceiling(num / denom)

      # Compute power for output
      current_power <- power_one_sided(p1, p2, n, alpha)

      # Results for one-sided sample size calculation
      result_n_one_sided <- list(sample_size = n,
                                 power = round(current_power, 4),
                                 conf.level = conf.level,
                                 alternative = alternative)
      return(result_n_one_sided)

    } else { # Compute the sample size for a two-sided test

      # Define power function for two-sided test
      power_two_sided <- function(p1, p2, n, alpha) {

        z_alpha <- qnorm(1 - alpha / 2)
        term1 <- (pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
          sqrt(pdisc - pdiff^2)
        term2 <- (-pdiff * sqrt(n) - z_alpha * sqrt(pdisc)) /
          sqrt(pdisc - pdiff^2)
        power <- pnorm(term1) + pnorm(term2)

        return(power)
      }

      # Calculate initial sample size using one-sided formula with alpha / 2
      z_alpha <- qnorm(1 - alpha / 2)
      z_beta <- qnorm(power)
      num <- (z_alpha * sqrt(pdisc) + z_beta * sqrt(pdisc - pdiff^2))^2
      denom <- pdiff^2
      n <- ceiling(num / denom)

      # Current power
      current_power <- power_two_sided(p1, p2, n, alpha)

      # If current power is above desired power:
      while (current_power >= power) {
        n <- n - 1
        current_power <- power_two_sided(p1, p2, n, alpha)
      }

      # If current power is below desired power:
      while (current_power < power) {
        n <- n + 1
        current_power <- power_two_sided(p1, p2, n, alpha)
      }

      # Results for sample size calculation (two-sided)
      result_n_two_sided <- list(sample_size = n,
                                 power = round(current_power, 4),
                                 conf.level = conf.level, alternative = alternative)
      return(result_n_two_sided)
    }
  }
}
