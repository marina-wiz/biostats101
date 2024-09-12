#' Calculate the Power and Sample Size for 2 Independent Proportions.
#'
#' This function calculates
#' 1. The power given independent proportions p1 and p2, the sample sizes for
#'    the two groups, the confidence level and the alternative hypothesis.
#' 2. The required sample size given independent proportions p1 and p2, the
#'    desired power, the confidence level and the alternative hypothesis.
#' 3. The required sample sizes n1 and n2 given the desired power, a sample
#'    size ratio n2/n1 (unbalanced designs), the confidence level and the
#'    alternative hypothesis.
#' 4. The sample size for the second group given independent proportions
#'    p1 and p2, the sample size for the first group (unbalanced designs), the
#'    desired power, the confidence level and the alternative hypothesis.
#' @param p1 Numeric, the proportion for the first group.
#' @param p2 Numeric, the proportion for the second group.
#' @param n1 Numeric, the sample size for the first group.
#' @param n2 Numeric, the sample size for the second group.
#' @param nratio Numeric, the sample size ratio (n2 / n1) for unbalanced
#' designs. Default is 1 when calculating sample sizes n1 and n2.
#' @param power Numeric, the desired power (1 - beta). Default is 0.8 when
#' calculating sample sizes n1 and n2 and when calculating n2 given n1.
#' @param conf.level Numeric, the confidence level (1 - alpha). Default is 0.95.
#' @param alternative Character, the type of alternative hypothesis.
#' Options are 'two.sided' (default) or 'one.sided'.
#' @param continuity Logical, indicating whether the continuity correction
#' should be applied. Default is FALSE.
#' @return A list with the following components:
#' - `n`: Total sample size (n1 + n2).
#' - `n1`: Sample size for the first group.
#' - `n2`: Sample size for the second group.
#' - `power`: The estimated power.
#' - `nratio`: The sample size ratio (n2 / n1), if applicable.
#' @references
#' Levin, B., & Chen, X. (1999). Is the one-half continuity correction used
#' once or twice to derive a well-known approximate sample size formula to
#' compare two independent binomial distributions?. *The American Statistician*,
#' 53(1), 62-66. https://doi.org/10.1080/00031305.1999.10474431.
#' Fleiss, J. L., Levin, B., & Paik, M. C. (2013). *Statistical methods for
#' rates and proportions*. John Wiley & Sons.
#' @examples
#' # Calculate the power for independent proportions given the sample sizes
#' power.2p.2n(p1 = 0.45, p2 = 0.6, n1 = 260, n2 = 130)
#'
#' # Calculate the sample size for independent proportions (default power = 0.8)
#' power.2p.2n(p1 = 0.45, p2 = 0.6)
#'
#' # Calculate n1 and n2 for independent proportions with ratio n2/n1
#' power.2p.2n(p1 = 0.44, p2 = 0.6, nratio = 2)
#'
#' # Calculate n2 given n1 for independent proportions
#' power.2p.2n(p1 = 0.44, p2 = 0.6, n1 = 108)
#' @export
power.2p.2n <- function(
    p1, p2, n1 = NULL, n2 = NULL, nratio = NULL, power = NULL,
    alternative = "two.sided", conf.level = 0.95, continuity = FALSE) {
  if (!is.null(power) && !is.null(n1) && !is.null(n2) && !is.null(nratio)) {
    stop("At least one of 'n1', 'n2', 'nratio', and 'power' must be NULL")
  }
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
  if (!is.null(conf.level)) {
    if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
      stop("'conf.level' must be numeric in [0, 1].")
    }
  }
  if (!is.null(power)) {
    if (!is.numeric(power) || conf.level <= 0 || conf.level >= 1) {
      stop("'power' must be numeric in [0, 1].")
    }
  }

  # Part 1. Estimate power given n1 & n2:
  if (!is.null(n1) && !is.null(n2) && is.null(nratio) && is.null(power)) {

    # Pooled proportion
    p_pool <- (n1 * p1 + n2 * p2) / (n1 + n2)

    # Standard deviations
    sigma_D <- sqrt((p1 * (1 - p1)) / n1 + (p2 * (1 - p2)) / n2)
    sigma_p <- sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))

    # Continuity correction term
    nratio = n2 / n1
    w1 = 1 / (1 + nratio)
    w2 = nratio / (1 + nratio)
    c <-  1 / (2 * (n1 + n2) * w1 * w2)

    # Z-scores for the confidence level
    alpha <- 1 - conf.level
    if (alternative == "two.sided") {
      z_alpha <- qnorm(1 - alpha/2)
      # z_beta <- qnorm(1 - alpha/2 / 2)
    } else if (alternative == "one.sided") {
      z_alpha <- qnorm(1 - alpha)
    } else {
      stop("Alternative must be either 'two.sided' or 'one.sided'.")
    }

    # Compute the power
    if (alternative == "one.sided" & p2 > p1) {
      power <- pnorm(((p2 - p1) - ifelse(continuity, c, 0)
                      - z_alpha * sigma_p) / sigma_D)
    } else if (alternative == "one.sided" & p1 < p2) {
      power <- pnorm((-(p2 - p1) - ifelse(continuity, c, 0)
                      - z_alpha * sigma_p) / sigma_D)
    } else {  # two-sided
      power <- pnorm(((p2 - p1) - ifelse(continuity, c, 0)
                      - z_alpha * sigma_p) / sigma_D) +
        pnorm((-(p2 - p1) - ifelse(continuity, c, 0)
               - z_alpha * sigma_p) / sigma_D)
    }

    # Return n, n1, n2, and power
    result <- list(n = n1 + n2, n1 = n1, n2 = n2, power = round(power, 4))
    return(result)
  }

  # Define two-sided power function (for parts 2 and 3)
  power_2p_2n_two_sided <- function(n1, n2) {
    z_alpha <- qnorm(1 - alpha/2)
    p_pool <- (n1 * p1 + n2 * p2) / (n1 + n2)
    sigma_D <- sqrt((p1 * (1 - p1)) / n1 + (p2 * (1 - p2)) / n2)
    sigma_p <- sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
    nratio = n2 / n1
    w1 = 1 / (1 + nratio)
    w2 = nratio / (1 + nratio)
    c <-  1 / (2 * (n1 + n2) * w1 * w2)
    alpha <- 1 - conf.level
    z_alpha <- qnorm(1 - alpha/2)
    power <- pnorm(((p2 - p1) - ifelse(continuity, c, 0)
                    - z_alpha * sigma_p) / sigma_D) +
      pnorm((-(p2 - p1) - ifelse(continuity, c, 0)
             - z_alpha * sigma_p) / sigma_D)
    return(power)
  }

  # Define one-sided power function (for parts 2 and 3)
  power_2p_2n_one_sided <- function(n1, n2) {
    alpha <- 1 - conf.level
    z_alpha <- qnorm(1 - alpha)
    p_pool <- (n1 * p1 + n2 * p2) / (n1 + n2)
    sigma_D <- sqrt((p1 * (1 - p1)) / n1 + (p2 * (1 - p2)) / n2)
    sigma_p <- sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
    nratio = n2 / n1
    w1 = 1 / (1 + nratio)
    w2 = nratio / (1 + nratio)
    c <-  1 / (2 * (n1 + n2) * w1 * w2)
    z_alpha <- qnorm(1 - alpha)
    sigma_D <- sqrt((p1 * (1 - p1)) / n1 + (p2 * (1 - p2)) / n2)
    sigma_p <- sqrt(p_pool * (1 - p_pool) * (1/n1 + 1/n2))
    if (p2 > p1) {
      power <- pnorm(((p2 - p1) - ifelse(continuity, c, 0)
                      - z_alpha * sigma_p) / sigma_D)
    } else if (p1 < p2) {
      power <- pnorm((-(p2 - p1) - ifelse(continuity, c, 0)
                      - z_alpha * sigma_p) / sigma_D)
    }
    return(power)
  }

  # Part 2. Estimate sample size given power and nratio:
  if (is.null(n1) && is.null(n2)) {
    if (is.null(power)) {
      power = 0.8
    }
    if (is.null(nratio)) {
      nratio = 1
    }
    if (!alternative %in% c("one.sided", "two.sided")) {
      stop("Alternative must be either 'two.sided' or 'one.sided'.")
    }

    if (alternative == "one.sided") {

      # Calculate z-values for confidence level and power
      alpha <- 1 - conf.level
      z_alpha <- qnorm(1 - alpha)
      z_beta <- qnorm(1 - power)

      # Pooled proportion
      p_pool <- (p1 + nratio * p2) / (1 + nratio)

      # Weights w1 and w2
      w1 <- 1 / (1 + nratio)
      w2 <- nratio / (1 + nratio)

      # Calculate the sample size using the formula provided
      n <- ((z_alpha * sqrt(p_pool * (1 - p_pool))
             - z_beta * sqrt(w2 * p1 * (1 - p1) + w1 * p2 * (1 - p2)))^2) /
        (w1 * w2 * (p2 - p1)^2)

      # If continuity correction is applied
      if (continuity) {
        n <- n / 4 * (1 + sqrt(1 + (2 / (n * w1 * w2 * abs(p2 - p1)))))^2
      }

      # Calculate n1 and n2 using the ratio
      n1 <- ceiling(n / (1 + nratio))
      n2 <- ceiling(nratio * n1)

      current_power <- power_2p_2n_one_sided(n1, n2)
    }

    if (alternative == "two.sided") {
      # Calculate z-values for confidence level and power
      alpha <- 1 - conf.level
      z_alpha <- qnorm(1 - alpha/2)
      z_beta <- qnorm(1 - power)

      # Pooled proportion
      p_pool <- (p1 + nratio * p2) / (1 + nratio)

      # Weights w1 and w2
      w1 <- 1 / (1 + nratio)
      w2 <- nratio / (1 + nratio)

      # Use one-sided sample size formula to get initial values for n1 and n2
      n <- ((z_alpha * sqrt(p_pool * (1 - p_pool))
             - z_beta * sqrt(w2 * p1 * (1 - p1) + w1 * p2 * (1 - p2)))^2) /
        (w1 * w2 * (p2 - p1)^2)

      # If continuity correction is applied
      if (continuity) {
        n <- n / 4 * (1 + sqrt(1 + (2 / (n * w1 * w2 * abs(p2 - p1)))))^2
      }

      # Calculate initial n1 and n2 using the ratio
      n1 <- ceiling(n / (1 + nratio))
      n2 <- ceiling(nratio * n1)

      # Calculate initial power
      current_power <- power_2p_2n_two_sided(n1, n2)

      # iteratively solve power equation

      # if current power is above desired power:
      if (current_power >= power) {
        n <- n - 1
        n1 <- ceiling(n / (1 + nratio))
        n2 <- ceiling(nratio * n1)
        current_power <- power_2p_2n_two_sided(n1, n2)
      }

      # if current power is below desired power:
      if (current_power < power) {
        n <- n + 1
        n1 <- ceiling(n / (1 + nratio))
        n2 <- ceiling(nratio * n2)
        current_power <- power_2p_2n_two_sided(n1, n2)
      }
    }
    # Return n, n1, n2, nratio, and power
    result <- list(n = n1 + n2, n1 = n1, n2 = n2, nratio = round(n2/n1, 4),
                   power = round(current_power, 4))
    return(result)
  }

  # 3. Estimate n2 given n1 (or n1 given n2)
  if (!is.null(n1) && is.null(n2) && is.null(nratio)) {

    if (is.null(power)) {
      power = 0.8
    }

    if (!alternative %in% c("one.sided", "two.sided")) {
      stop("Alternative must be either 'two.sided' or 'one.sided'.")
    }

    # initial value for n2 using one-sided equation and assuming nratio = 1.
    n2 = n1

    # Calculate z-values for confidence level and power
    alpha <- 1 - conf.level
    z_alpha <- qnorm(1 - alpha)
    z_beta <- qnorm(1 - power)

    # Pooled proportion and weights
    p_pool <- (p1 + p2) / 2
    w1 <- 1 / 2
    w2 <- 1 / 2

    # Calculate the sample size using the formula provided
    n <- ((z_alpha * sqrt(p_pool * (1 - p_pool))
           - z_beta * sqrt(w2 * p1 * (1 - p1) + w1 * p2 * (1 - p2)))^2) /
      (w1 * w2 * (p2 - p1)^2)

    # If continuity correction is applied
    if (continuity) {
      n <- n / 4 * (1 + sqrt(1 + (2 / (n * w1 * w2 * abs(p2 - p1)))))^2
    }

    # initial n2 value
    n2 = ceiling(n) - n1

    if (alternative == "two.sided") {

      # Calculate initial power
      current_power <- power_2p_2n_two_sided(n1, n2)

      # iteratively solve power equation

      # if current power is above desired power:
      while (current_power >= power) {
        n2 <- n2 - 1
        current_power <- power_2p_2n_two_sided(n1, n2)
      }

      # if current power is below desired power:
      while (current_power < power) {
        n2 <- n2 + 1
        current_power <- power_2p_2n_two_sided(n1, n2)
      }
    }

    if (alternative == "one.sided") {

      # Calculate initial power
      current_power <- power_2p_2n_one_sided(n1, n2)

      # iteratively solve power equation

      # if current power is above desired power:
      while (current_power >= power) {
        n2 <- n2 - 1
        current_power <- power_2p_2n_one_sided(n1, n2)
      }

      # if current power is below desired power:
      while (current_power < power) {
        n2 <- n2 + 1
        current_power <- power_2p_2n_one_sided(n1, n2)
      }
    }
    # Return n, n1, n2, and power
    result <- list(n = n1 + n2, n1 = n1, n2 = n2,
                   power = round(current_power, 4))
    return(result)
  }
  if (!is.null(n2) && is.null(n1) && is.null(nratio)) {

    if (is.null(power)) {
      power = 0.8
    }

    if (!alternative %in% c("one.sided", "two.sided")) {
      stop("Alternative must be either 'two.sided' or 'one.sided'.")
    }

    # initial value for n2 using one-sided equation and assuming nratio = 1.
    n1 = n2

    # Calculate z-values for confidence level and power
    alpha <- 1 - conf.level
    z_alpha <- qnorm(1 - alpha)
    z_beta <- qnorm(1 - power)

    # Pooled proportion and weights
    p_pool <- (p1 + p2) / 2
    w1 <- 1 / 2
    w2 <- 1 / 2

    # Calculate the sample size using the formula provided
    n <- ((z_alpha * sqrt(p_pool * (1 - p_pool))
           - z_beta * sqrt(w2 * p1 * (1 - p1) + w1 * p2 * (1 - p2)))^2) /
      (w1 * w2 * (p2 - p1)^2)

    # If continuity correction is applied
    if (continuity) {
      n <- n / 4 * (1 + sqrt(1 + (2 / (n * w1 * w2 * abs(p2 - p1)))))^2
    }

    # initial n2 value
    n1 = ceiling(n) - n2

    if (alternative == "two.sided") {

      # Calculate initial power
      current_power <- power_2p_2n_two_sided(n1, n2)

      # iteratively solve power equation

      # if current power is above desired power:
      while (current_power >= power) {
        n1 <- n1 - 1
        current_power <- power_2p_2n_two_sided(n1, n2)
      }

      # if current power is below desired power:
      while (current_power < power) {
        n1 <- n1 + 1
        current_power <- power_2p_2n_two_sided(n1, n2)
      }
    }

    if (alternative == "one.sided") {

      # Calculate initial power
      current_power <- power_2p_2n_one_sided(n1, n2)

      # iteratively solve power equation

      # if current power is above desired power:
      while (current_power >= power) {
        n1 <- n1 - 1
        current_power <- power_2p_2n_one_sided(n1, n2)
      }

      # if current power is below desired power:
      while (current_power < power) {
        n1 <- n1 + 1
        current_power <- power_2p_2n_one_sided(n1, n2)
      }
    }
    # Return n, n1, n2, and power
    result <- list(n = n1 + n2, n1 = n1, n2 = n2,
                   power = round(current_power, 4))
    return(result)
  }
}
