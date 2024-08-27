#' Create a Plot with Linear Regression and Confidence/Prediction Intervals
#'
#' Generates a `ggplot2` plot that includes the observed data, a simple linear
#' regression line, and both confidence and prediction interval bands.
#'
#' @name lm_plot
#' @param lm_model A simple linear regression model object created by `lm()`.
#' @param conf.level  A numeric value indicating the confidence level for the
#' confidence interval bands and prediction interval bands. Default is 0.95.
#' @param install_packages Logical, indicating whether required packages should
#' be automatically installed if not already available. Default is TRUE.
#' @return A `ggplot2` plot object featuring the observations,
#' the simple linear regression line, the confidence interval bands, and
#' the prediction interval bands.
#'
#' @importFrom stats na.omit pnorm predict qnorm qt sd setNames
#' @importFrom utils install.packages
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   x = rnorm(100, mean = 50, sd = 10),
#'   y = 3 + 0.5 * rnorm(100, mean = 50, sd = 10) + rnorm(100)
#' )
#'
#' # Run a regression model
#' my_model <- lm(y~x, data)
#'
#' # Create plot with regression line, confidence limits, and prediction limits
#' lm_plot(my_model)
#'
#' # Customize plot labels
#' lm_plot(my_model) + xlab("Your x-axis label") + ylab("Your y-axis label")
#' @export
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("fit", "lwr.PI", "upr.PI", "lwr.CI", "upr.CI",
                           "value", "linetype", "type"))
}
lm_plot <- function(lm_model, conf.level = 0.95, install_packages = TRUE) {

  # Define the packages needed
  required_packages <- c("dplyr", "tidyr", "ggplot2")

  # Function to install and load packages
  if (install_packages) {
    install_and_load <- function(pkg) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }
    }
    lapply(required_packages, install_and_load)
  }

  # Generate confidence and prediction intervals
  model_CI <- as.data.frame(predict(lm_model, interval = "confidence", level = conf.level))
  model_PI <- as.data.frame(suppressWarnings(
    predict(lm_model, interval = "prediction", level = conf.level)))

  # Create new data and pivot longer for the plot
  new_data <- as.data.frame(list(
    lm_model$model, fit = model_CI$fit,
    lwr.CI = model_CI$lwr, upr.CI = model_CI$upr,
    lwr.PI = model_PI$lwr, upr.PI = model_PI$upr))
  new_data_long <- new_data %>%
    pivot_longer(cols = c(fit, lwr.PI, upr.PI, lwr.CI, upr.CI),
                 names_to = "type", values_to = "value")

  # Generate labels based on the confidence level
  conf_label <- paste0(100 * conf.level, "% confidence limits")
  pred_label <- paste0(100 * conf.level, "% prediction limits")

  # Set line types based on the type of interval
  new_data_long$linetype <- factor(case_when(
    new_data_long$type == "fit" ~ "Fit",
    new_data_long$type %in% c("lwr.PI", "upr.PI") ~ pred_label,
    new_data_long$type %in% c("lwr.CI", "upr.CI") ~ conf_label
  ))

  # Get variable names
  y_var <- names(lm_model$model)[1]
  x_var <- names(lm_model$model)[2]

  # Plot the data
  plot <- ggplot(new_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(size = 1) +
    geom_line(data = new_data_long, aes(
      y = value, linetype = linetype, group = type)) +
    theme_bw() +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_linetype_manual(values = setNames(c(
      "solid", "dotted", "dashed"),
      c("Fit", pred_label, conf_label)
    ), guide = guide_legend(title = NULL))

  # Return the plot
  return(plot)
}
