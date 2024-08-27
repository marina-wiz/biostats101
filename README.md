# biostats101

## Description

This R package was designed to help beginners in biostatistics get started with ease. The package offers a set of user-friendly functions that fill the gaps in existing tools, making it easier for newcomers to perform essential biostatistical analyses without needing advanced programming skills. 

### Functions Overview

- **`mean_CI()`** Construct confidence intervals for the mean of a given variable. You can specify the confidence level and the alternative hypothesis.
- **`power.paired.prop()`** Calculate the power or sample size for paired proportions. You need to specify proportions $p_1$ and $p_2$, and either the power or the sample size. You can also specify the confidence level and the alternative hypothesis.
- **`power.2p.2n()`** Calculate the power or sample size(s) for independent proportions, for both balanced and unbalanced designs. You need to specify proportions $p_1$ and $p_2$. Additionally, you can specify:
	- The confidence level and alternative hypothesis.
	- The sample sizes $n_1$ and $n_2$ if you want to calculate the power. 
	- The sample size $n_1$ to calculate $n_2$ (or vice versa) for a desired power. 
	- The sample size ratio $n_2$ / $n_1$ to calculate both $n_1$ and $n_2$ for a desired power.
- **`lm_plot()`** Create a plot for a linear regression model that includes the line of best fit, confidence intervals, and prediction intervals. 

## Installation

You can install the development version of **biostats101** directly from GitHub using the `devtools` package. 
```r
# If you don't have `devtools` installed, you can install it with:
# install.packages("devtools")

# Then, install the package from GitHub:
devtools::install_github("marina-wiz/biostats101")
```

## Dependencies

This package has minimal dependencies:

- **`mean_CI()`**, **`power.paired.prop()`**, and **`power.2p.2n()`** do not require any additional R packages.
- **`lm_plot()`** requires the following R packages:
  - `dplyr`
  - `tidyr`
  - `ggplot2`

### Automatic Package Installation

By default, **`lm_plot()`** will check if these packages are installed and automatically install them if needed. You can also choose to skip the automatic installation by setting `install_packages = FALSE`.


## Usage

Here’s are examples of how to use the functions in **biostats101**:

### 1. `mean_CI`

```r
library(biostats101)

# Example data
values = c(5.2, 4.8, 6.3, 6.1, 7.2, 3.5, 4.9, 2.2, 3.7, 3.5, 8.9)

# Construct a 95% confidence interval for the mean
mean_CI(values, conf.level = 0.95, alternative = 'two.sided')
```

### 2. `power.paired.prop`
```r
library(biostats101)

# Calculate the power given the sample size for paired proportions
power.paired.prop(p1 = 0.1, p2 = 0.15, n = 900)

# Calculate the sample size given the power for paired proportions
power.paired.prop(p1 = 0.15, p2 = 0.1, n = 900, alternative='one.sided')
```

### 3. `power.2p.2n`
```r
library(biostats101)

# Calculate the power for independent proportions given the sample sizes
power.2p.2n(p1 = 0.45, p2 = 0.6, n1 = 260, n2 = 130)

# Calculate the sample size for independent proportions (default power = 0.8)
power.2p.2n(p1 = 0.45, p2 = 0.6)

# Calculate sample sizes for independent proportions given the nratio (n2/n1)
power.2p.2n(p1 = 0.44, p2 = 0.6, nratio = 2)

# Calculate the sample size n2 given sample size n1 for independent proportions 
power.2p.2n(p1 = 0.44, p2 = 0.6, n1 = 108)
```

### 4. `lm_plot`
```r
library(biostats101)

# Example dataset
mydata <- data.frame(
  x = rnorm(100, mean = 50, sd = 10),  
  y = 3 + 0.5 * rnorm(100, mean = 50, sd = 10) + rnorm(100) 
)

# Run a regression model
my_model <- lm(y ~ x, mydata)

# Create a plot with the line of best fit, confidence limits, and prediction limits
lm_plot(my_model) 

# Customize plot labels
lm_plot(my_model) + xlab("Your x-axis label") + ylab("Your y-axis label")
```

## References

The methods implemented in this package are based on the following works:
- Connor, R. J. (1987). Sample size for testing differences in proportions for the paired-sample design. *Biometrics*, 207-211. https://doi.org/10.2307/2531961.
- Fleiss, J. L., Levin, B., & Paik, M. C. (2013). *Statistical methods for rates and proportions*. John Wiley & Sons.
- Levin, B., & Chen, X. (1999). Is the one-half continuity correction used once or twice to derive a well-known approximate sample size formula to compare two independent binomial distributions?. *The American Statistician*, 53(1), 62-66. https://doi.org/10.1080/00031305.1999.10474431.
- McNemar, Q. (1947). Note on the sampling error of the difference between correlated proportions or percentages. *Psychometrika*, 12(2), 153-157. https://doi.org/10.1007/BF02295996.
