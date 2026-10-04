* **Script:** MAIC_glm.R, October 5, 2026
* **Version:** v1.0.0, October 5, 2026

# MAIC for Generalized Linear Models
 
This repository provides an R function for calculating the Modified Akaike

Information Criterion (MAIC) proposed in our study.
 
The function `MAIC_glm()` calculates the proposed MAIC for generalized linear

models fitted using the `glm()` function in R.
 
## Requirements
 
The following R package is required:
 
- `brglm2`
 
The package can be installed by
 
```r

install.packages("brglm2")

```
 
and loaded by
 
```r

library(brglm2)

```
 
## Files
 
- `MAIC_glm.R`: R function for calculating the proposed MAIC.

- `Example.R`: Example code illustrating the use of `MAIC_glm()`.

- `DataM_Logit_R0.csv`: Example data used in `Example.R`.
 
## Usage
 
First, load the required package and the function:
 
```r

library(brglm2)

source("MAIC_glm.R")

```
 
Suppose that `result` is a fitted model object returned by `glm()`, `Xs` is

the matrix of explanatory variables for the saturated model, and `ys` is the

response vector. The MAIC can then be calculated by
 
```r

MAIC_glm(result, Xs, ys)

```
 
The intercept should be included as a column of `Xs`.
 
See `Example.R` for a complete example.
 
## Supported Models
 
The current implementation supports the following distributions with their

canonical link functions:
 
- Gaussian distribution with the identity link

- Binomial distribution with the logit link

- Gamma distribution with the inverse link

- Poisson distribution with the log link

- Inverse Gaussian distribution with the inverse-squared link
 
An error is returned if a noncanonical link function or an unsupported

distribution is specified.
 
## Note on the Dispersion Parameter
 
The theoretical results underlying the proposed MAIC apply to generalized

linear models with a known dispersion parameter. For models with an unknown

dispersion parameter, `MAIC_glm()` calculates the MAIC by treating the

estimated dispersion parameter as if it were known.
 
## Reference
 
The manuscript describing the proposed MAIC is currently under preparation.
