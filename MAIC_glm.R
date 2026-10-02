##############################################################################
# MAIC_glm
# Script for calculating the modified AIC for generalized linear models
# with canonical links
#
# [Required package]
#	brglm2		Estimation by maximizing the penalized log-likelihood based on
#				Firth's method
#
# [Inputs]
#	object			obj		A fitted model object returned by the glm function
#	matrix[n, m]	Xs		Matrix of explanatory variables of the saturated model
#	numeric[n]		ys		Vector of response variable
#
# [Output]
#	numeric[1]		maic	MAIC
##
# [Details]
#	Calculates the Modified AIC for the generalized linear model fitted in obj.
#	The saturated model is fitted using Firth's method, and the resulting 
#	estimates are used to correct the AIC.
#	The intercept should be included as a column of Xs.
##############################################################################

MAIC_glm <- function (obj, Xs, ys) {

	#-------------------------------------------------------
	# Confirming Input
	#-------------------------------------------------------

	if (! "glm" %in% class(obj)) {

		# Return an error if the class is not 'glm'
		stop("The class of 'obj' must be 'glm'.")

	}


	#----------------------------------------------------------------------
	# Verifying the distribution and link function, and setting parameters
	#----------------------------------------------------------------------

	# Extract the required information from the object returned by glm
	fam <- obj$family
	nu <- NULL # dispersion parameter

	if (fam$family == "gaussian") {

		if (fam$link != "identity") {

			# Return an error if the link function is not 'identity' for 
			# the normal distribution
			stop("The distribution is 'gaussian', but the link function is not 'identity'.")
		}

		# dispersion parameter
		nu <- sigma(obj)^2

		phi <- nu
		a <- phi
		bf <- function (theta) { theta^2 / 2 }
		cf <- function (y, phi) { -(y^2 / phi + log(2 * pi * phi)) / 2 }

		mu <- function (theta) { theta }
		kappa2 <- function (theta, phi) { rep(phi, length(theta)) }
		kappa3 <- function (theta, phi) { rep(0, length(theta)) }
		kappa4 <- function (theta, phi) { rep(0, length(theta)) }

	} else if (fam$family == "binomial") {

		if (fam$link != "logit") {

			# Return an error if the link function is not 'logit' for 
			# the binominal distribution

			stop("The distribution is 'binomial', but the link function is not 'logit'.")

		}

		phi <- 1
		a <- phi
		bf <- function (theta) { log(1 + exp(theta)) }
		cf <- function (y, phi) { rep(0, length(y)) }

		mu <- function (theta) { 1 / (1 + exp(-theta)) }
		kappa2 <- function (theta, phi) {
			mu(theta) * (1 - mu(theta))
		}
		kappa3 <- function (theta, phi) {
			kappa2(theta, phi) * (1 - 2 * mu(theta))
		}
		kappa4 <- function (theta, phi) {
			k2 <- kappa2(theta, phi)
			k2 * (1 - 6 * k2)
		}

	} else if (fam$family == "Gamma") {

		if (fam$link != "inverse") {

			# Return an error if the link function is not 'inverse' for 
			# the Gamma distribution
			stop("The distribution is 'Gamma', but the link function is not 'inverse'.")

		}
		if (! is.numeric(nu)) {

			# Return an error if nu is not numeric
			stop("'nu' must be numeric.")

		}
		if (! (nu > 0)) {

			# Return an error if nu is not positive
			stop("'nu' must be positive.")

		}

		# dispersion parameter
		# Reciprocal for Gamma regression
		nu <- 1 / (summary(obj)$dispersion)

		phi <- nu
		a <- phi^(-1)
		bf <- function (theta) { -log(-theta) }
		cf <- function (y, phi) { phi * log(phi * y) - log(y) - log(gamma(phi)) }

		mu <- function (theta) { -1 / theta }
		kappa2 <- function (theta, phi) { 1 / (phi * theta^2) }
		kappa3 <- function (theta, phi) { -2 / (phi^2 * theta^3) }
		kappa4 <- function (theta, phi) { 6 / (phi^3 * theta^4) }

	} else if (fam$family == "poisson") {

		if (fam$link != "log") {

			# Return an error if the link function is not 'log' for 
			# the Poisson distribution
			stop("The distribution is 'poisson', but the link function is not 'log'.")

		}

		phi <- 1
		a <- phi
		bf <- function (theta) { exp(theta) }
		cf <- function (y, phi) { log(factorial(y)) } 

		nu <- NULL

		mu <- function (theta) { exp(theta) }
		kappa2 <- function (theta, phi) { mu(theta) }
		kappa3 <- function (theta, phi) { mu(theta) }
		kappa4 <- function (theta, phi) { mu(theta) }

	} else if (fam$family == "inverse.gaussian") {

		if (fam$link != "1/mu^2") {

			# Return an error if the link function is not '1/mu^2' for 
			# the inverse Gaussian distribution
			stop("The distribution is 'inverse.gaussian', but the link function is not '1/mu^2'.")

		}

		# dispersion parameter
		nu <- sigma(obj)^2

		phi <- nu
		a <- phi
		bf <- function (theta) { sqrt(-2 * theta) }
		cf <- function (y, phi) { -(log(2 * pi * phi * y^3) + (1 / (phi * y))) / 2 }

		mu <- function (theta) { 1 / sqrt(-2 * theta) }
		kappa2 <- function (theta, phi) { a * (-2 * theta)^(-3 / 2) }
		kappa3 <- function (theta, phi) { 3 * a^2 * (-2 * theta)^(-5 / 2) }
		kappa4 <- function (theta, phi) { 15 * a^3 * (-2 * theta)^(-7 / 2) }

	} else {

		# for distributions other than gaussian, binomial, Gamma, poisson, and inverse.gaussiann
		stop("The distribution must be 'gaussian', 'binomial', 'Gamma', 'poisson', or 'inverse.gaussian'.")

	}


	#-------------------------------------------------------
	# Calculating the values of functions in MAIC
	#-------------------------------------------------------

	# delta
	delta <- function (k3, W, a) {

		n <- length(k3)

		r <- rep(0, n)
		for (i in seq_len(n)) {

			w_i <- W[i, , drop = T]
			w_ii <- w_i[i]

			r[i] <- sum(k3 * (w_i^3 / w_ii))

		}

		r / (n * a^3)

	}


	#-------------------------------------------------------
	# Preliminary
	#-------------------------------------------------------

	n <- nrow(Xs)
	k_s <- ncol(Xs)

	d <- as.data.frame(cbind(Y = ys, Xs))

	#---------------------------------------------
	# Candidate model
	#---------------------------------------------

	X_c <- model.matrix(obj)

	eta_c <- predict(obj, newdata = d)	# linear predictor
	k2_c <- kappa2(eta_c, phi)
	k3_c <- kappa3(eta_c, phi)

	G2_c <- (t(X_c) %*% (k2_c * X_c)) / (n * a^2)
	G2_c_inv <- solve(G2_c)

	W_c <- X_c %*% G2_c_inv %*% t(X_c)

	del <- delta(k3_c, W_c, a)


	#---------------------------------------------
	# Saturated model
	#---------------------------------------------

	# Estimation of the saturated model using Firth's Method
	obj_F <- glm(Y ~ . + 0, data = d, family = fam, method = "brglm_fit",
				control = brglmControl(type = "MPL_Jeffreys"))

	X_s <- Xs

	theta_star <- predict(obj_F, newdata = d)	# linear predictor
	k2_s <- kappa2(theta_star, phi)
	k4_s <- kappa4(theta_star, phi)

	G2_s <- (t(X_s) %*% (k2_s * X_s)) / (n * a^2)
	G2_s_inv <- solve(G2_s)

	W_s <- X_s %*% G2_s_inv %*% t(X_s)

	h_num <- 2 * n * a^2 * k2_s
	h_den <- diag(W_s) * k4_s + h_num
	h_s <- h_num / h_den


	#---------------------------------------------
	# G2_star
	#---------------------------------------------

	v <- theta_star + (del / (2 * n))
	k2_v <- kappa2(v, phi)

	k2_star <- h_s * k2_v
	G2_star <- (t(X_c) %*% (k2_star * X_c)) / (n * a)


	#-------------------------------------------------------
	# Calculation for MAIC
	#-------------------------------------------------------

	# Log-likelihood
	ll <- logLik(obj)

	# Estimated generalized degrees of freedom
	G <- G2_star %*% G2_c_inv
	pen <- sum(diag(G))

	# MAIC
	maic <- -2 * as.numeric(ll) + 2 * pen

	# Return
	return(maic)

}

# EoF
