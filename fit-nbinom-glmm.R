fit_nbinom_glmm <- function(
    data,
    fit_null = TRUE,
    grad_tol = 1e-4,
    singular_tol = 1e-6
) {
    # --------------------------------------------------
    # 0. Validate inputs
    # --------------------------------------------------

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop("Package 'glmmTMB' is required.")
    }

    required <- c("id", "x", "y")

    if (!is.data.frame(data) ||
        !all(required %in% names(data))) {
        stop("data must be a data frame containing id, x, and y.")
    }

    if (!is.logical(fit_null) ||
        length(fit_null) != 1L ||
        is.na(fit_null)) {
        stop("fit_null must be TRUE or FALSE.")
    }

    if (!is.numeric(grad_tol) ||
        length(grad_tol) != 1L ||
        !is.finite(grad_tol) ||
        grad_tol <= 0) {
        stop("grad_tol must be positive and finite.")
    }

    if (!is.numeric(singular_tol) ||
        length(singular_tol) != 1L ||
        !is.finite(singular_tol) ||
        singular_tol < 0) {
        stop("singular_tol must be nonnegative and finite.")
    }

    if (anyNA(data[, required]) ||
        !is.numeric(data$x) ||
        !is.numeric(data$y) ||
        any(!is.finite(data$x)) ||
        any(!is.finite(data$y)) ||
        any(data$y < 0) ||
        any(data$y != floor(data$y))) {
        stop("x must be finite and y must contain nonnegative integer counts.")
    }

    data$id <- factor(data$id)

    if (nlevels(data$id) < 2L) {
        stop("At least two clusters are required.")
    }

    # --------------------------------------------------
    # 1. Helper: fitting diagnostics
    # --------------------------------------------------

    fit_diagnostics <- function(fit) {

        ll <- tryCatch(
            as.numeric(stats::logLik(fit)),
            error = function(e) NA_real_
        )

        likelihood_ok <- (
            length(ll) == 1L &&
            is.finite(ll)
        )

        convergence_code <- fit$fit$convergence

        optimizer_ok <- isTRUE(
            convergence_code == 0L
        )

        hessian_ok <- isTRUE(
            fit$sdr$pdHess
        )

        gradient <- tryCatch(
            fit$obj$gr(fit$fit$par),
            error = function(e) NA_real_
        )

        gradient_ok <- (
            length(gradient) > 0L &&
            all(is.finite(gradient))
        )

        max_abs_gradient <- if (gradient_ok) {
            max(abs(gradient))
        } else {
            NA_real_
        }

        gradient_norm <- if (gradient_ok) {
            sqrt(sum(gradient^2))
        } else {
            NA_real_
        }

        stationary <- isTRUE(
            gradient_ok &&
            max_abs_gradient <= grad_tol
        )

        # Optimizer-reported convergence
        converged <- optimizer_ok && likelihood_ok

        # Separate numerical screening criterion.
        # This is not proof of a global maximum.
        numerically_acceptable <- (
            likelihood_ok &&
            hessian_ok &&
            stationary
        )

        list(
            convergence_code = convergence_code,
            optimizer_message = fit$fit$message,
            optimizer_ok = optimizer_ok,
            hessian_ok = hessian_ok,
            likelihood_ok = likelihood_ok,
            max_abs_gradient = max_abs_gradient,
            gradient_norm = gradient_norm,
            stationary = stationary,
            converged = converged,
            numerically_acceptable = numerically_acceptable
        )
    }

    # --------------------------------------------------
    # 2. Fit alternative model: Gamma >= 0
    # --------------------------------------------------

    fit_alt <- glmmTMB::glmmTMB(
        y ~ x + (1 + x | id),
        family = glmmTMB::nbinom2(link = "log"),
        data = data
    )

    alt_diag <- fit_diagnostics(fit_alt)

    beta_alt <- unname(
        glmmTMB::fixef(fit_alt)$cond
    )

    Gamma_alt <- as.matrix(
        glmmTMB::VarCorr(fit_alt)$cond$id
    )

    # Symmetrize and remove glmmTMB-specific attributes
    Gamma_alt <- (
        Gamma_alt + t(Gamma_alt)
    ) / 2

    Gamma_alt <- matrix(
        as.numeric(Gamma_alt),
        nrow = 2L,
        ncol = 2L,
        dimnames = list(
            c("(Intercept)", "x"),
            c("(Intercept)", "x")
        )
    )

    theta_alt <- unname(
        stats::sigma(fit_alt)
    )

    # --------------------------------------------------
    # 3. Covariance singularity diagnostics
    # --------------------------------------------------

    if (all(is.finite(Gamma_alt))) {

        eigenvalues <- eigen(
            Gamma_alt,
            symmetric = TRUE,
            only.values = TRUE
        )$values

        min_eigenvalue <- min(eigenvalues)
        max_eigenvalue <- max(eigenvalues)

        near_singular <- (
            min_eigenvalue <=
                singular_tol * max(1, max_eigenvalue)
        )

        eigenvalue_ratio <- if (max_eigenvalue > 0) {
            min_eigenvalue / max_eigenvalue
        } else {
            NA_real_
        }

    } else {

        min_eigenvalue <- NA_real_
        max_eigenvalue <- NA_real_
        near_singular <- NA
        eigenvalue_ratio <- NA_real_
    }

    alt_diag$min_eigenvalue <- min_eigenvalue
    alt_diag$max_eigenvalue <- max_eigenvalue
    alt_diag$eigenvalue_ratio <- eigenvalue_ratio
    alt_diag$near_singular <- near_singular

    alternative <- list(
        fit = fit_alt,
        beta = beta_alt,
        Gamma = Gamma_alt,
        theta = theta_alt,
        logLik = as.numeric(stats::logLik(fit_alt)),
        converged = alt_diag$converged,
        diagnostics = alt_diag
    )

    # --------------------------------------------------
    # 4. Fit null model: Gamma = 0
    # --------------------------------------------------

    null <- NULL

    if (fit_null) {

        fit0 <- glmmTMB::glmmTMB(
            y ~ x,
            family = glmmTMB::nbinom2(link = "log"),
            data = data
        )

        null_diag <- fit_diagnostics(fit0)

        beta_null <- unname(
            glmmTMB::fixef(fit0)$cond
        )

        theta_null <- unname(
            stats::sigma(fit0)
        )

        null <- list(
            fit = fit0,
            beta = beta_null,
            Gamma = matrix(
                0,
                nrow = 2L,
                ncol = 2L
            ),
            theta = theta_null,
            logLik = as.numeric(stats::logLik(fit0)),
            converged = null_diag$converged,
            diagnostics = null_diag
        )
    }

    # --------------------------------------------------
    # 5. Return
    # --------------------------------------------------

    list(
        alternative = alternative,
        null = null
    )
}
