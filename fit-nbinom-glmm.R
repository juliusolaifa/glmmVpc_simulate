fit_nbinom_glmm <- function(
    data,
    fit_null = TRUE
) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop("Package 'glmmTMB' is required.")
    }

    required <- c("id", "x", "y")

    if (!all(required %in% names(data))) {
        stop("data must contain id, x, and y.")
    }

    data$id <- factor(data$id)

    # Alternative: random intercept and slope
    fit_alt <- glmmTMB::glmmTMB(
        y ~ x + (1 + x | id),
        family = glmmTMB::nbinom2(link = "log"),
        data = data
    )

    # Extract alternative-model parameters
    beta_hat <- unname(
        glmmTMB::fixef(fit_alt)$cond
    )

    Gamma_hat <- as.matrix(
        glmmTMB::VarCorr(fit_alt)$cond$id
    )

    theta_hat <- unname(
        stats::sigma(fit_alt)
    )

    ll_alt <- as.numeric(
        stats::logLik(fit_alt)
    )

    # Convergence diagnostics
    alt_ok <- (
        fit_alt$fit$convergence == 0 &&
        isTRUE(fit_alt$sdr$pdHess) &&
        is.finite(ll_alt)
    )

    result <- list(
        alternative = list(
            fit = fit_alt,
            beta = beta_hat,
            Gamma = Gamma_hat,
            theta = theta_hat,
            logLik = ll_alt,
            converged = alt_ok
        ),
        null = NULL
    )

    if (fit_null) {

        # Null: Gamma = 0
        fit0 <- glmmTMB::glmmTMB(
            y ~ x,
            family = glmmTMB::nbinom2(link = "log"),
            data = data
        )

        ll_null <- as.numeric(
            stats::logLik(fit0)
        )

        null_ok <- (
            fit0$fit$convergence == 0 &&
            isTRUE(fit0$sdr$pdHess) &&
            is.finite(ll_null)
        )

        result$null <- list(
            fit = fit0,
            beta = unname(
                glmmTMB::fixef(fit0)$cond
            ),
            Gamma = matrix(0, 2, 2),
            theta = unname(
                stats::sigma(fit0)
            ),
            logLik = ll_null,
            converged = null_ok
        )
    }

    result
}
