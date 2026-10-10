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

    # --------------------------------------------------
    # Helper: extract fitting diagnostics
    # --------------------------------------------------

    fit_diagnostics <- function(fit, grad_tol = 1e-4) {

        ll <- as.numeric(stats::logLik(fit))
    
        convergence_code <- fit$fit$convergence
    
        optimizer_ok <- isTRUE(convergence_code == 0L)
    
        hessian_ok <- isTRUE(fit$sdr$pdHess)
    
        likelihood_ok <- (
            length(ll) == 1L &&
            is.finite(ll)
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
    
        stationary <- isTRUE(
            gradient_ok &&
            max_abs_gradient <= grad_tol
        )
    
        list(
            convergence_code = convergence_code,
            optimizer_message = fit$fit$message,
            optimizer_ok = optimizer_ok,
            hessian_ok = hessian_ok,
            likelihood_ok = likelihood_ok,
            max_abs_gradient = max_abs_gradient,
            stationary = stationary,
            converged = optimizer_ok && likelihood_ok,
            numerically_acceptable = likelihood_ok &&
                hessian_ok &&
                stationary
        )
    }
    # --------------------------------------------------
    # 1. Alternative model
    # --------------------------------------------------

    fit_alt <- glmmTMB::glmmTMB(
        y ~ x + (1 + x | id),
        family = glmmTMB::nbinom2(link = "log"),
        data = data
    )

    alt_diag <- fit_diagnostics(fit_alt)

    Gamma_hat <- as.matrix(
        glmmTMB::VarCorr(fit_alt)$cond$id
    )

    Gamma_hat <- (
        Gamma_hat + t(Gamma_hat)
    ) / 2

    alt <- list(
        fit = fit_alt,

        beta = unname(
            glmmTMB::fixef(fit_alt)$cond
        ),

        Gamma = Gamma_hat,

        theta = unname(
            stats::sigma(fit_alt)
        ),

        logLik = as.numeric(
            stats::logLik(fit_alt)
        ),

        converged = alt_diag$converged,

        diagnostics = alt_diag
    )

    # --------------------------------------------------
    # 2. Null model
    # --------------------------------------------------

    null <- NULL

    if (fit_null) {

        fit0 <- glmmTMB::glmmTMB(
            y ~ x,
            family = glmmTMB::nbinom2(link = "log"),
            data = data
        )

        null_diag <- fit_diagnostics(fit0)

        null <- list(
            fit = fit0,

            beta = unname(
                glmmTMB::fixef(fit0)$cond
            ),

            Gamma = matrix(0, 2, 2),

            theta = unname(
                stats::sigma(fit0)
            ),

            logLik = as.numeric(
                stats::logLik(fit0)
            ),

            converged = null_diag$converged,

            diagnostics = null_diag
        )
    }

    # --------------------------------------------------
    # 3. Return results
    # --------------------------------------------------

    list(
        alternative = alt,
        null = null
    )
}
