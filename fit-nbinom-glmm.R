
fit_nbinom_glmm <- function(
    data,
    fit_null = TRUE,
    grad_tol = 1e-4,
    singular_tol = 1e-6,
    refit = TRUE
) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop("Package 'glmmTMB' is required.")
    }

    required <- c("id", "x", "y")

    if (!is.data.frame(data) ||
        !all(required %in% names(data))) {
        stop("data must contain id, x, and y.")
    }

    if (!is.logical(fit_null) || length(fit_null) != 1L ||
        is.na(fit_null)) {
        stop("fit_null must be TRUE or FALSE.")
    }

    if (!is.logical(refit) || length(refit) != 1L ||
        is.na(refit)) {
        stop("refit must be TRUE or FALSE.")
    }

    if (!is.numeric(grad_tol) || length(grad_tol) != 1L ||
        !is.finite(grad_tol) || grad_tol <= 0) {
        stop("grad_tol must be positive and finite.")
    }

    if (!is.numeric(singular_tol) ||
        length(singular_tol) != 1L ||
        !is.finite(singular_tol) || singular_tol < 0) {
        stop("singular_tol must be nonnegative and finite.")
    }

    if (anyNA(data[, required]) ||
        !is.numeric(data$x) ||
        !is.numeric(data$y) ||
        any(!is.finite(data$x)) ||
        any(!is.finite(data$y)) ||
        any(data$y < 0) ||
        any(data$y != floor(data$y))) {
        stop("Invalid x or y values.")
    }

    data$id <- factor(data$id)

    if (nlevels(data$id) < 2L) {
        stop("At least two clusters are required.")
    }

    # --------------------------------------------------
    # 1. Extract numerical diagnostics
    # --------------------------------------------------

    fit_diagnostics <- function(fit) {

        ll <- tryCatch(
            as.numeric(stats::logLik(fit)),
            error = function(e) NA_real_
        )

        objective <- tryCatch(
            as.numeric(fit$fit$objective),
            error = function(e) NA_real_
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

        convergence_code <- fit$fit$convergence

        optimizer_ok <- isTRUE(convergence_code == 0L)
        hessian_ok <- isTRUE(fit$sdr$pdHess)
        likelihood_ok <- length(ll) == 1L && is.finite(ll)
        objective_ok <- length(objective) == 1L &&
            is.finite(objective)

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
            objective_ok = objective_ok,
            objective = objective,
            objective_logLik = if (objective_ok) {
                -objective
            } else {
                NA_real_
            },
            max_abs_gradient = max_abs_gradient,
            gradient_norm = gradient_norm,
            stationary = stationary,
            converged = optimizer_ok && likelihood_ok,
            numerically_acceptable = (
                optimizer_ok &&
                likelihood_ok &&
                hessian_ok &&
                stationary
            )
        )
    }

    # --------------------------------------------------
    # 2. Fit alternative: default optimizer
    # --------------------------------------------------

    fit_alt_default <- glmmTMB::glmmTMB(
        y ~ x + (1 + x | id),
        family = glmmTMB::nbinom2(link = "log"),
        data = data
    )

    diag_default <- fit_diagnostics(fit_alt_default)

    fit_alt <- fit_alt_default
    alt_diag <- diag_default
    selected_optimizer <- "nlminb"

    refit_diagnostics <- list(
        nlminb = diag_default
    )

    # --------------------------------------------------
    # 3. Conditional refit using BFGS
    # --------------------------------------------------

    needs_refit <- (
        !diag_default$optimizer_ok ||
        !diag_default$hessian_ok ||
        !diag_default$likelihood_ok ||
        !diag_default$stationary
    )

    if (refit && needs_refit) {

        fit_alt_bfgs <- tryCatch(
            glmmTMB::glmmTMB(
                y ~ x + (1 + x | id),
                family = glmmTMB::nbinom2(link = "log"),
                data = data,
                control = glmmTMB::glmmTMBControl(
                    optimizer = stats::optim,
                    optArgs = list(method = "BFGS")
                )
            ),
            error = function(e) NULL
        )

        if (!is.null(fit_alt_bfgs)) {

            diag_bfgs <- fit_diagnostics(fit_alt_bfgs)

            refit_diagnostics$bfgs <- diag_bfgs

            # Prefer numerically acceptable fits.
            # Among equally acceptable candidates,
            # prefer the smaller objective.

            choose_bfgs <- FALSE

            if (diag_bfgs$numerically_acceptable &&
                !diag_default$numerically_acceptable) {

                choose_bfgs <- TRUE

            } else if (
                diag_bfgs$numerically_acceptable ==
                    diag_default$numerically_acceptable &&
                diag_bfgs$objective_ok &&
                (
                    !diag_default$objective_ok ||
                    diag_bfgs$objective <
                        diag_default$objective
                )
            ) {

                choose_bfgs <- TRUE
            }

            if (choose_bfgs) {
                fit_alt <- fit_alt_bfgs
                alt_diag <- diag_bfgs
                selected_optimizer <- "BFGS"
            }
        }
    }

    alt_diag$selected_optimizer <- selected_optimizer

    # --------------------------------------------------
    # 4. Extract alternative estimates
    # --------------------------------------------------

    beta_alt <- unname(
        glmmTMB::fixef(fit_alt)$cond
    )

    Gamma_alt <- as.matrix(
        glmmTMB::VarCorr(fit_alt)$cond$id
    )

    Gamma_alt <- (Gamma_alt + t(Gamma_alt)) / 2

    Gamma_alt <- matrix(
        as.numeric(Gamma_alt),
        nrow = 2L,
        ncol = 2L,
        dimnames = list(
            c("(Intercept)", "x"),
            c("(Intercept)", "x")
        )
    )

    theta_alt <- unname(stats::sigma(fit_alt))

    # --------------------------------------------------
    # 5. Covariance singularity diagnostics
    # --------------------------------------------------

    if (all(is.finite(Gamma_alt))) {

        eigenvalues <- eigen(
            Gamma_alt,
            symmetric = TRUE,
            only.values = TRUE
        )$values

        min_eigenvalue <- min(eigenvalues)
        max_eigenvalue <- max(eigenvalues)

        eigenvalue_ratio <- if (max_eigenvalue > 0) {
            min_eigenvalue / max_eigenvalue
        } else {
            NA_real_
        }

        near_singular <- (
            min_eigenvalue <=
                singular_tol * max(1, max_eigenvalue)
        )

    } else {

        min_eigenvalue <- NA_real_
        max_eigenvalue <- NA_real_
        eigenvalue_ratio <- NA_real_
        near_singular <- NA
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
        logLik = tryCatch(
            as.numeric(stats::logLik(fit_alt)),
            error = function(e) NA_real_
        ),
        objective = alt_diag$objective,
        objective_logLik = alt_diag$objective_logLik,
        converged = alt_diag$converged,
        diagnostics = alt_diag,
        refit_diagnostics = refit_diagnostics
    )

    # --------------------------------------------------
    # 6. Fit null model
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
            beta = unname(glmmTMB::fixef(fit0)$cond),
            Gamma = matrix(0, 2L, 2L),
            theta = unname(stats::sigma(fit0)),
            logLik = tryCatch(
                as.numeric(stats::logLik(fit0)),
                error = function(e) NA_real_
            ),
            objective = null_diag$objective,
            objective_logLik = null_diag$objective_logLik,
            converged = null_diag$converged,
            diagnostics = null_diag
        )
    }

    # --------------------------------------------------
    # 7. Return
    # --------------------------------------------------

    list(
        alternative = alternative,
        null = null
    )
}
