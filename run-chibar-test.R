run_chibar_test <- function(
    fit,
    X,
    Z,
    family = "nbinom",
    rel.tol = 1e-7,
    subdivisions = 200L
) {
    # Check convergence
    if (!isTRUE(fit$alternative$converged) ||
        !isTRUE(fit$null$converged)) {
        stop("Both models must converge.")
    }

    # 1. Theoretical Fisher information under H0
    fim <- glmm_fim_theoretical(
        X = X,
        Z = Z,
        beta = fit$null$beta,
        family = family,
        dispersion = fit$null$theta
    )

    # 2. Exact chi-bar-square weights
    weights_result <- chibar_weights_theoretical(
        information = fim$I_metric,
        rel.tol = rel.tol,
        subdivisions = subdivisions
    )

    # 3. Likelihood ratio test
    lrt <- glmm_chibar_lrt(
        logLik_null = fit$null$logLik,
        logLik_alt = fit$alternative$logLik,
        weights = weights_result$w
    )

    list(
        statistic = lrt$statistic,
        p_value = lrt$p_value,
        weights = weights_result$w,
        integration_error = weights_result$abs.error,
        I_metric = fim$I_metric
    )
}
