run_chibar_test <- function(
    fit,
    X,
    Z,
    weights_fun,
    family = "nbinom"
) {
    if (!isTRUE(fit$alternative$converged) ||
        !isTRUE(fit$null$converged)) {
        stop("Both models must converge.")
    }

    # Boundary Fisher information at null estimates
    fim <- glmm_fim_theoretical(
        X = X,
        Z = Z,
        beta = fit$null$beta,
        family = family,
        dispersion = fit$null$theta
    )

    # Exact weights from nuisance-adjusted information
    weights <- weights_fun(fim$I_metric)

    # LRT statistic and p-value
    test <- glmm_chibar_lrt(
        logLik_null = fit$null$logLik,
        logLik_alt = fit$alternative$logLik,
        weights = weights
    )

    list(
        statistic = test$statistic,
        p_value = test$p_value,
        weights = weights,
        I_metric = fim$I_metric
    )
}
