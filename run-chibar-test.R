run_chibar_test <- function(
    fit,
    X,
    Z,
    family
) {
    if (!fit$alternative$converged ||
        !fit$null$converged) {
        stop("Both models must converge.")
    }

    # Fisher information under H0: Gamma = 0
    fim <- glmm_fim_theoretical(
        X = X,
        Z = Z,
        beta = fit$null$beta,
        family = family,
        dispersion = fit$null$theta
    )

    # Nuisance-adjusted information metric
    I_metric <- fim$I_metric

    # Exact chi-bar-square weights
    weights <- chibar_weights_theoretical(
        information = I_metric
    )

    # LRT and p-value
    lrt <- glmm_chibar_lrt(
        logLik_null = fit$null$logLik,
        logLik_alt = fit$alternative$logLik,
        weights = weights
    )

    lrt
}
