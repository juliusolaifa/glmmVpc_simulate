run_one_replicate <- function(
    scenario,
    family
) {
    # Generate one dataset
    dat <- simulate_nbinom_glmm(
        n_clusters = scenario$n_clusters,
        cluster_size = scenario$cluster_size,
        beta = scenario$beta,
        Gamma = scenario$Gamma,
        theta = scenario$theta
    )

    # Fit alternative and null models
    fit <- fit_nbinom_glmm(dat)

    # Check convergence
    if (!fit$alternative$converged ||
        !fit$null$converged) {
        return(list(
            converged = FALSE,
            lrt = NA_real_,
            p_value = NA_real_
        ))
    }

    # Construct cluster-level design matrices
    X <- split(
        model.matrix(~ x, data = dat),
        dat$id
    )

    Z <- X

    # Calculate LRT
    test <- run_chibar_test(
        fit = fit,
        X = X,
        Z = Z,
        family = family
    )

    list(
        converged = TRUE,
        beta_hat = fit$alternative$beta,
        Gamma_hat = fit$alternative$Gamma,
        theta_hat = fit$alternative$theta,
        lrt = test$statistic,
        p_value = test$p_value
    )
}
