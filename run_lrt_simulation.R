run_lrt_simulation <- function(
    scenario,
    R = 1000L,
    alpha = 0.05,
    weights = c("oracle", "plugin"),
    seed = 123,
    rel.tol = 1e-7,
    subdivisions = 200L,
    singular_tol = 1e-6,
    progress = TRUE
) {

    weights <- match.arg(
        weights,
        choices = c("oracle", "plugin"),
        several.ok = TRUE
    )

    weights <- unique(weights)

    stopifnot(
        length(R) == 1L,
        is.finite(R),
        R >= 1,
        R == as.integer(R),
        length(alpha) == 1L,
        is.finite(alpha),
        alpha > 0,
        alpha < 1,
        length(seed) == 1L,
        is.finite(seed),
        length(singular_tol) == 1L,
        is.finite(singular_tol),
        singular_tol >= 0
    )

    required <- c(
        "n_clusters",
        "cluster_size",
        "beta",
        "Gamma",
        "theta"
    )

    if (!all(required %in% names(scenario))) {
        stop("Missing required scenario parameters.")
    }

    # Cluster-specific design matrices
    make_design <- function(dat) {

        X_all <- model.matrix(~ x, data = dat)

        idx <- split(
            seq_len(nrow(dat)),
            dat$id
        )

        X <- lapply(
            idx,
            function(i) X_all[i, , drop = FALSE]
        )

        list(X = X, Z = X)
    }

    # Preserve the caller's random-number state
    had_seed <- exists(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
    )

    if (had_seed) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
    }

    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
        )) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    set.seed(seed)

    # -------------------------------------------------
    # 1. Oracle weights: compute once per scenario
    # -------------------------------------------------

    oracle_weights <- NULL

    if ("oracle" %in% weights) {

        # Use the scenario's fixed design.
        # This does not consume simulation RNG draws.
        m <- scenario$cluster_size
        n <- scenario$n_clusters

        x_design <- rep(
            seq(0, 1, length.out = m),
            times = n
        )

        dat_design <- data.frame(
            id = factor(rep(seq_len(n), each = m)),
            x = x_design
        )

        design <- make_design(dat_design)

        fim_true <- glmm_fim_theoretical(
            X = design$X,
            Z = design$Z,
            beta = scenario$beta,
            family = "nbinom",
            dispersion = scenario$theta
        )

        oracle_weights <- chibar_weights_theoretical(
            information = fim_true$I_metric,
            rel.tol = rel.tol,
            subdivisions = subdivisions
        )$w
    }

    # -------------------------------------------------
    # 2. Run one replicate
    # -------------------------------------------------

    one_replicate <- function(r) {

        out <- list(
            replicate = r,
            converged = FALSE,
            singular = NA,
            statistic = NA_real_,
            p_oracle = NA_real_,
            p_plugin = NA_real_,
            error = NA_character_
        )

        tryCatch({

            # Generate data
            dat <- simulate_nbinom_glmm(
                n_clusters = scenario$n_clusters,
                cluster_size = scenario$cluster_size,
                beta = scenario$beta,
                Gamma = scenario$Gamma,
                theta = scenario$theta
            )

            # Fit null and alternative
            fit <- fit_nbinom_glmm(dat)

            if (!isTRUE(fit$alternative$converged) ||
                !isTRUE(fit$null$converged)) {
                stop("Model convergence check failed.")
            }

            out$converged <- TRUE

            # Flag covariance matrices close to singularity
            G <- fit$alternative$Gamma

            ev <- eigen(
                (G + t(G)) / 2,
                symmetric = TRUE,
                only.values = TRUE
            )$values

            out$singular <- (
                min(ev) <= singular_tol * max(1, max(ev))
            )

            # LRT statistic
            T <- 2 * (
                fit$alternative$logLik -
                fit$null$logLik
            )

            if (!is.finite(T) || T < -1e-8) {
                stop("Invalid likelihood ratio statistic.")
            }

            T <- max(0, T)
            out$statistic <- T

            # Oracle calibration
            if ("oracle" %in% weights) {

                test <- glmm_chibar_lrt(
                    logLik_null = fit$null$logLik,
                    logLik_alt = fit$alternative$logLik,
                    weights = oracle_weights
                )

                out$p_oracle <- test$p_value
            }

            # Plug-in calibration
            if ("plugin" %in% weights) {

                design <- make_design(dat)

                fim <- glmm_fim_theoretical(
                    X = design$X,
                    Z = design$Z,
                    beta = fit$null$beta,
                    family = "nbinom",
                    dispersion = fit$null$theta
                )

                w <- chibar_weights_theoretical(
                    information = fim$I_metric,
                    rel.tol = rel.tol,
                    subdivisions = subdivisions
                )$w

                test <- glmm_chibar_lrt(
                    logLik_null = fit$null$logLik,
                    logLik_alt = fit$alternative$logLik,
                    weights = w
                )

                out$p_plugin <- test$p_value
            }

            out

        }, error = function(e) {

            out$error <- conditionMessage(e)
            out
        })
    }

    # -------------------------------------------------
    # 3. Repeat simulations
    # -------------------------------------------------

    results <- vector("list", R)

    for (r in seq_len(R)) {

        results[[r]] <- one_replicate(r)

        if (progress && (r %% 100L == 0L || r == R)) {
            message("Completed ", r, " / ", R)
        }
    }

    results <- do.call(
        rbind,
        lapply(results, as.data.frame)
    )

    rownames(results) <- NULL

    # -------------------------------------------------
    # 4. Summarize results
    # -------------------------------------------------

    summarize_method <- function(p) {

        valid <- is.finite(p)
        n_valid <- sum(valid)

        if (n_valid == 0L) {
            return(data.frame(
                n_valid = 0L,
                rejection_rate = NA_real_,
                mcse = NA_real_
            ))
        }

        rate <- mean(p[valid] < alpha)

        data.frame(
            n_valid = n_valid,
            rejection_rate = rate,
            mcse = sqrt(rate * (1 - rate) / n_valid)
        )
    }

    summary <- do.call(
        rbind,
        lapply(weights, function(method) {

            p <- results[[paste0("p_", method)]]

            cbind(
                method = method,
                summarize_method(p)
            )
        })
    )

    rownames(summary) <- NULL

    # Diagnostics
    diagnostics <- list(
        n_simulations = R,
        n_converged = sum(results$converged),
        convergence_rate = mean(results$converged),
        n_singular = sum(results$singular, na.rm = TRUE),
        n_errors = sum(!is.na(results$error))
    )

    list(
        scenario = scenario,
        alpha = alpha,
        summary = summary,
        diagnostics = diagnostics,
        results = results,
        oracle_weights = oracle_weights
    )
}
