run_lrt_simulation <- function(
    scenario,
    R = 1000L,
    alpha = 0.05,
    weights = c("oracle", "plugin"),
    seed = 123L,
    inclusion = c("finite", "screened", "strict"),
    grad_tol = 1e-4,
    singular_tol = 1e-6,
    lrt_tol = 1e-6,
    rel.tol = 1e-7,
    subdivisions = 200L,
    progress = TRUE
) {

    # --------------------------------------------------
    # 0. Validate inputs
    # --------------------------------------------------

    inclusion <- match.arg(inclusion)

    weights <- match.arg(
        weights,
        choices = c("oracle", "plugin"),
        several.ok = TRUE
    )

    weights <- unique(weights)

    required <- c(
        "n_clusters",
        "cluster_size",
        "beta",
        "Gamma",
        "theta"
    )

    if (!is.list(scenario) ||
        !all(required %in% names(scenario))) {
        stop(
            "scenario must contain n_clusters, cluster_size, ",
            "beta, Gamma, and theta."
        )
    }

    if (length(R) != 1L ||
        !is.finite(R) ||
        R < 1 ||
        R != floor(R)) {
        stop("R must be a positive integer.")
    }

    R <- as.integer(R)

    if (length(alpha) != 1L ||
        !is.finite(alpha) ||
        alpha <= 0 ||
        alpha >= 1) {
        stop("alpha must lie strictly between 0 and 1.")
    }

    if (length(seed) != 1L ||
        !is.finite(seed) ||
        seed != floor(seed)) {
        stop("seed must be an integer.")
    }

    if (length(lrt_tol) != 1L ||
        !is.finite(lrt_tol) ||
        lrt_tol < 0) {
        stop("lrt_tol must be nonnegative and finite.")
    }

    if (!is.logical(progress) ||
        length(progress) != 1L ||
        is.na(progress)) {
        stop("progress must be TRUE or FALSE.")
    }

    if (!is.numeric(scenario$beta) ||
        length(scenario$beta) != 2L ||
        any(!is.finite(scenario$beta))) {
        stop("scenario$beta must contain two finite values.")
    }

    if (!is.matrix(scenario$Gamma) ||
        !identical(dim(scenario$Gamma), c(2L, 2L)) ||
        any(!is.finite(scenario$Gamma))) {
        stop("scenario$Gamma must be a finite 2 x 2 matrix.")
    }

    if (!isTRUE(all.equal(
        scenario$Gamma,
        t(scenario$Gamma),
        tolerance = 1e-10
    ))) {
        stop("scenario$Gamma must be symmetric.")
    }

    if (min(eigen(
        scenario$Gamma,
        symmetric = TRUE,
        only.values = TRUE
    )$values) < -1e-10) {
        stop("scenario$Gamma must be positive semidefinite.")
    }

    if (length(scenario$theta) != 1L ||
        !is.finite(scenario$theta) ||
        scenario$theta <= 0) {
        stop("scenario$theta must be positive.")
    }

    if (length(scenario$n_clusters) != 1L ||
        !is.finite(scenario$n_clusters) ||
        scenario$n_clusters < 2 ||
        scenario$n_clusters != floor(scenario$n_clusters)) {
        stop("n_clusters must be an integer >= 2.")
    }

    if (length(scenario$cluster_size) != 1L ||
        !is.finite(scenario$cluster_size) ||
        scenario$cluster_size < 2 ||
        scenario$cluster_size != floor(scenario$cluster_size)) {
        stop("cluster_size must be an integer >= 2.")
    }

    # --------------------------------------------------
    # 1. Preserve the user's RNG state
    # --------------------------------------------------

    had_seed <- exists(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
    )

    if (had_seed) {
        old_seed <- get(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
        )
    }

    on.exit({
        if (had_seed) {
            assign(
                ".Random.seed",
                old_seed,
                envir = .GlobalEnv
            )
        } else if (exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
        )) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    set.seed(seed)

    # --------------------------------------------------
    # 2. Construct balanced design matrices
    # --------------------------------------------------

    n <- as.integer(scenario$n_clusters)
    m <- as.integer(scenario$cluster_size)

    x <- seq(0, 1, length.out = m)

    X_one <- cbind(
        "(Intercept)" = 1,
        x = x
    )

    Z_one <- X_one

    X_list <- replicate(
        n,
        X_one,
        simplify = FALSE
    )

    Z_list <- replicate(
        n,
        Z_one,
        simplify = FALSE
    )

    # --------------------------------------------------
    # 3. Oracle chi-bar-square weights
    # --------------------------------------------------

    oracle_weights <- NULL
    oracle_information <- NULL
    oracle_error <- NULL

    if ("oracle" %in% weights) {

        oracle_result <- tryCatch({

            fim <- glmm_fim_theoretical(
                X = X_list,
                Z = Z_list,
                beta = scenario$beta,
                family = "nbinom",
                dispersion = scenario$theta
            )

            cw <- chibar_weights_theoretical(
                information = fim$I_metric,
                rel.tol = rel.tol,
                subdivisions = subdivisions
            )

            list(
                weights = cw$w,
                information = fim$I_metric
            )

        }, error = function(e) {
            oracle_error <<- conditionMessage(e)
            NULL
        })

        if (!is.null(oracle_result)) {
            oracle_weights <- oracle_result$weights
            oracle_information <- oracle_result$information
        }
    }

    # --------------------------------------------------
    # 4. Helpers
    # --------------------------------------------------

    safe_diagnostic <- function(fit, name, default = NA) {

        if (is.null(fit) ||
            is.null(fit$diagnostics) ||
            is.null(fit$diagnostics[[name]])) {
            return(default)
        }

        fit$diagnostics[[name]]
    }

    valid_weights <- function(w) {

        is.numeric(w) &&
            length(w) == 4L &&
            all(is.finite(w)) &&
            all(w >= -1e-8) &&
            abs(sum(w) - 1) <= 1e-6
    }

    mixture_pvalue <- function(T, w) {

        if (!is.finite(T) || !valid_weights(w)) {
            return(NA_real_)
        }

        if (T <= 0) {
            return(1)
        }

        w <- pmax(w, 0)
        w <- w / sum(w)

        # For T > 0, the chi-square_0 component
        # contributes zero survival probability.
        p <- sum(
            w[2:4] * stats::pchisq(
                T,
                df = 1:3,
                lower.tail = FALSE
            )
        )

        min(1, max(0, p))
    }

    empty_record <- function(r) {

        data.frame(
            replicate = r,
            statistic = NA_real_,
            p_oracle = NA_real_,
            p_plugin = NA_real_,
            reject_oracle = NA,
            reject_plugin = NA,
            usable = FALSE,
            alt_optimizer_ok = NA,
            null_optimizer_ok = NA,
            alt_hessian_ok = NA,
            null_hessian_ok = NA,
            alt_stationary = NA,
            null_stationary = NA,
            alt_numerically_acceptable = NA,
            null_numerically_acceptable = NA,
            alt_gradient = NA_real_,
            null_gradient = NA_real_,
            near_singular = NA,
            min_eigenvalue = NA_real_,
            max_eigenvalue = NA_real_,
            alt_convergence_code = NA_integer_,
            null_convergence_code = NA_integer_,
            failure_reason = NA_character_,
            stringsAsFactors = FALSE
        )
    }

    # --------------------------------------------------
    # 5. Replicate loop
    # --------------------------------------------------

    records <- vector("list", R)

    for (r in seq_len(R)) {

        rec <- empty_record(r)

        tryCatch({

            # ------------------------------------------
            # Simulate observations
            # ------------------------------------------

            dat <- simulate_nbinom_glmm(
                n_clusters = n,
                cluster_size = m,
                beta = scenario$beta,
                Gamma = scenario$Gamma,
                theta = scenario$theta
            )

            # ------------------------------------------
            # Fit alternative and null models
            # ------------------------------------------

            fit <- fit_nbinom_glmm(
                data = dat,
                fit_null = TRUE,
                grad_tol = grad_tol,
                singular_tol = singular_tol
            )

            if (r == bad[1]) {

                alt <- fit$alternative$fit
                null <- fit$null$fit
            
                print(data.frame(
                    model = c("alternative", "null"),
                    reported_logLik = c(
                        as.numeric(logLik(alt)),
                        as.numeric(logLik(null))
                    ),
                    optimizer_objective = c(
                        alt$fit$objective,
                        null$fit$objective
                    ),
                    hessian_ok = c(
                        alt$sdr$pdHess,
                        null$sdr$pdHess
                    ),
                    convergence_code = c(
                        alt$fit$convergence,
                        null$fit$convergence
                    )
                ))
            }

            alt <- fit$alternative
            null <- fit$null

            # ------------------------------------------
            # Record diagnostics regardless of whether
            # the likelihood ratio is usable
            # ------------------------------------------

            rec$alt_optimizer_ok <-
                safe_diagnostic(alt, "optimizer_ok")

            rec$null_optimizer_ok <-
                safe_diagnostic(null, "optimizer_ok")

            rec$alt_hessian_ok <-
                safe_diagnostic(alt, "hessian_ok")

            rec$null_hessian_ok <-
                safe_diagnostic(null, "hessian_ok")

            rec$alt_stationary <-
                safe_diagnostic(alt, "stationary")

            rec$null_stationary <-
                safe_diagnostic(null, "stationary")

            rec$alt_numerically_acceptable <-
                safe_diagnostic(
                    alt,
                    "numerically_acceptable"
                )

            rec$null_numerically_acceptable <-
                safe_diagnostic(
                    null,
                    "numerically_acceptable"
                )

            rec$alt_gradient <-
                safe_diagnostic(
                    alt,
                    "max_abs_gradient",
                    NA_real_
                )

            rec$null_gradient <-
                safe_diagnostic(
                    null,
                    "max_abs_gradient",
                    NA_real_
                )

            rec$near_singular <-
                safe_diagnostic(
                    alt,
                    "near_singular"
                )

            rec$min_eigenvalue <-
                safe_diagnostic(
                    alt,
                    "min_eigenvalue",
                    NA_real_
                )

            rec$max_eigenvalue <-
                safe_diagnostic(
                    alt,
                    "max_eigenvalue",
                    NA_real_
                )

            rec$alt_convergence_code <-
                safe_diagnostic(
                    alt,
                    "convergence_code",
                    NA_integer_
                )

            rec$null_convergence_code <-
                safe_diagnostic(
                    null,
                    "convergence_code",
                    NA_integer_
                )

            # ------------------------------------------
            # Calculate likelihood-ratio statistic
            # ------------------------------------------

            ll_alt <- alt$logLik
            ll_null <- null$logLik

            if (!is.finite(ll_alt) ||
                !is.finite(ll_null)) {
                rec$failure_reason <- "nonfinite_likelihood"
                return_value <- NULL
            } else {

                T_raw <- 2 * (ll_alt - ll_null)

                rec$statistic <- T_raw

                if (!is.finite(T_raw) ||
                    T_raw < -lrt_tol) {

                    rec$failure_reason <- "invalid_lrt"

                } else {

                    # Only clamp a tiny negative value
                    # attributable to numerical error.
                    T <- max(0, T_raw)

                    rec$statistic <- T

                    # ----------------------------------
                    # Inclusion criteria
                    # ----------------------------------

                    acceptable <- switch(
                        inclusion,

                        finite = TRUE,

                        screened = isTRUE(
                            rec$alt_numerically_acceptable
                        ) && isTRUE(
                            rec$null_numerically_acceptable
                        ),

                        strict = isTRUE(
                            rec$alt_optimizer_ok
                        ) && isTRUE(
                            rec$null_optimizer_ok
                        ) && isTRUE(
                            rec$alt_numerically_acceptable
                        ) && isTRUE(
                            rec$null_numerically_acceptable
                        )
                    )

                    if (!acceptable) {

                        rec$failure_reason <-
                            "excluded_by_inclusion_rule"

                    } else {

                        rec$usable <- TRUE

                        # ------------------------------
                        # Oracle p-value
                        # ------------------------------

                        if ("oracle" %in% weights &&
                            valid_weights(oracle_weights)) {

                            rec$p_oracle <- mixture_pvalue(
                                T,
                                oracle_weights
                            )

                            rec$reject_oracle <-
                                rec$p_oracle < alpha
                        }

                        # ------------------------------
                        # Plug-in p-value
                        # ------------------------------

                        if ("plugin" %in% weights) {

                            plugin_result <- tryCatch({

                                fim <- glmm_fim_theoretical(
                                    X = X_list,
                                    Z = Z_list,
                                    beta = null$beta,
                                    family = "nbinom",
                                    dispersion = null$theta
                                )

                                cw <- chibar_weights_theoretical(
                                    information = fim$I_metric,
                                    rel.tol = rel.tol,
                                    subdivisions = subdivisions
                                )

                                cw$w

                            }, error = function(e) {
                                rec$failure_reason <-
                                    paste0(
                                        "plugin_weights: ",
                                        conditionMessage(e)
                                    )
                                NULL
                            })

                            if (valid_weights(plugin_result)) {

                                rec$p_plugin <- mixture_pvalue(
                                    T,
                                    plugin_result
                                )

                                rec$reject_plugin <-
                                    rec$p_plugin < alpha
                            }
                        }
                    }
                }
            }

        }, error = function(e) {

            rec$failure_reason <- paste0(
                "fit_or_simulation: ",
                conditionMessage(e)
            )
        })

        records[[r]] <- rec

        if (progress &&
            (r %% max(1L, floor(R / 10L)) == 0L ||
             r == R)) {
            message(
                sprintf(
                    "Completed %d / %d replicates",
                    r,
                    R
                )
            )
        }
    }

    # --------------------------------------------------
    # 6. Combine replicate records
    # --------------------------------------------------

    results <- do.call(
        rbind,
        records
    )

    # --------------------------------------------------
    # 7. Summarize rejection rates
    # --------------------------------------------------

    summarize_method <- function(method) {

        p_name <- paste0("p_", method)
        reject_name <- paste0("reject_", method)

        valid <- !is.na(results[[p_name]])

        n_valid <- sum(valid)

        if (n_valid == 0L) {
            return(data.frame(
                method = method,
                n_valid = 0L,
                rejection_rate = NA_real_,
                mcse = NA_real_
            ))
        }

        rate <- mean(
            results[[reject_name]][valid]
        )

        data.frame(
            method = method,
            n_valid = n_valid,
            rejection_rate = rate,
            mcse = sqrt(
                rate * (1 - rate) / n_valid
            )
        )
    }

    summary <- do.call(
        rbind,
        lapply(weights, summarize_method)
    )

    # --------------------------------------------------
    # 8. Diagnostic summary
    # --------------------------------------------------

    count_true <- function(x) {
        sum(x %in% TRUE)
    }

    diagnostics <- list(
        n_requested = R,
        n_usable_lrt = count_true(results$usable),
        n_excluded = R - count_true(results$usable),

        n_alt_optimizer_warnings = sum(
            results$alt_optimizer_ok %in% FALSE
        ),

        n_null_optimizer_warnings = sum(
            results$null_optimizer_ok %in% FALSE
        ),

        n_alt_hessian_failures = sum(
            results$alt_hessian_ok %in% FALSE
        ),

        n_null_hessian_failures = sum(
            results$null_hessian_ok %in% FALSE
        ),

        n_alt_gradient_failures = sum(
            results$alt_stationary %in% FALSE
        ),

        n_null_gradient_failures = sum(
            results$null_stationary %in% FALSE
        ),

        n_near_singular = count_true(
            results$near_singular
        ),

        failure_reasons = table(
            results$failure_reason,
            useNA = "ifany"
        )
    )

    # --------------------------------------------------
    # 9. Return
    # --------------------------------------------------

    list(
        scenario = scenario,
        R = R,
        alpha = alpha,
        inclusion = inclusion,
        weights_requested = weights,
        oracle_weights = oracle_weights,
        oracle_information = oracle_information,
        oracle_error = oracle_error,
        results = results,
        summary = summary,
        diagnostics = diagnostics
    )
}
