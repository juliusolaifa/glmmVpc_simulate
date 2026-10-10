set.seed(123)

dat <- simulate_nbinom_glmm(
    n_clusters = 100,
    cluster_size = 5,
    beta = c(0.5, 0.2),
    Gamma = matrix(0, 2, 2),
    theta = 2
)

fit <- fit_nbinom_glmm(dat)

# Construct cluster-specific matrices
X_all <- model.matrix(~ x, data = dat)

idx <- split(seq_len(nrow(dat)), dat$id)

X <- lapply(
    idx,
    function(i) X_all[i, , drop = FALSE]
)

Z <- X

# Run chi-bar-square LRT
test <- run_chibar_test(
    fit = fit,
    X = X,
    Z = Z
)

test$statistic
test$p_value
test$weights
test$integration_error
